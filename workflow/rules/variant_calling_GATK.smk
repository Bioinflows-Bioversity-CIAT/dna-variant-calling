# Rule for running GATK HaplotypeCaller
rule haplotype_caller:
    input:
        bam = rules.map_reads.output,  # BAM file
        bai = rules.samtools_index.output,  # BAM index file
        ref = rules.copy_reference.output,  # Reference genome
        genome_dict = rules.create_dict.output  # Genome dictionary
    output:
        gvcf = 'results/{plate}/variant_calling/GATK/{ref}/HaplotypeCaller/intervals/{chrom}/{sample}_{chrom}.g.vcf.gz'  # GVCF output
    log:
        'results/{plate}/variant_calling/GATK/{ref}/log/HaplotypeCaller/{chrom}/{sample}_{chrom}.log'  # Log file
    params:
        extra = get_GATK_HaplotypeCaller_params(),  # Extra parameters
        intervals = lambda wildcards: f"{wildcards.chrom}"  # Intervals (chromosome)
    threads: resources["GATK"]["HaplotypeCaller"]['threads']  # Number of threads
    resources:
        mem_mb = resources["GATK"]["HaplotypeCaller"]['mem']  # Memory allocation
    wrapper:
        "v3.10.2/bio/gatk/haplotypecaller"  # Wrapper

# Rule for running GATK CombineGVCFs
rule combine_gvcfs:
    input:
        gvcfs = get_gvcfs_by_sample,  # List of GVCF files
        ref = rules.copy_reference.output  # Reference genome
    output:
        gvcf = "results/{plate}/variant_calling/GATK/{ref}/CombineGVCFs/{sample}.g.vcf.gz"  # Combined GVCF output
    log:
        'results/{plate}/variant_calling/GATK/{ref}/log/CombineGVCFs/{sample}.log'  # Log file
    params:
        extra = get_GATK_CombineGVCFs_params()  # Extra parameters
    resources:
        mem_mb = resources['GATK']['CombineGVCFs']['mem']  # Memory allocation
    wrapper:
        "v3.10.2/bio/gatk/combinegvcfs"  # Wrapper

# Rule for running GATK GenomicsDBImport
rule genomics_db_import:
    input:
        gvcfs = get_gvcfs_DB  # List of GVCF files
    output:
        db = directory("results/{plate}/variant_calling/GATK/{ref}/DB/{chrom}")  # Genomic database output
    log:
        'results/{plate}/variant_calling/GATK/{ref}/log/GenomicsDBImport/{chrom}.log'  # Log file
    params:
        extra = get_GenomicsDBImport_params(),  # Extra parameters
        intervals = lambda wildcards: "{interval}".format(interval=wildcards.chrom)  # Intervals per chromosome
    threads: 4  # Number of threads
    resources:
        mem_mb = 34000  # Memory allocation
    wrapper:
        "v3.10.2/bio/gatk/genomicsdbimport"  # Wrapper

# Rule for running GATK GenotypeGVCFs
rule genotype_gvcfs:
    input:
        genomicsdb = rules.genomics_db_import.output.db,  # Genomic database
        ref = rules.copy_reference.output  # Reference genome
    output:
        vcf = "results/{plate}/variant_calling/GATK/{ref}/GenotypeGVCFs/{chrom}/{interval_i}-{interval_e}.vcf.gz"  # Genotyped VCF output
    log:
        'results/{plate}/variant_calling/GATK/{ref}/log/GenotypeGVCFs/{chrom}/{interval_i}-{interval_e}.log'  # Log file
    params:
        extra = get_GenotypeGVCFs_params(),  # Extra parameters
        intervals = lambda wildcards: "{chrom}:{interval_i}-{interval_e}".format(
            chrom=wildcards.chrom,
            interval_i=wildcards.interval_i,
            interval_e=wildcards.interval_e
        )  # Intervals (chromosome and range)
    resources:
        mem_mb = 1024  # Memory allocation
    wrapper:
        "v3.10.2/bio/gatk/genotypegvcfs"  # Wrapper

# Rule for concatenating VCF files using bcftools
rule bcftools_concat:
    input:
        calls = get_interval_raw_vcfs,  # List of VCF files
        fai = rules.genome_faidx.output  # Reference index file
    output:
        vcf = "results/{plate}/variant_calling/GATK/{ref}/{plate}.raw.vcf.gz"  # Concatenated VCF output
    log:
        'results/{plate}/variant_calling/GATK/{ref}/log/bcftools_merge/{plate}.log'  # Log file
    params:
        uncompressed_bcf = False,  # Output format
        extra = "-Oz"  # Extra parameters
    threads: 30  # Number of threads
    resources:
        mem_mb = 10240  # Memory allocation
    wrapper:
        "v3.14.0/bio/bcftools/concat"  # Wrapper

# Rule for selecting variants using GATK SelectVariants
rule select_variants:
    input:
        vcf = "results/{plate}/variant_calling/GATK/{ref}/{plate}.raw.vcf.gz",  # Input VCF
        ref = rules.copy_reference.output  # Reference genome
    output:
        biallelic_snp_vcf = "results/{plate}/variant_calling/GATK/{ref}/{plate}.biallelic.snp.vcf.gz"  # Output VCF
    log:
        'results/{plate}/variant_calling/GATK/{ref}/log/SelectVariants/{plate}.log'  # Log file
    params:
        extra = get_SelectVariants_params()  # Extra parameters
    resources:
        mem_mb = 10240  # Memory allocation
    wrapper:
        "v3.14.0/bio/gatk/selectvariants"  # Wrapper

# Rule for filtering variants using GATK VariantFiltration
rule variant_filtering:
    input:
        biallelic_snp_vcf = "results/{plate}/variant_calling/GATK/{ref}/{plate}.biallelic.snp.vcf.gz",  # Input VCF
        ref = rules.copy_reference.output  # Reference genome
    output:
        hardfiltered_vcf = "results/{plate}/variant_calling/GATK/{ref}/{plate}.hardfiltered.snp.vcf.gz"  # Filtered VCF output
    log:
        'results/{plate}/variant_calling/GATK/{ref}/log/VariantFiltration/{plate}.log'  # Log file
    params:
        filters = {"myfilter": "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"}  # Filters
    resources:
        mem_mb = 10240  # Memory allocation
    wrapper:
        "v3.14.0/bio/gatk/variantfiltration"  # Wrapper
