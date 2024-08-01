rule haplotype_caller:
    input:
        # Single or list of bam files
        bam = rules.map_reads.output,
        # Index file for BAM
        bai = rules.samtools_index.output,
        # Reference genome file
        ref = rules.copy_reference.output,
        # Genome dictionary file
        genome_dict = rules.create_dict.output
    output:
        # Output GVCF file for the specified sample and chromosome
        gvcf = 'results/{plate}/variant_calling/GATK/{ref}/HaplotyeCaller/
        intervals/{chrom}/{sample}_{chrom}.g.vcf.gz'
    log:
        # Log file for HaplotypeCaller process
        'results/{plate}/variant_calling/GATK/{ref}/log/HaplotyeCaller/{chrom}/{sample}_{chrom}.log'
    params:
        # Extra parameters for GATK HaplotypeCaller
        extra = get_GATK_HaplotypeCaller_params(),
        # Interval for HaplotypeCaller, typically a chromosome
        intervals = lambda wildcards: f"{wildcards.chrom}"
    threads: resources["GATK"]["HaplotypeCaller"]['threads']
    resources:
        # Number of threads to use, defined in the resources dictionary
        mem_mb = resources["GATK"]["HaplotypeCaller"]['mem']
    wrapper:
        # Memory allocation for HaplotypeCaller
        "v3.10.2/bio/gatk/haplotypecaller"

rule combine_gvcfs:
    input:
        # List of GVCF files for a given sample
        gvcfs = get_gvcfs_by_sample,
        # Reference genome file
        ref = rules.copy_reference.output
    output:
        # Output combined GVCF file for the specified sample
        gvcf ="results/{plate}/variant_calling/GATK/{ref}/CombineGVCFs/{sample}.g.vcf.gz"
    log:
        # Log file for CombineGVCFs process
        'results/{plate}/variant_calling/GATK/{ref}/log/CombineGVCFs/{sample}.log'
    params:
        # Extra parameters for GATK CombineGVCFs
        extra = get_GATK_CombineGVCFs_params() # Optional
    resources:
        # Memory allocation for CombineGVCFs
        mem_mb = resources['GATK']['CombineGVCFs']['mem']
    wrapper:
        # Wrapper for GATK CombineGVCFs
        "v3.10.2/bio/gatk/combinegvcfs"

rule genomics_db_import:
    input:
        # List of GVCF files for the genomic database import
        gvcfs = get_gvcfs_DB
    output:
        # Directory for the genomic database
        db = directory("results/{plate}/variant_calling/GATK/{ref}/DB/{chrom}")
    log:
        # Log file for GenomicsDBImport process
        'results/{plate}/variant_calling/GATK/{ref}/log/GenomicsDBImport/{chrom}.log'
    params:
        # Extra parameters for GATK GenomicsDBImport
        extra = get_GenomicsDBImport_params(),  # Optional
        # Interval for GenomicsDBImport, typically a chromosome
        intervals = lambda wildcards: "{interval}".format(interval = wildcards.chrom)
    threads: 4
    resources:
        # Memory allocation for GenomicsDBImport
        mem_mb = 34000
    wrapper:
        # Wrapper for GATK GenomicsDBImport
        "v3.10.2/bio/gatk/genomicsdbimport"

rule genotype_gvcfs:
    input:
        # Genomic database directory from GenomicsDBImport
        genomicsdb = rules.genomics_db_import.output.db,
        # Reference genome file
        ref = rules.copy_reference.output
    output:
        # Output VCF file with genotyped variants for the specified intervals
        vcf = "results/{plate}/variant_calling/GATK/{ref}/GenotypeGVCFs/{chrom}/{interval_i}-{interval_e}.vcf.gz"
    log:
        # Log file for GenotypeGVCFs process
        'results/{plate}/variant_calling/GATK/{ref}/log/GenotypeGVCFs/{chrom}/{interval_i}-{interval_e}.log'
    params:
        # Extra parameters for GATK GenotypeGVCFs
        extra = get_GenotypeGVCFs_params(), # Optional
        # Interval for GenotypeGVCFs, specified as chromosome and range
        intervals = lambda wildcards: "{chrom}:{interval_i}-{interval_e}".format(										chrom = wildcards.chrom, 																			interval_i = wildcards.interval_i,																	interval_e = wildcards.interval_e
        )
    resources:
        # Memory allocation for GenotypeGVCFs
        mem_mb = 1024
    wrapper:
        # Wrapper for GATK GenotypeGVCFs
        "v3.10.2/bio/gatk/genotypegvcfs"

rule bcftools_concat:
    input:
    	# A function or list of VCF files to be concatenated
        calls = get_interval_raw_vcfs,
        # The reference index file
        fai = rules.genome_faidx.output 
    output:
    	# The output concatenated VCF file
        vcf = "results/{plate}/variant_calling/GATK/{ref}/{plate}.raw.vcf.gz"
    log:
    	# Log file for the concatenation process
        'results/{plate}/variant_calling/GATK/{ref}/log/bcftools_merge/{plate}.log'
    params:
        # Parameter to control if the output should be uncompressed BCF format
        uncompressed_bcf = False,
        # Optional parameters for bcftools concat (excluding -o, which is handled by Snakemake)
        extra = "-Oz"
    threads: 30
    resources:
        # Memory allocation for the rule
        mem_mb = 10240
    wrapper:
        # The wrapper for bcftools concat
        "v3.10.2/bio/bcftools/concat"

rule select_variants:
    input:
        # Input raw VCF file from bcftools concat
        vcf = "results/{plate}/variant_calling/GATK/{ref}/{plate}.raw.vcf.gz"",
        # Reference genome file
        ref = rules.copy_reference.output
    output:
        # Output VCF file with selected biallelic variants
        biallelic_vcf = "results/{plate}/variant_calling/GATK/{ref}/{plate}.biallelic.vcf.gz"
    log:
        # Log file for the SelectVariants process
        'results/{plate}/variant_calling/GATK/{ref}/log/SelectVariants/{plate}.log'
    params:
        # Optional extra parameters for GATK SelectVariants
        extra = get_GATK_SelectVariants_params() # Optional
    resources:
        # Memory allocation for the rule
        mem_mb = 1024,
    wrapper:
        # The wrapper for GATK SelectVariants
        "v3.14.0/bio/gatk/selectvariants"

rule variant_filtering:
    input:
        # Input VCF file from SelectVariants
        biallelic_vcf = "results/{plate}/variant_calling/GATK/{ref}/{plate}.biallelic.vcf.gz",
        # Reference genome file
        ref = rules.copy_reference.output
    output:
        # Output hard-filtered VCF file
        hardfiltered_vcf = "results/{plate}/variant_calling/GATK/{ref}/{plate}.hardfiltered.vcf.gz"
    log:
        # Log file for the VariantFiltration process
        'results/{plate}/variant_calling/GATK/{ref}/log/VariantFiltration/{plate}.log'
    params:
        # Filters for VariantFiltration (could be expanded for clarity or flexibility)
        filters = {"myfilter": "QD = 2.0 || MQ = 40.0 || FS = 60.0 || SOR = 3.0 || MQR = -12.5 || RPR = -8.0"},
        # Optional extra parameters for GATK VariantFiltration
        extra = get_GATK_VariantFiltration_params() # Optional
    resources:
        # Memory allocation for the rule
        mem_mb = 1024
    wrapper:
        # The wrapper for GATK VariantFiltration
        "v3.14.0/bio/gatk/variantfiltration"