rule haplotype_caller_gvcf:
    input:
        # single or list of bam files
        bam = rules.map_reads.output,
        bai = rules.samtools_index.output,
        ref = rules.copy_reference.output,
        genome_dict = rules.create_dict.output
    output:
        gvcf='results/{plate}/variant_calling/GATK/{ref}/HaplotyeCaller/intervals/{chrom}/{sample}_{chrom}.g.vcf.gz'
    log:
        'results/{plate}/variant_calling/GATK/{ref}/log/HaplotyeCaller/{chrom}/{sample}_{chrom}.log'
    params:
        extra=get_GATK_HaplotypeCaller_params(),
        intervals = lambda wildcards: f"{wildcards.chrom}"

    threads: resources["GATK"]["HaplotypeCaller"]['threads']
    resources:
        mem_mb=resources["GATK"]["HaplotypeCaller"]['mem'],
    wrapper:
        "v3.10.2/bio/gatk/haplotypecaller"

rule combine_by_sample_gvcfs:
    input:
        gvcfs = get_gvcfs_by_sample,
        ref = rules.copy_reference.output,
    output:
        gvcf ="results/{plate}/variant_calling/GATK/{ref}/CombineGVCFs/{sample}.g.vcf.gz",
    log:
        'results/{plate}/variant_calling/GATK/{ref}/log/CombineGVCFs/{sample}.log'
    params:
        extra = get_GATK_CombineGVCFs_params(),  
    resources:
        mem_mb=resources['GATK']['CombineGVCFs']['mem'],
    wrapper:
        "v3.10.2/bio/gatk/combinegvcfs"


rule genomics_db_import:
    input:
        gvcfs=get_gvcfs_DB,
    output:
        db=directory("results/{plate}/variant_calling/GATK/{ref}/DB/{chrom}"),
    log:
        'results/{plate}/variant_calling/GATK/{ref}/log/GenomicsDBImport/{chrom}.log'
    params:
        extra= get_GenomicsDBImport_params(),  # optional
        intervals = lambda wildcards: "{interval}".format(interval = wildcards.chrom)
    threads: 4
    resources:
        mem_mb=34000,
    wrapper:
        "v3.10.2/bio/gatk/genomicsdbimport"


rule genotype_gvcfs:
    input:
        genomicsdb = rules.genomics_db_import.output.db,
        ref=rules.copy_reference.output,
    output:
        vcf = "results/{plate}/variant_calling/GATK/{ref}/GenotypeGVCFs/{chrom}/{interval_i}-{interval_e}.vcf.gz"
    log:
        'results/{plate}/variant_calling/GATK/{ref}/log/GenotypeGVCFs/{chrom}/{interval_i}-{interval_e}.log'
    params:
        extra=get_GenotypeGVCFs_params(),
        intervals = lambda wildcards: "{chrom}:{interval_i}-{interval_e}".format(chrom = wildcards.chrom,
            interval_i = wildcards.interval_i,
            interval_e = wildcards.interval_e)
    resources:
        mem_mb=1024
    wrapper:
        "v3.10.2/bio/gatk/genotypegvcfs"


rule bcftools_concat:
    input:
        calls = get_interval_raw_vcfs,
        fai = rules.genome_faidx.output
    output:
        vcf = "results/{plate}/variant_calling/GATK/{ref}/{plate}.raw.vcf.gz"
    log:
        'results/{plate}/variant_calling/GATK/{ref}/log/bcftools_merge/{plate}.log'
    params:
        uncompressed_bcf=False,
        extra="-Oz",  # optional parameters for bcftools concat (except -o)
    threads: 30
    resources:
        mem_mb=10240,
    wrapper:
        "v3.10.2/bio/bcftools/concat"