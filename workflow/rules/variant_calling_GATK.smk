
rule haplotype_caller_gvcf:
    input:
        # single or list of bam files
        bam = "results/{plate}/mapping/{ref}/bwa/mapping/{sample}.sorted.bam",
        bai = "results/{plate}/mapping/{ref}/bwa/mapping/{sample}.sorted.bam.bai",
        ref = "resources/{ref}.fasta",
        genome_dict = "resources/{ref}.dict",
    output:
        gvcf='results/{plate}/mapping/{ref}/GATK/gvcf/intervals/{sample}/{sample}_{interval}.g.vcf.gz'
    log:
        'results/{plate}/mapping/{ref}/GATK/gvcf/intervals/{sample}/{sample}_{interval}.log'
    params:
        extra=get_GATK_HaplotypeCaller_params(),
        intervals = lambda wildcards: f"{wildcards.interval}"

    threads: resources["GATK"]["HaplotypeCaller"]['threads']
    resources:
        mem_mb=resources["GATK"]["HaplotypeCaller"]['mem'],
    wrapper:
        "file:///home/scruz/software/snakemake-wrappers/bio/gatk/haplotypecaller"

rule combine_by_sample_gvcfs:
    input:
        gvcfs = get_gvcfs_by_sample,
        ref = "resources/{ref}.fasta",
    output:
        gvcf ="results/{plate}/mapping/{ref}/GATK/gvcf/{sample}.g.vcf.gz",
    log:
        "results/{plate}/mapping/{ref}/GATK/gvcf/{sample}_comb.log",
    params:
        extra = get_GATK_CombineGVCFs_params(),  
    resources:
        mem_mb=resources['GATK']['CombineGVCFs']['mem'],
    wrapper:
        "file:///home/scruz/software/snakemake-wrappers/bio/gatk/combinegvcfs"




