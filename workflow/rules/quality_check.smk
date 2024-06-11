rule fastqc_raw:
    input:
        "results/{plate}/demultiplexing/stacks/{sample}.{group}.fq.gz"
    output:
        html="results/{plate}/quality_check/{ref}/multiqc_report_per_sample/{sample}/{sample}.{group}.html",
        zip=temp("results/{plate}/quality_check/{ref}/multiqc_report_per_sample/{sample}/{sample}.{group}_fastqc.zip")
    params:
        extra = "--quiet"
    log:
        "results/{plate}/quality_check/{ref}/logs/fatqc/{sample}.{group}_raw.log"
    threads: 1
    resources:
        mem_mb = 1024
    wrapper:
        "v3.11.0/bio/fastqc"

rule fastqc_trimmed:
    input:
        "results/{plate}/trimming/trimmomatic/paired/{sample}_trimmed.{group}.fastq.gz"
    output:
        html="results/{plate}/quality_check/{ref}/multiqc_report_per_sample/{sample}/{sample}_trimmed.{group}.html",
        zip=temp("results/{plate}/quality_check/{ref}/multiqc_report_per_sample/{sample}/{sample}_trimmed.{group}_fastqc.zip")
    params:
        extra = "--quiet"
    log:
        "results/{plate}/quality_check/{ref}/logs/fatqc/{sample}_trimmed.{group}_raw.log"
    threads: 1
    resources:
        mem_mb = 1024
    wrapper:
        "v3.11.0/bio/fastqc"

rule samtools_stats:
    input:
        bam="results/{plate}/mapping/bwa/{ref}/{sample}.sorted.bam"
    output:
        "results/{plate}/quality_check/{ref}/multiqc_report_per_sample/{sample}/{sample}.sorted.bam.stats",
    log:
        "results/{plate}/quality_check/{ref}/logs/samtools/{sample}_samtools.log",
    wrapper:
        "v3.11.0/bio/samtools/stats"

rule qualimap_bamqc:
    input:
        bam="results/{plate}/mapping/bwa/{ref}/{sample}.sorted.bam",
    output:
        directory("results/{plate}/quality_check/{ref}/multiqc_report_per_sample/{sample}/{sample}")
    conda:
        "../envs/NGSEP.yaml"  
    params:
        mem_size="4G",
    log:
        "results/{plate}/quality_check/{ref}/logs/qualimap/{sample}_qualimap.log",
    shell:
        """
        {config[qualimap][path]} bamqc -bam {input.bam} -outdir {output} --java-mem-size={params.mem_size} 2> {log}
        """

rule bcf_stats:
    input:
        "results/{plate}/variant_calling/NGSEP/{ref}/second_variant_call_plate/{sample}_bwa_NGSEP.vcf.gz",
    output:
        "results/{plate}/quality_check/{ref}/multiqc_report_per_sample/{sample}/{sample}.vcf.stats",
    log:
        "results/{plate}/quality_check/{ref}/logs/bcftools/{sample}.bcftools.stats.log",
    params:
        "",
    wrapper:
        "v3.11.0/bio/bcftools/stats"

rule multiqc:
    input:
        fastqc1r1="results/{plate}/quality_check/{ref}/multiqc_report_per_sample/{sample}/{sample}.1_fastqc.zip",
        fastqc1r2="results/{plate}/quality_check/{ref}/multiqc_report_per_sample/{sample}/{sample}.2_fastqc.zip",
        fastqc2r1="results/{plate}/quality_check/{ref}/multiqc_report_per_sample/{sample}/{sample}_trimmed.1_fastqc.zip",
        fastqc2r2="results/{plate}/quality_check/{ref}/multiqc_report_per_sample/{sample}/{sample}_trimmed.2_fastqc.zip",
        samtools="results/{plate}/quality_check/{ref}/multiqc_report_per_sample/{sample}/{sample}.sorted.bam.stats",
        qualimap="results/{plate}/quality_check/{ref}/multiqc_report_per_sample/{sample}/{sample}",
        bcftools="results/{plate}/quality_check/{ref}/multiqc_report_per_sample/{sample}/{sample}.vcf.stats",
        config = "/home/mnarvaez/pipeline/dna-variant-calling/data/multiqc_config.yaml"
    output:
        html="results/{plate}/quality_check/{ref}/multiqc_report_per_sample/{sample}/multiqc_report.html",
        data_dir=directory("results/{plate}/quality_check/{ref}/multiqc_report_per_sample/{sample}/multiqc_data")
    params:
        extra=lambda wildcards: "-i 'Quality report for sample {sample}'".format(sample=wildcards.sample),
        use_input_files_only=True,
    log:
        "results/{plate}/quality_check/{ref}/logs/multiqc/{sample}.multiqc.stats.log"
    threads: 1
    resources:
        mem_mb=1024
    wrapper:
        "v3.11.0/bio/multiqc"

rule merged_vcf_stats:
    input:
        vcf = "results/{plate}/variant_calling/NGSEP/{ref}/{plate}_merged.vcf.gz",
        ref = rules.copy_reference.output,
        genome_dict = rules.create_dict.output,
        index = rules.genome_faidx.output,
        index_vcf = rules.merge_vcfs_by_plate.output.index
    output:
        "results/{plate}/quality_check/{ref}/viariantqc/variantqc.html",
    log:
        "results/{plate}/quality_check/{ref}/viariantqc/variantqc.log",
    conda:
        "../envs/GATK_plus_variantqc.yaml"
    params:
        "",
    shell:
        """
        java -jar {config[DISCVR-Seq][path]} VariantQC -R {input.ref} -V {input.vcf} -O {output}
        """