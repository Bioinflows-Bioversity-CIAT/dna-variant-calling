rule single_sample_variant_detector:
    input:
        bam="results/{plate}/mapping/{ref}/bwa/mapping/{sample}.sorted.bam",
        ref="resources/{ref}.fasta"
    output:
        "results/{plate}/mapping/{ref}/NGSEP/first_variant_calling/{sample}_bwa_NGSEP.vcf.gz"
    params:
        params = " ".join(["-{param} {value}".format(param=i_param, value = config['NGSEP']['SingleSampleVariantsDetector'][i_param]) for i_param in config['NGSEP']['SingleSampleVariantsDetector'].keys()]),
        mem = "-Xmx3g"
    resources:
         mem_mb=3000
    conda:
        "../envs/NGSEP.yaml"
    log:
        "results/{plate}/mapping/{ref}/NGSEP/first_variant_calling/{sample}_bwa_NGSEP.log"
    shell:
        """
        java {params.mem} -jar {config[NGSEP][path]} \
            SingleSampleVariantsDetector {params.params} \
            -r {input.ref} \
            -i {input.bam} \
            -o {output} 2> {log} &&
            cat {output}.vcf | bgzip > {output} &&
            rm {output}.vcf
        """

rule merge_variants_by_plate:
    input:
        vcfs = get_sample_vcfs_by_plate_merge_variants,
        ref_list = "resources/{ref}.fasta.fai"
    params:
        mem = "-Xmx40g"
    resources:
         mem_mb=40000
    conda:
        "../envs/NGSEP.yaml"
    output:
        'results/{plate}/mapping/{ref}/NGSEP/vcf/{plate}_merged_variants.vcf'
    log:
        'results/{plate}/mapping/{ref}/NGSEP/vcf/{plate}_merged_variants.log'
    shell:
        """
        java {params.mem} -jar {config[NGSEP][path]} \
            MergeVariants -s {input.ref_list} -o {output} {input.vcfs} 2> {log}
        """
        
rule single_sample_variant_detector_two:
    input:
        bam="results/{plate}/mapping/{ref}/bwa/mapping/{sample}.sorted.bam",
        ref="resources/{ref}.fasta",
        known_variants = 'results/{plate}/mapping/{ref}/NGSEP/vcf/{plate}_merged_variants.vcf'
    output:
        "results/{plate}/mapping/{ref}/NGSEP/vcf/second_variant_call_plate/{sample}_bwa_NGSEP.vcf.gz"
        
    params:
        params = " ".join(["-{param} {value}".format(param=i_param, value = config['NGSEP']['SingleSampleVariantsDetector'][i_param]) for i_param in config['NGSEP']['SingleSampleVariantsDetector'].keys()]),
        mem = "-Xmx3g",
    conda:
        "../envs/NGSEP.yaml"
    resources:
         mem_mb=3000
    log:
        "results/{plate}/mapping/{ref}/NGSEP/vcf/second_variant_call_plate/{sample}_bwa_NGSEP.log"
    shell:
        """
        java {params.mem} -jar {config[NGSEP][path]} \
            SingleSampleVariantsDetector {params.params} -knownVariants {input.known_variants}\
            -sampleId {wildcards.sample} \
            -r {input.ref} \
            -i {input.bam} \
            -o {output} 2> {log} &&
            cat {output}.vcf | bgzip > {output} &&
            rm {output}.vcf
        """
        
rule merge_vcfs_by_plate:
    input:
        vcfs = get_sample_vcfs_by_plate_merge_vcfs,
        ref_list = "resources/{ref}.fasta.fai"
    params:
         mem = "-Xmx40g"
    resources:
         mem_mb=40000
    conda:
        "../envs/NGSEP.yaml"
    output:
        'results/{plate}/mapping/{ref}/NGSEP/vcf/{plate}_merged.vcf.gz'
    log:
        'results/{plate}/mapping/{ref}/NGSEP/vcf/{plate}_merged.log'

    shell:
        """
        java {params.mem} -jar {config[NGSEP][path]} \
            VCFMerge  -s {input.ref_list} -o {output}.vcf {input.vcfs} 2> {log} && \
            cat {output}.vcf | bgzip > {output} &&
            rm {output}.vcf
        """