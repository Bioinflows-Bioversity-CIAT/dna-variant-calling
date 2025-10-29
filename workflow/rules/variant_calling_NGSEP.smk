rule single_sample_variant_detector:
    input:
        bam=rules.map_reads.output,
        ref=rules.copy_reference.output
    output:
        temp(f"{basedir}/results/{{plate}}/variant_calling/NGSEP/{{ref}}/first_variant_calling/{{sample}}_bwa_NGSEP.vcf.gz")
    params:
        params = " ".join(["-{param} {value}".format(param=i_param, value = config['NGSEP']['SingleSampleVariantsDetector'][i_param]) for i_param in config['NGSEP']['SingleSampleVariantsDetector'].keys()]),
        mem = "-Xmx3g"
    resources:
         mem_mb=3000
    conda:
        "../envs/NGSEP.yaml"
    log:
        f"{basedir}/log/variant_calling/NGSEP/first_variant_calling/{{plate}}/{{ref}}/{{sample}}_NGSEP.log"
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
        ref_list =  f"{basedir}/resources/{{ref}}/{{ref}}.fasta.fai"
    params:
        mem = "-Xmx40g"
    resources:
         mem_mb=40000
    conda:
        "../envs/NGSEP.yaml"
    output:
        temp(f'{basedir}/results/{{plate}}/variant_calling/NGSEP/{{ref}}/{{plate}}_merged_variants.vcf')
    log:
        f"{basedir}/log/variant_calling/NGSEP/merge_variants/{{plate}}/{{ref}}/merge_variants_NGSEP.log"
    shell:
        """
        java {params.mem} -jar {config[NGSEP][path]} \
            MergeVariants -s {input.ref_list} -o {output} {input.vcfs} 2> {log}
        """
        
rule single_sample_variant_detector_two:
    input:
        bam=rules.map_reads.output,
        ref=rules.copy_reference.output,
        known_variants = rules.merge_variants_by_plate.output
    output:
        f"{basedir}/results/{{plate}}/variant_calling/NGSEP/{{ref}}/second_variant_call_plate/{{sample}}_bwa_NGSEP.vcf.gz"
    params:
        params = " ".join(["-{param} {value}".format(param=i_param, value = config['NGSEP']['SingleSampleVariantsDetector'][i_param]) for i_param in config['NGSEP']['SingleSampleVariantsDetector'].keys()]),
        mem = "-Xmx3g",
    conda:
        "../envs/NGSEP.yaml"
    resources:
         mem_mb=3000
    log:
        f"{basedir}/log/variant_calling/NGSEP/second_variant_call_plate/{{plate}}/{{ref}}/{{sample}}_NGSEP.log"
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
        ref_list =  f"{basedir}/resources/{{ref}}/{{ref}}.fasta.fai"
    params:
         mem = "-Xmx40g"
    resources:
         mem_mb=40000
    conda:
        "../envs/NGSEP.yaml"
    output:
        f"{basedir}/results/{{plate}}/variant_calling/NGSEP/{{ref}}/{{plate}}_merged.vcf.gz"
    log:
        f"{basedir}/log/variant_calling/NGSEP/merge_vcfs_by_plate/{{plate}}/{{ref}}/merge_vcfs_NGSEP.log"
    shell:
        """
        java {params.mem} -jar {config[NGSEP][path]} \
            VCFMerge  -s {input.ref_list} -o {output}.vcf {input.vcfs} 2> {log} && \
            cat {output}.vcf | bgzip > {output} &&
            rm {output}.vcf
        """