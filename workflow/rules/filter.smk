rule first_filtering:
    input:
         f"{basedir}/results/{{plate}}/variant_calling/NGSEP/{{ref}}/{{plate}}_merged.vcf.gz"
    output:
         temp(f"{basedir}/results/{{plate}}/variant_calling/NGSEP/{{ref}}/{{plate}}_merged_Q{{qual}}.vcf.gz")
    log:
        f'{basedir}/results/{{plate}}/variant_calling/NGSEP/{{ref}}/{{plate}}_merged_Q{{qual}}.log'
    params:
        mem = "-Xmx40g"
    conda:
        "../envs/NGSEP.yaml"
    shell:
        """
        java {params.mem} -jar {config[NGSEP][path]} \
            VCFFilter -i {input} -q {wildcards.qual} 2> {log}| \
            bgzip > {output}
        """

rule sort_plate_vcfs:
    input:
        vcf = f"{basedir}/results/{{plate}}/variant_calling/NGSEP/{{ref}}/{{plate}}_merged_Q{{qual}}.vcf.gz"
    output:
        vcf = f"{basedir}/results/{{plate}}/variant_calling/NGSEP/{{ref}}/{{plate}}_merged_Q{{qual}}_sorted.vcf.gz",
        index = f"{basedir}/results/{{plate}}/variant_calling/NGSEP/{{ref}}/{{plate}}_merged_Q{{qual}}_sorted.vcf.gz.tbi"
    conda:
        "ngs"
    resources:
        tmpdir = get_big_temp
    shell:
        """
        if [[ ! -f "{input.vcf}.tbi" ]]
        then
            tabix -p vcf {input.vcf}
        fi
        bcftools sort -T {resources.tmpdir} -Oz -o {output.vcf} {input.vcf} && \
        tabix -p vcf {output.vcf}
        """
        
rule snps_plate_vcfs:
    input:
        vcf = f"{basedir}/results/{{plate}}/variant_calling/NGSEP/{{ref}}/{{plate}}_merged_Q{{qual}}_sorted.vcf.gz",
    output:
        snps_vcf = f"{basedir}/results/{{plate}}/variant_calling/NGSEP/{{ref}}/{{plate}}_merged_Q{{qual}}_Dp{{dp}}_MAF{{maf}}_sorted.vcf.gz"
    conda:
        "ngs"
    resources:
        tmpdir = get_big_temp
    threads: 30
    shell:
        """
        bcftools view -v snps -m 2 -M 2 {input.vcf} | \
        grep -v scaffold | \
        bcftools +fill-tags | \
        bcftools filter -S . -i '(FORMAT/GQ)>={wildcards.qual} & (FORMAT/DP)>={wildcards.dp} & MAF >={wildcards.maf}' | \
        bgzip > {output.snps_vcf}
        """