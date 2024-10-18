rule copy_reference:
    input:
        fasta = get_reference_fasta
    output:
        f"{base_dir}/resources/{{ref}}/{{ref}}.fasta"
    shell:
        """
        cp {input.fasta} {output}
        """

rule genome_faidx:
    input:
        path = rules.copy_reference.output
    output:
        f"{base_dir}/resources/{{ref}}/{{ref}}.fasta.fai"
    cache: True
    conda:
        "../envs/ngs.yaml"
    shell:
        """
        samtools faidx {input.path}
        """

rule bwa_index:
    input:
        rules.copy_reference.output
    output:
         idx=multiext("resources/{ref}.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        f"{base_dir}/logs/{{ref}}_bwa_index.log"
    params:
        algorithm="is",
    wrapper:
        "v4.7.2/bio/bwa/index"

rule get_intervals:
    input:
        fai = rule.genome_faidx.output
    output:
        intervals = f"{base_dir}/resources/{{ref}}/{{ref}}_intervals.txt"
    params:
        l = config['GATK']['interval_length'] # NOT USED
    run:
        with open(output.intervals, "w") as out:
            with open(input.fai, "r") as f:
                for line in f:
                    line = line.strip().split('\t')
                    chrom = line[0]
                    print(chrom, file=out)
rule create_dict:
    input:
        rules.copy_reference.output
    output:
        f"{base_dir}/resources/{{ref}}/{{ref}}.dict"
    log:
        f"{base_dir}/resources/{{ref}}/{{ref}}_dict.log",
    resources:
        mem_mb=1024,
    wrapper:
        "v4.7.2/bio/picard/createsequencedictionary"