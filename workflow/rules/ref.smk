rule copy_reference:
    input:
        fasta = get_reference_fasta
    output:
        f"{basedir}/resources/{{ref}}/{{ref}}.fasta"
    shell:
        """
        cp {input.fasta} {output}
        """

checkpoint genome_faidx:
    input:
        path = rules.copy_reference.output
    output:
        f"{basedir}/resources/{{ref}}/{{ref}}.fasta.fai"
    cache: True
    wrapper:
        "v4.7.2/bio/samtools/faidx"
    

rule bwa_index:
    input:
        rules.copy_reference.output
    output:
         multiext(f"{basedir}/resources/{{ref}}/{{ref}}.fasta", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")
    log:
        f"{basedir}/log/reference/{{ref}}_bwa_index.log"
    params:
        algorithm="is",
    wrapper:
        "v5.9.0/bio/bwa-mem2/index"

checkpoint get_intervals:
    input:
        fai =  f"{basedir}/resources/{{ref}}/{{ref}}.fasta.fai"
    output:
        intervals = f"{basedir}/resources/{{ref}}/{{ref}}_intervals.txt"
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
        f"{basedir}/resources/{{ref}}/{{ref}}.dict"
    log:
        f"{basedir}/log/reference/{{ref}}_dict.log",
    resources:
        mem_mb=1024,
    wrapper:
        "v4.7.2/bio/picard/createsequencedictionary"
