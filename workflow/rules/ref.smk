rule copy_reference:
    input:
        fasta = get_reference_fasta
    output:
        "resources/{ref}/{ref}.fasta"
    shell:
        """
        cp {input.fasta} {output}
        """

checkpoint genome_faidx:
    input:
        path = rules.copy_reference.output
    output:
        "resources/{ref}/{ref}.fasta.fai"
    cache: True
    wrapper:
        "v3.10.2/bio/samtools/faidx"

rule bwa_index:
    input:
        "resources/{ref}/{ref}.fasta",
    output:
         idx=multiext("resources/{ref}.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "logs/{ref}_bwa_index.log"
    params:
        algorithm="bwtsw",
    wrapper:
        "v3.10.2/bio/bwa/index"

checkpoint get_intervals:
    input:
        fai = "resources/{ref}/{ref}.fasta.fai"
    output:
        intervals = "resources/{ref}/{ref}_intervals.txt"
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
        "resources/{ref}/{ref}.fasta"
    output:
        "resources/{ref}/{ref}.dict"
    log:
        "resources/{ref}/{ref}_dict.log",
    params:
    resources:
        mem_mb=1024,
    wrapper:
        "v3.10.2/bio/picard/createsequencedictionary"