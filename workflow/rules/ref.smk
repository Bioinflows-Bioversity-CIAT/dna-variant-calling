rule copy_reference:
    input:
        fasta = get_reference_fasta
    output:
        "resources/{ref}.fasta"
    shell:
        """
        cp {input.fasta} {output}
        """

rule genome_faidx:
    input:
        path = "resources/{ref}.fasta"
    output:
        "resources/{ref}.fasta.fai"
    cache: True
    wrapper:
        "file:///home/scruz/software/snakemake-wrappers/bio/samtools/faidx"

rule bwa_index:
    input:
        "resources/{ref}.fasta",
    output:
         idx=multiext("resources/{ref}.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "logs/{ref}_bwa_index.log"
    params:
        algorithm="bwtsw",
    wrapper:
        "file:///home/scruz/software/snakemake-wrappers/bio/bwa/index"

rule get_intervals:
    input:
        fai = "resources/{ref}.fasta.fai"
    output:
        intervals = "resources/{ref}_intervals.txt"
    params:
        l = config['GATK']['interval_length']
    run:
        with open(output.intervals, "w") as out:
            with open(input.fai, "r") as f:
                for line in f:
                    line = line.strip().split('\t')
                    chrom = line[0]
                    length = int(line[1])
                    start = 1
                    end = params.l
                    while (start <= length):
                        if(end > length):
                            end = length
                        interval_text = chrom +":"+str(start)+"-"+str(end)
                        print(interval_text, file=out)
                        start += params.l
                        end += params.l
rule create_dict:
    input:
        "resources/{ref}.fasta"
    output:
        "resources/{ref}.dict"
    log:
        "resources/{ref}_dict.log",
    params:
    resources:
        mem_mb=1024,
    wrapper:
        "file:///home/scruz/software/snakemake-wrappers/bio/picard/createsequencedictionary"