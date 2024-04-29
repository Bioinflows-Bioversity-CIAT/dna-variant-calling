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

