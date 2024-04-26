rule genome_faidx:
    input:
        path = "resources/{ref}.fasta"
    output:
        "resources/{ref}.fasta.fai"
    cache: True
    wrapper:
        "v3.9.0/bio/samtools/faidx"

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
        "v3.9.0/bio/bwa/index"