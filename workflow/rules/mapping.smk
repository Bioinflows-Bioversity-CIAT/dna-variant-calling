rule map_reads:
    input:
        reads = get_trimmed_reads,
        idx = rules.bwa_index.output,
    output:
        "results/{plate}/mapping/{ref}/bwa/mapping/{sample}.sorted.bam"
    log:
        "results/{plate}/logs/bwa/{ref}/{sample}.log"
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=get_read_group,
        sorting="samtools",
        sort_order="coordinate",
    resources:
        tmpdir = get_big_temp
    threads: resources['bwa_mem']['threads']
    wrapper:
        "file:///home/scruz/software/snakemake-wrappers/bio/bwa/mem"              

