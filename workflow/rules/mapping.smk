rule map_reads:
    input:
        reads = get_trimmed_reads,
        idx = rules.bwa_index.output,
    output:
        "results/{plate}/mapping/bwa/{ref}/{sample}.sorted.bam"
    log:
        "results/{plate}/mapping/bwa/{ref}/log/bwa_{sample}.log"
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=get_read_group,
        sorting="samtools",
        sort_order="coordinate",
    resources:
        tmpdir = get_big_temp
    threads: resources['bwa_mem']['threads']
    wrapper:
        "v3.10.2/bio/bwa/mem"
        

rule samtools_index:
    input:
        rules.map_reads.output,
    output:
        "results/{plate}/mapping/bwa/{ref}/{sample}.sorted.bam.bai"
    log:
        "results/{plate}/mapping/bwa/{ref}/log/index_bam_{sample}.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v3.10.2/bio/samtools/index"