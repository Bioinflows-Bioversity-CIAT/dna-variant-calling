# Rule for mapping reads using BWA
rule map_reads:
    input:
        reads = get_trimmed_reads,  # Trimmed reads
        idx = rules.bwa_index.output,  # BWA index
    output:
        "results/{plate}/mapping/bwa/{ref}/{sample}.sorted.bam"  # Sorted BAM output
    log:
        "results/{plate}/mapping/bwa/{ref}/log/bwa_{sample}.log"  # Log file
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],  # BWA index base
        extra=get_read_group,  # Extra parameters
        sorting="samtools",  # Sorting method
        sort_order="coordinate",  # Sort order
    resources:
        tmpdir = get_big_temp  # Temporary directory
    threads: resources['bwa_mem']['threads']  # Number of threads
    wrapper:
        "v3.10.2/bio/bwa/mem"  # Wrapper

# Rule for indexing BAM files using samtools
rule samtools_index:
    input:
        rules.map_reads.output,  # BAM file
    output:
        "results/{plate}/mapping/bwa/{ref}/{sample}.sorted.bam.bai"  # BAM index output
    log:
        "results/{plate}/mapping/bwa/{ref}/log/index_bam_{sample}.log",  # Log file
    params:
        extra="",  # Optional parameters
    threads: 4  # Number of threads
    wrapper:
        "v3.10.2/bio/samtools/index"  # Wrapper
