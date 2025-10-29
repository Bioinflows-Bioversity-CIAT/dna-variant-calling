rule map_reads2:
    input:
        reads = get_trimmed_reads,
        idx = rules.bwa_index.output,
        ref = rules.copy_reference.output,
    output:
        f"{basedir}/results/{{plate}}/mapping/bwa/{{ref}}/{{sample}}.sam"
    log:
        f"{basedir}/log/mapping/{{plate}}/{{ref}}/{{sample}}_bwa_mem.log"
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=get_read_group,
        sorting="samtools",
        sort_order="coordinate",
    resources:
        tmpdir = get_big_temp
    threads: resources['bwa_mem']['threads']
    conda:
        "ngs"
    shell:
        """
        bwamem -t {threads} \
        -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:Illumina' \
        {input.ref} \
        {input.reads} > {output} 2> {log}

        """
rule map_reads:
    input:
        reads = get_trimmed_reads,
        idx = rules.bwa_index.output,
    output:
        f"{basedir}/results/{{plate}}/mapping/bwa/{{ref}}/{{sample}}.sorted.bam"
    log:
        f"{basedir}/log/mapping/{{plate}}/{{ref}}/{{sample}}_bwa_mem.log"
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort="samtools",  # Can be 'none', 'samtools', or 'picard'.
        sort_order="coordinate",  # Can be 'coordinate' (default) or 'queryname'.
        sort_extra="",  # Extra args for samtools/picard sorts.
    threads: resources['bwa_mem']['threads']
    wrapper:
        "v5.9.0/bio/bwa-mem2/mem"        


rule samtools_index:
    input:
        rules.map_reads.output,
    output:
        f"{basedir}/results/{{plate}}/mapping/bwa/{{ref}}/{{sample}}.sorted.bam.bai"
    log:
        f"{basedir}/log/mapping/{{plate}}/{{ref}}/{{sample}}_index_bam.log"
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v4.7.2/bio/samtools/index"
