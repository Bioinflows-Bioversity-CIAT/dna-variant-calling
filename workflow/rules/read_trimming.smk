rule trim_reads_pe:
    input:
        unpack(get_sample_fastq)
    output:
        r1="results/{plate}/trimming/trimmomatic/paired/{sample}_trimmed.1.fastq.gz",
        r2="results/{plate}/trimming/trimmomatic/paired/{sample}_trimmed.2.fastq.gz",
        r1_unpaired=temp("results/{plate}/trimming/trimmomatic/paired/{sample}.1.unpaired.fastq.gz"),
        r2_unpaired=temp("results/{plate}/trimming/trimmomatic/paired/{sample}.2.unpaired.fastq.gz"),
        trimlog="results/{plate}/trimming/trimmomatic/paired/{sample}.trimlog.txt"
    params:
        **config["trimmomatic"]["pe"],
        extra=lambda w, output: "-trimlog {output}".format(output = output.trimlog)
    log:
        "results/{plate}/trimming/trimmomatic/log/paired/{sample}.log"
    threads:
        resources['trimmomatic']['threads']
    resources:
        mem_mb=resources['trimmomatic']['mem']
    wrapper:
        "v3.10.2/bio/trimmomatic/pe"
        
rule trim_reads_se:
    input:
        get_sample_fastq,
    output:
        "results/{plate}/trimming/trimmomatic/single/{sample}_trimmed.fastq.gz"
    params:
        **config["trimmomatic"]["pe"]
    log:
        "results/{plate}/trimming/trimmomatic/log/single/{sample}_trimmed.log"
    threads:
        resources['trimmomatic']['threads']
    resources:
        mem_mb=resources['trimmomatic']['mem']
    wrapper:
        "v3.10.2/bio/trimmomatic/se"