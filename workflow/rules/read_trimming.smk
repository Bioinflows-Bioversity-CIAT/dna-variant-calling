rule trim_reads_pe:
    input:
        unpack(get_sample_fastq)
    output:
        r1=f"{base_dir}/results/{{plate}}/trimming/trimmomatic/paired/{{sample}}.1.fastq.gz",
        r2=f"{base_dir}/results/{{plate}}/trimming/trimmomatic/paired/{{sample}}.2.fastq.gz",
        r1_unpaired=temp(f"{base_dir}/results/{{plate}}/trimming/trimmomatic/paired/{{sample}}.1.unpaired.fastq.gz"),
        r2_unpaired=temp(f"{base_dir}/results/{{plate}}/trimming/trimmomatic/paired/{{sample}}.2.unpaired.fastq.gz"),
        trimlog=temp(f"{base_dir}/results/{{plate}}/trimming/trimmomatic/paired/{{sample}}.trimlog.txt")
    params:
        **config["trimmomatic"]["pe"],
        extra=lambda w, output: "-trimlog {output}".format(output = output.trimlog)
    log:
        f"{base_dir}/log/trimming/{{plate}}/trimmomatic/paired/{{sample}}.log"
    threads:
        resources['trimmomatic']['threads']
    resources:
        mem_mb=resources['trimmomatic']['mem']
    wrapper:
        "v4.7.2/bio/trimmomatic/pe"
        
rule trim_reads_se:
    input:
        get_sample_fastq,
    output:
        f"{base_dir}/results/{{plate}}/trimming/trimmomatic/single/{{sample}}.fastq.gz"
    params:
        **config["trimmomatic"]["pe"]
    log:
        f"{base_dir}/log/trimming/{{plate}}/trimmomatic/single/{{sample}}.log"
    threads:
        resources['trimmomatic']['threads']
    resources:
        mem_mb=resources['trimmomatic']['mem']
    wrapper:
        "v4.7.2/bio/trimmomatic/se"