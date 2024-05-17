rule haplotype_caller:
    input:
        ref = "resources/{ref}.fasta",
        bam = "results/{plate}/mapping/{ref}/bwa/mapping/{sample}.sorted.bam",
        interval = ""
    output:
    log:
    threads:
    params:
    conda:
    shell:
        """
        
        """