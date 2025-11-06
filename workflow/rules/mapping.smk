rule minimap:
    input:
        fastq = get_sample_fastq,
        ref = f"{basedir}/resources/{{ref}}/{{ref}}.fasta"
    output:
        sam = temp(f"{basedir}/results/{{plate}}/mapping/minimap2/{{ref}}/{{sample}}.sam")
    conda:
        "anchorwave"
    shell:
        """
        minimap2 -ax map-ont {input.ref} {input.fastq} > {output.sam}
        """

rule sam_bam:
    input:
        sam = f"{basedir}/results/{{plate}}/mapping/minimap2/{{ref}}/{{sample}}.sam"
    output:
        bam = temp(f"{basedir}/results/{{plate}}/mapping/minimap2/{{ref}}/{{sample}}.bam"),
        bam_sort = temp(f"{basedir}/results/{{plate}}/mapping/minimap2/{{ref}}/{{sample}}_sort.bam"),
        bam_index = temp(f"{basedir}/results/{{plate}}/mapping/minimap2/{{ref}}/{{sample}}_sort.bam.bai")
    conda:
        "anchorwave"
    shell:
        """
        samtools view -S -b {input.sam} > {output.bam} && \
        samtools sort -o {output.bam_sort} {output.bam} && \
        samtools index {output.bam_sort}
        """ 

rule rehead_bam_file:
    input:
        bam = f"{basedir}/results/{{plate}}/mapping/minimap2/{{ref}}/{{sample}}_sort.bam"
    output:
        bam = f"{basedir}/results/{{plate}}/mapping/minimap2/{{ref}}/{{sample}}_rehead.bam",
        index = f"{basedir}/results/{{plate}}/mapping/minimap2/{{ref}}/{{sample}}_rehead.bam.bai"
    conda:
        "anchorwave"
    shell:
        """
        samtools addreplacerg \
        -r "@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tPL:ONT" \
        --output-fmt BAM \
        -o {output.bam} {input.bam} && \
        samtools index {output.bam}
        """
