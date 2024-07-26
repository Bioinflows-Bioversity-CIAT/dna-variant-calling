rule fastqc:
    input:
        unpack(get_sample_fastq),
    output:
        html=f"{basedir}/results/{{plate}}/qc/fastqc/{{sample}}.html",
        zip=f"{basedir}/results/{{plate}}/qc/fastqc/{{sample}}.zip",
    params: "--quiet"
    threads: resources['fastqc']['sample']
    wrapper:
        "file:///home/scruz/software/snakemake-wrappers/bio/fastqc"
        
rule fastqc_library:
    input:
        get_library_fastqc
    output:
        html=f"{basedir}/results/{{plate}}/qc/fastqc_run/{{sq_unit}}-{{group}}.html",
        zip=f"{basedir}/results/{{plate}}/qc/fastqc_run/{{sq_unit}}-{{group}}.zip",
    params: "--quiet"
    log:
        f"{basedir}/results/{{plate}}/logs/fastqc/{{sq_unit}}_{{group}}.log",
    threads: resources['fastqc']['library']
    wrapper:
        "file:///home/scruz/software/snakemake-wrappers/bio/fastqc"


rule samtools_stats:
    input:
        bam=f"{basedir}/results/{{plate}}/mapping/bwa/{{ref}}/{{sample}}.sorted.bam"
    output:
        f"{basedir}/results/{{plate}}/mapping/bwa/{{ref}}/samtools-stats/{{sample}}.txt"
    wrapper:
        "file:///home/scruz/software/snakemake-wrappers/bio/samtools/stats"
        
rule qual_stats_read_pos:
    input:
        bam=f"{basedir}/results/{{plate}}/mapping/bwa/{{ref}}/{{sample}}.sorted.bam",
        ref=rules.copy_reference.output
    output:
        f"{basedir}/results/{{plate}}/mapping/bwa/{{ref}}/readpos_stats/{{sample}}_readpos.stats"
    params:
        mem = "-Xmx3g"
    conda:
        "../envs/NGSEP.yaml"  
    shell:
        """
        java {params.mem} -jar {config[NGSEP][path]} BasePairQualStats \
            -o {output} -r {input.ref} {input.bam}
        """
    
rule plate_readpos_qc:
    input:
        readpos_stats = get_readpos_files
    output:
        merged_file = f"{basedir}/results/{{plate}}/mapping/bwa/{{ref}}/stats/readpos/plateQC_stats.csv",
        plot_file = f"{basedir}/results/{{plate}}/mapping/bwa/{{ref}}/stats/readpos/plateQC_readpos.pdf",
    run:
        import pandas as pd
        import matplotlib.pyplot as plt
        import numpy as np
        
        formated_outs = list()
        for readpos_out_path in input.readpos_stats:
            # magic -14 refers to file extension "_readpos.stats"
            sample = readpos_out_path.split('/')[-1][:-14]
            # ignore last 3 lines
            readpos_out = pd.read_csv(readpos_out_path, skipfooter=3, sep='\t', header=None, engine='python')
            if readpos_out.shape[0] > 1:
                readpos_out['sample'] = sample
                readpos_out['mismatch_mm'] = (readpos_out[1]/readpos_out[3])*100
                readpos_out['mismatch_um'] = (readpos_out[2]/readpos_out[4])*100
                formated_outs.append(readpos_out)
            else:
                print("check sample {sample} seems to be empty the readpos stats file".format(sample = sample))
        merged_readpos = pd.concat(formated_outs, ignore_index = True)
        means = merged_readpos.groupby(0, as_index = False).mean(numeric_only=True)
        fig = plt.figure(figsize=(15,8))
        ax = fig.add_subplot(111)
        ax.plot(means[0], means['mismatch_mm'], label='Multialignments')
        ax.plot(means[0], means['mismatch_um'], label='Unique alignments')
        ax.set_xlabel('Read Position')
        ax.set_ylabel('Missmatch Percentage (/%)')
        xticks = list(np.arange(0, 151,2))
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks, rotation=45, size=8)
        plt.grid(axis='y', color='0.9')
        plt.grid(axis='x', color='0.8')
        plt.legend()
        plt.savefig(output.plot_file)     
        merged_readpos.to_csv(output.merged_file, header = None, index = False)
        
rule plate_efficiencies_stats:
    input:
        unpack(get_plate_efficiencies)
    output:
        plate_qc = f'{basedir}/results/{{plate}}/mapping/bwa/{{ref}}/stats/plateQC.json'
    run:
        import os
        import json
        
        plate_qc_data = dict()
        # get total counts of stacks output
        data=dict()
        flag=False
        with open(input.stacks_log,'r') as f:
            for line in f:
                if line.startswith('BEGIN total_raw_read_counts'):
                    flag=True
                elif line.strip().endswith('END total_raw_read_counts'):
                    flag=False
                elif flag:
                    i_data = line[:-1].split('\t')
                    data[i_data[0]] = i_data[1:]
        
        plate_qc_data['total_reads'] = int(data['Total Sequences'][0])
        plate_qc_data['reads_not_found'] = plate_qc_data['total_reads'] - int(data['Retained Reads'][0])
        # aggregate mapping stats of all samples
        
        plate_qc_data['Gbp_non_adapters'] = 0
        plate_qc_data['Gbp_mapped'] = 0
        
        flag=False
        for i_file in input.samtools_stats:
            samtools_data = dict()
            with open(i_file,'r') as f:
                    for line in f:
                        if line.startswith('SN'):
                            flag=True
                            i_data = line[3:].split(':')
                            samtools_data[i_data[0]] = i_data[1]
                        elif not line.startswith('SN'):
                            flag=False
                            
            plate_qc_data['Gbp_non_adapters'] += int(samtools_data["total length"].split("\t")[1])
            plate_qc_data['Gbp_mapped'] += int(samtools_data["bases mapped (cigar)"].split("\t")[1])
        plate_qc_data['Gbp_non_adapters'] = plate_qc_data['Gbp_non_adapters']/1000000000
        plate_qc_data['Gbp_mapped'] = plate_qc_data['Gbp_mapped']/1000000000
        with open(output.plate_qc, "w") as write_file:
            json.dump(plate_qc_data, write_file, indent=4)
        
rule multiqc:
    input: 
        get_multiqc_files
    output:
        html=f"{basedir}/results/{{plate}}/mapping/bwa/{{ref}}/stats/multiqc.html",
        data_dir=directory(f"{basedir}/results/{{plate}}/mapping/bwa/{{ref}}/stats/multiqc_data")
    params: 
        use_input_files_only=True
    log:
        f"{basedir}/results/{{plate}}/logs/{{ref}}_multiqc.log",
    wrapper:
        "file:///home/scruz/software/snakemake-wrappers/bio/multiqc"

rule plate_mapping_stats:
    input:
        mapping_stats = f'{basedir}/results/{{plate}}/mapping/bwa/{{ref}}/stats/multiqc_data'
    output:
        plate_stats_plot = f'{basedir}/results/{{plate}}/mapping/bwa/{{ref}}/stats/plate_stats_plot.pdf',
    run:
        import matplotlib.pyplot as plt
        import pandas as pd
        import seaborn as sns
        import numpy as np
        barcodes = pd.read_csv(config['barcodes'])
        
        mapping_stats = pd.read_csv(input.mapping_stats + "/multiqc_general_stats.txt", sep='\t')
        mapping_stats.rename(columns = {'Sample': 'line_id'}, inplace = True)
        mapping_stats['log_MQ0'] = mapping_stats['Samtools: stats_mqc-generalstats-samtools_stats-reads_MQ0_percent']+1
        mapping_responses = ['log_MQ0',
                            'Samtools: stats_mqc-generalstats-samtools_stats-reads_mapped_percent',
                            'Samtools: stats_mqc-generalstats-samtools_stats-error_rate']
                            
        
        plate_data = sample_units.merge(barcodes, on = 'barcode')
        plate_data = plate_data.merge(mapping_stats, on = 'line_id' )
        
        
        fig, axs = plt.subplots(3, 1, figsize=(8, 10))
        for response, ax in zip(mapping_responses, axs):
            average = pd.pivot_table(data= plate_data, columns='well_r', index = 'well_c', values = response)
            sns.heatmap(average, cmap='viridis', ax = ax)
            ax.set_title(response)
        plt.suptitle(wildcards.plate + " ref:" + wildcards.ref)
        plt.tight_layout()
        plt.savefig(output.plate_stats_plot)

rule vcf_summary:
    input:
        snps_vcf = f"{basedir}/results/{{plate}}/variant_calling/NGSEP/{{ref}}/{{plate}}_merged_Q{{qual}}_Dp{{dp}}_MAF{{maf}}_sorted.vcf.gz"
    output:
        summary = f"{basedir}/results/{{plate}}/mapping/bwa/{{ref}}/{{plate}}_merged_Q{{qual}}_Dp{{dp}}_MAF{{maf}}_sorted_summary.txt",
    conda:
        'ngs'
    shell:
        """
        bcftools stats -s - -d 0,100,2 {input.snps_vcf} > {output.summary}
        """

rule plate_snp_stats:
    input:
        snp_stats = f"{basedir}/results/{{plate}}/mapping/bwa/{{ref}}/{{plate}}_merged_Q{{qual}}_Dp{{dp}}_MAF{{maf}}_sorted_summary.txt",
    output:
        plate_stats_plot = f"{basedir}/results/{{plate}}/mapping/bwa/{{ref}}/stats/{{plate}}_merged_Q{{qual}}_Dp{{dp}}_MAF{{maf}}_stats_plot.pdf"
        
    run:
        import matplotlib.pyplot as plt
        import pandas as pd
        import seaborn as sns
        import numpy as np

        barcodes = pd.read_csv(config['barcodes'])
        
        f = open(input.snp_stats, "r")
        lines = f.readlines()

        general_lines = [line.replace('\n', '').split('\t') for line in lines if re.match('^SN', line) != None]
        n_snps = int(general_lines[3][-1])

        taxa_lines = [line.replace('\n', '').split('\t') for line in lines if re.match('^PSC', line) != None]
        taxa_df = pd.DataFrame(taxa_lines)
        taxa_df.columns = ['PSC','id','line_id','nRefHom','nNonRefHom','nHets','nTransitions','nTransversions','nIndels','average_depth','nSingletons','nHapRef','nHapAlt','nMissing']
        taxa_df['total_variants'] = taxa_df['nRefHom'].astype(int) + taxa_df['nNonRefHom'].astype(int) + taxa_df['nHets'].astype(int)
        taxa_df['HO'] = taxa_df['nHets'].astype(int)/taxa_df['total_variants']
        taxa_df['TvTs'] =taxa_df['nTransitions'].astype(int) / taxa_df['nTransversions'].astype(int) 
        taxa_df['average_depth'] = taxa_df['average_depth'].astype(float)
        taxa_df['log_total_variants'] = np.log(taxa_df['total_variants']+1)

        plate_data = sample_units.merge(barcodes, on = 'barcode')
        plate_data = plate_data.merge(taxa_df, on = 'line_id' )
        
        
        snp_responses = ['log_total_variants', 'average_depth', 'HO', 'TvTs']
        
        fig, axs = plt.subplots(4, 1, figsize=(8, 10))
        for response, ax in zip(snp_responses, axs):
            print(response)
            average = pd.pivot_table(data= plate_data, columns='well_r', index = 'well_c', values = response)
            sns.heatmap(average, cmap='viridis', ax = ax)
            ax.set_title(response)
        plt.suptitle(wildcards.plate + " ref:" + wildcards.ref)
        plt.tight_layout()
        plt.savefig(output.plate_stats_plot)


rule bcf_stats:
    input:
        f"{basedir}/results/{{plate}}/variant_calling/NGSEP/{{ref}}/second_variant_call_plate/{{sample}}_bwa_NGSEP.vcf.gz"
    output:
        f'{basedir}/results/{{plate}}/mapping/bwa/{{ref}}/stats/by_sample/bcftools/{{sample}}.vcf.stats'
    log:
        f'{basedir}/results/{{plate}}/mapping/bwa/{{ref}}/stats/by_sample/bcftools/log/{{sample}}.vcf.stats.log'
    params:
        "",
    wrapper:
        "file:///home/scruz/software/snakemake-wrappers/bio/bcftools/stats"


rule multiqc_by_sample:
    input:
        samtools=f"{basedir}/results/{{plate}}/mapping/bwa/{{ref}}/samtools-stats/{{sample}}.txt",
        bcftools=f'{basedir}/results/{{plate}}/mapping/bwa/{{ref}}/stats/by_sample/bcftools/{{sample}}.vcf.stats'
    output:
        html=f"{basedir}/results/{{plate}}/mapping/bwa/{{ref}}/stats/by_sample/multiqc/{{sample}}/{{sample}}.multiqc_report.html",
        data_dir=directory(f"{basedir}/results/{{plate}}/mapping/bwa/{{ref}}/stats/by_sample/multiqc/{{sample}}/multiqc_data")
    params:
        extra=lambda wildcards: "-i 'Quality report for sample {sample}'".format(sample=wildcards.sample),
        use_input_files_only=True,
    log:
        f"{basedir}/results/{{plate}}/mapping/bwa/{{ref}}/stats/by_sample/multiqc/{{sample}}/{{sample}}.multiqc.log",
    threads: 4
    resources:
        mem_mb=1024
    wrapper:
        "file:///home/scruz/software/snakemake-wrappers/bio/multiqc"

