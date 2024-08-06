# ---------------------------------------------------------------------------------------------
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version
from yaml import safe_load
import glob as glob
import random
import string

# ---------------------------------------------------------------------------------------------
min_version("5.18.0")

# ---------------------------------------------------------------------------------------------
# Config params
configfile: "config/config.yaml"

# ---------------------------------------------------------------------------------------------
# Resources params
with open(config["resources_config"], "r") as f:
    resources = safe_load(f)

# ---------------------------------------------------------------------------------------------
# Load and validate sequencing units
sequencing_units = pd.read_table(config["sequencing_units"], sep="\t")
validate(sequencing_units, schema="../schemas/units.schema.yaml")
sequencing_units['plate'] = sequencing_units.plate.astype(str)

# Load and validate sample units
sample_units = pd.read_table(config["sample_units"], sep="\t")
validate(sample_units, schema="../schemas/samples.schema.yaml")
sample_units['plate'] = sample_units.plate.astype(str)
sample_units = sample_units.merge(sequencing_units[['sequencing_unit_id', 'plate', 'fq1', 'fq2']], on="plate")

# Load and validate references
references = pd.read_table(config['references'], sep = '\t')
validate(references, schema="../schemas/references.schema.yaml")
references.set_index("ref_name", inplace=True)

# Set sequencing_units index
sequencing_units.set_index("plate", inplace=True)

# Print the config for debugging
print(config)



# ---------------------------------------------------------------------------------------------
# UTIL functions

def get_demultiplex_fastqs(wildcards):
    
    # Summary: Retrieve the paths to fastq files for a given sequencing unit.
    # Args: 
    #       wildcards (dict): Snakemake wildcards, should contain 'plate'.
    # Returns:
    #       dict: A dictionary containing paths to fastq files, with keys 'r1' and 'r2' if paired-end, otherwise only 'r1'.

    fastqs = sequencing_units.loc[wildcards.plate, ['fq1','fq2']]
    if not pd.isna(fastqs.fq2):
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    else:
        return {"r1": fastqs.fq1}

def get_rawread_folder_name(wildcards):
    
    # Summary: Get the folder name of raw reads based on the fastq file path.
    # Args: 
    #       wildcards (dict): Snakemake wildcards, should contain 'plate'.
    # Returns:
    #       str: The name of the folder containing raw reads.

    fastqs = sequencing_units.loc[wildcards.plate, ['fq1','fq2']]
    raw_read_folder = fastqs.fq1.split('/')[-2]
    return raw_read_folder

def get_demultiplex_params(wildcards):
    
    # Summary: Construct the command-line parameters for demultiplexing.
    # Args: 
    #       wildcards (dict): Snakemake wildcards, should contain 'plate'.
    # Returns:
    #       str: The command-line parameters for demultiplexing.

    fastqs = sequencing_units.loc[wildcards.plate, ['fq1','fq2']]
    if not pd.isna(fastqs.fq2):
        return "-1 {R1} -2 {R2}".format(R1=fastqs.fq1, R2=fastqs.fq2)
    else:
        return "-f {R1}".format(R1=fastqs.fq1)

def get_sample_fastq(wildcards):
    
    # Summary: Get the paths to sample fastq files after demultiplexing.
    # Args: 
    #       wildcards (dict): Snakemake wildcards, should contain 'plate' and 'sample'.
    # Returns:
    #       dict or str: A dictionary with paths to paired fastq files if paired-end, otherwise a string path to the single fastq file.
    
    checkpoint_output = checkpoints.demultiplex.get(**wildcards).output.outdir
    fastqs = sequencing_units.loc[wildcards.plate, ['fq1','fq2']]
    
    if not pd.isna(fastqs.fq2):
        sample_fastq = glob.glob(checkpoint_output + "/{sample}.[1,2].fq.gz".format(sample=wildcards.sample))
        return {"r1": sample_fastq[0], "r2": sample_fastq[1]}
    else:
        sample_fastq = glob.glob(checkpoint_output + "/{sample}.fq.gz".format(sample=wildcards.sample))
        return sample_fastq

def get_trimmed_reads(wildcards):
    
    # Summary: Get the paths to trimmed reads.
    # Args: 
    #       wildcards (dict): Snakemake wildcards, should contain 'plate' and 'sample'.
    # Returns:
    #       list: A list of paths to trimmed reads.
    
    fastqs = sequencing_units.loc[wildcards.plate, ['fq1','fq2']]
    if not pd.isna(fastqs.fq2):
        return expand(
            "results/{plate}/trimming/trimmomatic/paired/{sample}.{group}.fastq.gz",
            group=[1, 2],
            **wildcards
        )
    else:
        return expand(
            "results/{plate}/trimming/trimmomatic/single/{sample}.fastq.gz",
            **wildcards
        )

def get_read_group(wildcards):
    
    # Summary: Construct the read group string for alignment.
    # Args: 
    #       wildcards (dict): Snakemake wildcards, should contain 'sample'.
    # Returns:
    #       str: The read group string.
    
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample,
        platform="ILLUMINA",
    )

def get_big_temp(wildcards):
    
    # Summary: Get a large temporary directory path for jobs requiring more space.
    # Args: 
    #       wildcards (dict): Snakemake wildcards.
    # Returns:
    #       str: The path to a large temporary directory.
    
    if config['bigtmp']:
        if config['bigtmp'].endswith("/"):
            return config['bigtmp'] + "".join(random.choices(string.ascii_uppercase, k=12)) + "/"
        else:
            return config['bigtmp'] + "/" + "".join(random.choices(string.ascii_uppercase, k=12)) + "/"
    else:
        return tempfile.gettempdir()

def get_reference_fasta(wildcards):
    
    # Summary: Get the path to the reference fasta file.
    # Args: 
    #       wildcards (dict): Snakemake wildcards, should contain 'ref'.
    # Returns:
    #       str: The path to the reference fasta file.

    return references.loc[wildcards.ref, "ref_path"]

def get_sample_vcfs_by_plate_merge_variants(wildcards):
    
    # Summary: Get the paths to sample VCFs for the first variant calling.
    # Args: 
    #       wildcards (dict): Snakemake wildcards, should contain 'plate' and 'ref'.
    # Returns:
    #       list: A list of paths to sample VCFs.

    checkpoint_output = checkpoints.demultiplex.get(**wildcards).output.outdir
    sample_list = glob.glob(checkpoint_output + "/*[!rem]*.fq.gz")
    sample_names = list(set([s.split('/')[-1].split('.')[0] for s in sample_list]))

    vcfs = expand("results/{plate}/variant_calling/NGSEP/{ref}/first_variant_calling/{sample}_bwa_NGSEP.vcf.gz",
        plate=wildcards.plate,
        ref=wildcards.ref,
        sample=sample_names)

    return vcfs

def get_sample_vcfs_by_plate_merge_vcfs(wildcards):
    
    # Summary: Get the paths to sample VCFs for the second variant calling.
    # Args: 
    #       wildcards (dict): Snakemake wildcards, should contain 'plate' and 'ref'.
    # Returns:
    #       list: A list of paths to sample VCFs.

    checkpoint_output = checkpoints.demultiplex.get(**wildcards).output.outdir
    sample_list = glob.glob(checkpoint_output + "/*[!rem]*.fq.gz")
    sample_names = list(set([s.split('/')[-1].split('.')[0] for s in sample_list]))

    vcfs = expand("results/{plate}/variant_calling/NGSEP/{ref}/second_variant_call_plate/{sample}_bwa_NGSEP.vcf.gz",
        plate=wildcards.plate,
        ref=wildcards.ref,
        sample=sample_names)

    return vcfs

def get_GATK_HaplotypeCaller_params():
    
    # Summary: Construct the parameters for GATK HaplotypeCaller.
    # Returns:
    #       str: The command-line parameters for GATK HaplotypeCaller.
    
    annot = ' '.join(["-G {p}".format(p=p) for p in config["GATK"]['HaplotypeCaller']['G']])
    kmers = ' '.join(["-kmer-size {p}".format(p=p) for p in config["GATK"]['HaplotypeCaller']['kmer-size']])
    extra = ' '.join(["{p}".format(p=p) for p in config["GATK"]['HaplotypeCaller']['extra']])
    min_base_qual = "--min-base-quality-score {p}".format(p=config["GATK"]['HaplotypeCaller']['min_base-qual'])
    params = ' '.join([min_base_qual, annot, kmers, extra])
    return params

def get_GATK_CombineGVCFs_params():
    
    # Summary: Construct the parameters for GATK CombineGVCFs.
    # Returns:
    #       str: The command-line parameters for GATK CombineGVCFs.

    annot = ' '.join(["-G {p}".format(p=p) for p in config["GATK"]['CombineGVCFs']['G']])
    return annot

def get_gvcfs_DB(wildcards):
    
    # Summary: Get the paths to GVCFs for the GenomicsDBImport step.
    # Args: 
    #       wildcards (dict): Snakemake wildcards, should contain 'plate' and 'ref'.
    # Returns:
    #       list: A list of paths to GVCFs.
    
    checkpoint_output = checkpoints.demultiplex.get(**wildcards).output.outdir
    sample_list = glob.glob(checkpoint_output + "/*[!rem]*.fq.gz")
    sample_names = list(set([s.split('/')[-1].split('.')[0] for s in sample_list]))
    gvcfs_list = list()

    for isample in sample_names:
        gvcf = 'results/{plate}/variant_calling/GATK/{ref}/CombineGVCFs/{sample}.g.vcf.gz'.format(
            sample=isample, **wildcards)
        gvcfs_list.append(gvcf)
    return gvcfs_list

def get_gvcfs_by_sample(wildcards):
    
    # Summary: Get the paths to GVCFs for a sample, split by intervals.
    # Args: 
    #       wildcards (dict): Snakemake wildcards, should contain 'plate', 'ref', and 'sample'.
    # Returns:
    #       list: A list of paths to GVCFs.
    
    intervals = pd.read_csv("resources/{ref}/{ref}_intervals.txt".format(**wildcards), header=None)
    gvcfs = ['results/{plate}/variant_calling/GATK/{ref}/HaplotyeCaller/intervals/{interval}/{sample}_{interval}.g.vcf.gz'.format(
        interval=i[0], **wildcards) for n, i in intervals.iterrows()]
    return gvcfs

def get_GenomicsDBImport_params():

    # Summary: Construct the parameters for GATK GenomicsDBImport.
    # Returns:
    #       str: The command-line parameters for GATK GenomicsDBImport.

    params = ""
    for param in config['GATK']['GenomicsDBImport'].keys():
        if not isinstance(config['GATK']['GenomicsDBImport'][param], list):
            params += "--{param} {value} ".format(param=param, value=config['GATK']['GenomicsDBImport'][param])
        else:
            for option in config['GATK']['GenomicsDBImport'][param]:
                params += "{option} ".format(option=option)
    return params

def get_GenotypeGVCFs_params():
    
    # Summary: Construct the parameters for GATK GenotypeGVCFs.
    # Returns:
    #       str: The command-line parameters for GATK GenotypeGVCFs.
    
    params = ""
    for param in config['GATK']['GenotypeGVCFs'].keys():
        if param != "interval_length":
            if not isinstance(config['GATK']['GenotypeGVCFs'][param], list):
                params += "--{param} {value} ".format(param=param, value=config['GATK']['GenotypeGVCFs'][param])
            else:
                for option in config['GATK']['GenotypeGVCFs'][param]:
                    params += "{option} ".format(option=option)
    return params

def create_intervals(start, end, interval_size):
    
    # Summary: Create intervals for genomic regions.
    # Args: 
    #       start (int): Start position of the genomic region.
    #       end (int): End position of the genomic region.
    #       interval_size (int): Size of each interval.
    # Returns:
    #       list: A list of tuples, each representing an interval (start, end).
    
    intervals = []
    current = start

    while current < end:
        next_value = min(current + interval_size, end)
        intervals.append((current, next_value - 1))
        current = next_value
    
    if current == end:
        pass
    else:
        intervals.append((current, end))
    return intervals

def get_interval_raw_vcfs(wildcards):
    
    # Summary: Get the paths to raw VCFs for each interval based on the reference genome.
    # Args: 
    #       wildcards (dict): Snakemake wildcards, should contain 'plate', 'ref'.
    # Returns:
    #       list: A list of paths to raw VCFs for each interval.

    fai = "resources/{ref}/{ref}.fasta.fai".format(**wildcards)

    with open(fai, 'r') as file:
        lines = file.readlines()
    file.close()
    
    intervals_list = list()
    for line in lines:
        data = line.strip().split('\t')
        chrom = data[0]
        length = int(data[1])
        intervals = create_intervals(1, length, int(config['GATK']['GenotypeGVCFs']['interval_length']))

        vcfs = ["results/{plate}/variant_calling/GATK/{ref}/GenotypeGVCFs/{chrom}/{interval_i}-{interval_e}.vcf.gz".format(
            chrom=chrom,
            interval_i=str(i[0]),
            interval_e=str(i[1]),
            **wildcards) for i in intervals]

        intervals_list.extend(vcfs)
    print(intervals_list)
    return intervals_list

def get_SelectVariants_params():
    
    # Summary: Construct the parameters for GATK SelectVariants.
    # Returns:
    #       str: The command-line parameters for GATK SelectVariants.

    params = ""
    for param in config['GATK']['SelectVariants'].keys():
        if param != "interval_length":
            if isinstance(config['GATK']['SelectVariants'][param], list):
                for option in config['GATK']['SelectVariants'][param]:
                    params += "{option} ".format(option=option)
            else:
                params += "--{param} {value} ".format(param=param, value=config['GATK']['SelectVariants'][param])
    return params


# ---------------------------------------------------------------------------------------------
wildcard_constraints:
    sq_unit = "|".join(sequencing_units['sequencing_unit_id'].unique()),
    plate="|".join(sample_units['plate'].unique()),
    sample="|".join(sample_units['line_id'].unique()),
    ref="|".join(references.index),
    interval_i=r"\d+",
    interval_e=r"\d+"
