import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version
from yaml import safe_load
import glob as glob
import random
import string

min_version("5.18.0")

# Config params
configfile: "config/config.yaml"


# Resources params
with open(config["resources_config"], "r") as f:
    resources = safe_load(f)


sequencing_units = pd.read_table(config["sequencing_units"], sep="\t")
validate(sequencing_units, schema="../schemas/units.schema.yaml")

sequencing_units['plate'] = sequencing_units.plate.astype(str)

sample_units = pd.read_table(config["sample_units"], sep="\t")
validate(sample_units, schema="../schemas/samples.schema.yaml")
sample_units['plate'] = sample_units.plate.astype(str)

sample_units = sample_units.merge(sequencing_units[['sequencing_unit_id', 'plate', 'fq1', 'fq2']], on="plate")


references = pd.read_table(config['references'], sep = '\t')
validate(references, schema="../schemas/references.schema.yaml")
references.set_index("ref_name", inplace = True)


sequencing_units.set_index("plate", inplace = True)

print(config)

# UTIL functions
def get_demultiplex_fastqs(wildcards):
    fastqs = sequencing_units.loc[wildcards.plate, ['fq1','fq2']]
    if not pd.isna(fastqs.fq2):
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    else:
        return {"r1": fastqs.fq1}

def get_rawread_folder_name(wildcards):
    fastqs = sequencing_units.loc[wildcards.plate, ['fq1','fq2']]
    raw_read_folder = fastqs.fq1.split('/')[-2]
    return raw_read_folder

def get_demultiplex_params(wildcards):
    fastqs = sequencing_units.loc[wildcards.plate, ['fq1','fq2']]
    if not pd.isna(fastqs.fq2):
        return "-1 {R1} -2 {R2}".format(R1=fastqs.fq1, R2=fastqs.fq2)
    else:
        return "-f {R1}".format(R1=fastqs.fq1)
    
def get_sample_fastq(wildcards):
    checkpoint_output = checkpoints.demultiplex.get(**wildcards).output.outdir
    fastqs = sequencing_units.loc[wildcards.plate, ['fq1','fq2']]
    
    if not pd.isna(fastqs.fq2):
        sample_fastq = sorted(glob.glob(checkpoint_output + "/{sample}.[1,2].fq.gz".format(sample = wildcards.sample)))
        return {"r1": sample_fastq[0], "r2": sample_fastq[1]}
    else:
        sample_fastq = glob.glob(checkpoint_output + "/{sample}.fq.gz".format(sample = wildcards.sample))
        return sample_fastq

def get_trimmed_reads(wildcards):
    fastqs = sequencing_units.loc[wildcards.plate, ['fq1','fq2']]
    if not pd.isna(fastqs.fq2):
        return expand(
            "results/{plate}/trimming/trimmomatic/paired/{sample}_trimmed.{group}.fastq.gz",
            group = [1,2],
            **wildcards
        )
    else:
        return expand(
            "results/{plate}/trimming/trimmomatic/single/{sample}_trimmed.fastq.gz",
            **wildcards
        )

def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample,
        platform="ILLUMINA",
    )

def get_big_temp(wildcards):
    """Sets a temp dir for rules that need more temp space that is typical on some cluster environments. Defaults to system temp dir."""
    if config['bigtmp']:
        if config['bigtmp'].endswith("/"):
            return config['bigtmp'] + "".join(random.choices(string.ascii_uppercase, k=12)) + "/"
        else:
            return config['bigtmp'] + "/" + "".join(random.choices(string.ascii_uppercase, k=12)) + "/"
    else:
        return tempfile.gettempdir()

def get_reference_fasta(wildcards):
    return(references.loc[wildcards.ref, "ref_path"])


def get_sample_vcfs_by_plate_merge_variants(wildcards):
    checkpoint_output = checkpoints.demultiplex.get(**wildcards).output.outdir
    sample_list = glob_wildcards(os.path.join(checkpoint_output, "{sample}.1.fq.gz")).sample
    sample_list = [i for i in sample_list if '.rem' not in i]
    vcfs = expand("results/{plate}/variant_calling/NGSEP/{ref}/first_variant_calling/{sample}_bwa_NGSEP.vcf.gz",
        plate = wildcards.plate,
        ref = wildcards.ref,
        sample = sample_list )    
    
    return vcfs

def get_sample_vcfs_by_plate_merge_vcfs(wildcards):
    checkpoint_output = checkpoints.demultiplex.get(**wildcards).output.outdir
    sample_list = glob_wildcards(os.path.join(checkpoint_output, "{sample}.1.fq.gz")).sample
    sample_list = [i for i in sample_list if '.rem' not in i]
    vcfs = expand("results/{plate}/variant_calling/NGSEP/{ref}/second_variant_call_plate/{sample}_bwa_NGSEP.vcf.gz",
        plate = wildcards.plate,
        ref = wildcards.ref,
        sample = sample_list )
    return vcfs

def get_GATK_HaplotypeCaller_params():
    # Annotation params
    annot = ' '.join(["-G {p}".format(p=p) for p in config["GATK"]['HaplotypeCaller']['G']])
    kmers = ' '.join(["-kmer-size {p}".format(p=p) for p in config["GATK"]['HaplotypeCaller']['kmer-size']])
    extra = ' '.join(["{p}".format(p=p) for p in config["GATK"]['HaplotypeCaller']['extra']])
    min_base_qual = "--min-base-quality-score {p}".format(p=config["GATK"]['HaplotypeCaller']['min_base-qual'])
    params = ' '.join([min_base_qual, annot, kmers, extra])
    return params

def get_GATK_CombineGVCFs_params():
    # Annotation params
    annot = ' '.join(["-G {p}".format(p=p) for p in config["GATK"]['CombineGVCFs']['G']])
    return annot

def get_gvcfs_DB(wildcards):
    checkpoint_output = checkpoints.demultiplex.get(**wildcards).output.outdir
    sample_list = glob.glob(checkpoint_output + "/*[!rem]*.fq.gz")
    sample_names = list(set([s.split('/')[-1].split('.')[0] for s in sample_list]))

    gvcfs_list = list()
    for isample in sample_names:
        gvcf = 'results/{plate}/variant_calling/GATK/{ref}/CombineGVCFs/{sample}.g.vcf.gz'.format(
            sample = isample,**wildcards)
        gvcfs_list.append(gvcf)
    return gvcfs_list

def get_gvcfs_by_sample(wildcards):
    intervals = pd.read_csv("resources/{ref}/{ref}_intervals.txt".format(**wildcards), header = None)
    gvcfs = ['results/{plate}/variant_calling/GATK/{ref}/HaplotyeCaller/intervals/{interval}/{sample}_{interval}.g.vcf.gz'.format(
        interval = i[0],
        **wildcards) for n,i in intervals.iterrows()]
    return gvcfs

def get_GenomicsDBImport_params():
    params = ""
    for param in config['GATK']['GenomicsDBImport'].keys():
        if type(config['GATK']['GenomicsDBImport'][param]) != list:
            params += "--{param} {value} ".format(param = param,
                value = config['GATK']['GenomicsDBImport'][param]
            )
        else:
            for option in config['GATK']['GenomicsDBImport'][param]:
                params += "{option} ".format(option = option)
    return params

def get_GenotypeGVCFs_params():
    params = ""
    for param in config['GATK']['GenotypeGVCFs'].keys():
        if param != "interval_length":
            if type(config['GATK']['GenotypeGVCFs'][param]) != list:
                params += "--{param} {value} ".format(param = param,
                    value = config['GATK']['GenotypeGVCFs'][param])
            else:
                for option in config['GATK']['GenotypeGVCFs'][param]:
                    params += "{option} ".format(option = option)
    return params

def create_intervals(start, end, interval_size):
    intervals = []
    current = start

    while current < end:
        next_value = min(current + interval_size, end)
        intervals.append((current, next_value - 1))
        current = next_value

    # Handle the remaining value if it exists
    if current == end:
        pass
    else:
        intervals.append((current, end))

    return intervals


def get_interval_raw_vcfs(wildcards):
    
    # Read genome fai to infer the intervals
    #fai = "resources/{ref}/{ref}.fasta.fai".format(**wildcards)
    fai = "resources/{ref}/{ref}.fasta.fai".format(**wildcards)
    # Open the file in read mode
    with open(fai, 'r') as file:
        # Read all lines from the file
        lines = file.readlines()
    file.close()
    
    intervals_list = list()
    # Print each line
    for line in lines:
        data = line.strip().split('\t')
        chrom = data[0]
        length = int(data[1])
        
        intervals = create_intervals(1,length, int(config['GATK']['GenotypeGVCFs']['interval_length']))

        vcfs = ["results/{plate}/variant_calling/GATK/{ref}/GenotypeGVCFs/{chrom}/{interval_i}-{interval_e}.vcf.gz".format(
            chrom = chrom,
            interval_i = str(i[0]),
            interval_e = str(i[1]),
            **wildcards
        ) for i in intervals]

        intervals_list.extend(vcfs)
    print(intervals_list)
    return intervals_list


wildcard_constraints:
    sq_unit = "|".join(sequencing_units['sequencing_unit_id'].unique()),
    plate="|".join(sample_units['plate'].unique()),
    sample="|".join(sample_units['line_id'].unique()),
    ref = "|".join(references.index),
    interval_i = "\d+",
    interval_e = "\d+",


####################

def get_input_all():
    checkpoint_output = checkpoints.demultiplex.get(plate = "44", ref = "v21").output.outdir
    sample_list = glob.glob(checkpoint_output + "/*[!rem]*.fq.gz")
    sample_names = list(set([s.split('/')[-1].split('.')[0] for s in sample_list]))
    bams = expand("results/{plate}/mapping/bwa/{ref}/{sample}.sorted.bam",
        plate = "44",
        ref = "v21",
        sample = sample_names )
    return bams

def get_sample_multiqc():
    stacks_dir = "~/pipelines/modules/dna-variant-calling/results/44/demultiplexing/stacks"
    sample_list = [file for file in os.listdir(stacks_dir) if file.endswith(".fq.gz")]
    sample_names = [file.split('.')[0] for file in sample_list]
    return sample_names

def get_sample_qc_by_plate(wildcards):
    checkpoint_output = checkpoints.demultiplex.get(**wildcards).output.outdir
    sample_list = glob_wildcards(os.path.join(checkpoint_output, "{sample}.1.fq.gz")).sample
    sample_list = [i for i in sample_list if '.rem' not in i]
    qcs = expand("results/{plate}/quality_check/{ref}/multiqc_report_per_sample/{sample}/multiqc_report.html",
        plate = wildcards.plate,
        ref = wildcards.ref,
        sample = sample_list )
    return qcs


