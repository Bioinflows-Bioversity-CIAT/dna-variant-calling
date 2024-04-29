import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version
from yaml import safe_load
import glob as glob

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
        sample_fastq = glob.glob(checkpoint_output + "/{sample}.[1,2].fq.gz".format(sample = wildcards.sample))
        return {"r1": sample_fastq[0], "r2": sample_fastq[1]}
    else:
        sample_fastq = glob.glob(checkpoint_output + "/{sample}.fq.gz".format(sample = wildcards.sample))
        return sample_fastq

def get_trimmed_reads(wildcards):
    fastqs = sequencing_units.loc[wildcards.plate, ['fq1','fq2']]
    if not pd.isna(fastqs.fq2):
        return expand(
            "results/{plate}/reads/trimmed/paired/{sample}.{group}.fastq.gz",
            group = [1,2],
            **wildcards
        )
    else:
        return expand(
            "results/{plate}/reads/trimmed/single/{sample}.fastq.gz",
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

def get_sample_by_plate_test(ref, plate):
    plate_df = sample_units[sample_units['plate'] == plate]
    vcfs = list()
    for n, row in plate_df.iterrows():
        i_vcf = 'results/{plate}/mapping/{ref}/bwa/mapping/{sample}.sorted.bam'.format(plate=row.plate, 
                                                                                 sample=row.line_id,
                                                                                 ref=ref)
        vcfs.append(i_vcf)
    return vcfs

def get_reference_fasta(wildcards):
    return(references.loc[wildcards.ref, "ref_path"])