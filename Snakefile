

import pandas as pd
import os
import glob

configfile: 'config.yaml'

# take information from config
PATH_FASTQ = config["path"]["fastq"]
PATH_FUSION = config["path"]["fusion"]
PATH_QC = config["path"]["qc"]
PATH_LOG = config["path"]["log"]
PATH_BAM = config['path']['bam']

REFLIB = config['star_fusion']['reflib']
INDEX = config['star_fusion']['index']

# some options are taken from environmental variables
if config['Use_global'] is True:
    for env in config['environmental']:
        if config['environmental'][env] in os.environ:
            globals()[env] = os.environ[config['environmental'][env]]
            print(env + ' is overrived by enviromental variable ' + config['environmental'][env] + ': ' + globals()[env])


# if samples table exists, take fastq listed in the table,
# otherwise, it will just take every fastqs in PATH_FASTQ
if os.path.isfile('samples.tsv'):
    samples = pd.read_csv('samples.tsv', sep='\t')
    IDs = samples.IDs
    print("obtaining samples from config.yaml")
else:
    print("obtaining samples from the path : " + PATH_FASTQ)
    Temp = glob.glob(PATH_FASTQ+"/*.fastq.gz")
    IDs = [os.path.splitext(os.path.splitext(os.path.basename(tmp))[0])[0] for tmp in Temp]

print('found ' + str(len(IDs)) + ' fastq files to be processed')

rule all:
    input:
        expand(PATH_FUSION+"/{sample}.fusion_candidate", sample=IDs)

include: "rules/star_fusion.snakefile"
