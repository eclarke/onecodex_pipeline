from pprint import pprint
import glob
import onecodex_api as ocx

configfile: "snakemake_config.yaml"

config_default = {
    'data_fp': 'data_files/paired/',
    'output_fp': 'data_files/one_codex/',
    'api_key': os.environ.get("ONE_CODEX_API_KEY")
}

config_default['samples'] = glob.glob(config_default['data_fp'] + "*.fastq")

update_config(config_default, config)
config = config_default

all_targets = [
    config['output_fp'] + os.path.basename(s).strip(".fastq") + ".ocx.tsv"
    for s in config['samples']
]

rule all:
    input:
        all_targets

rule get_taxa:
    input:
        config['data_fp'] + "{sample}.fastq"
    output:
        config['output_fp'] + "{sample}.ocx.tsv"
    run:
        ocx.get_taxa_in_sample(input.sample, config['output_fp'], config['api_key'])

