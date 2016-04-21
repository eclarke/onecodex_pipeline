from pprint import pprint
import csv
import glob

from snakemake.utils import R

import onecodex_api as ocx

configfile: "snakemake_config.yaml"

# Identify all samples based on data_fp config option
config['samples'] = [
    os.path.basename(s)
    for s in glob.glob(config['data_fp'] + "*.fastq")
]

# Dynamically generate targets based on samples
all_targets = [
    config['output_fp'] + os.path.basename(s) + ".taxa.tsv"
    for s in config['samples']
]


rule all:
    input:
        all_samples = config['output_fp'] + "all_samples.tsv",
        all_taxa = config['output_fp'] + "all_taxa.tsv",
        sample_info = config['output_fp'] + "sample_summary.tsv"


rule get_taxa:
    input:
        config['data_fp'] + "{sample}",
    output:
        config['output_fp'] + "{sample}.taxa.tsv"
    run:
        ocx.get_taxa_in_sample(os.path.basename(input[0]), config['output_fp'], config['api_key'])

        
rule get_sample_info:
    input:
        expand(config['data_fp'] + "{sample}", sample=config['samples'])
    output:
        config['output_fp'] + "sample_summary.tsv"
    run:
        with open(output[0], 'w') as out:
            writer = csv.DictWriter(
                out,
                ['sample_filename','reference_id','n_reads','p_mapped'],
                extrasaction='ignore',
                delimiter='\t'
            )
            writer.writeheader()
            for s in input:
                s = os.path.basename(s)
                writer.writerows(
                    ocx.get_ocx_analysis_for_sample(s, config['api_key']))

                
rule aggregate_taxa:
    input:
        all_targets
    output:
        all_samples = config['output_fp'] + "all_samples.tsv",
        all_taxa = config['output_fp'] + "all_taxa.tsv"
    script:
        "aggregate_samples.R"


