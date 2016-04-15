import os
import glob

from ruffus import *
import ruffus.cmdline as cmdline

DATA_FP="/media/8TB_PLAYGROUND/home/ecl/ext/100_SCID/103_Virome/data_files/"

fwd_files = glob.glob(DATA_FP + "*R1.fastq.gz")
rev_files = glob.glob(DATA_FP + "*R2.fastq.gz")

starting_files = list(zip(fwd_files, rev_files))

parser = cmdline.get_argparse(description="Pairs reads using PEAR")
options = parser.parse_args()
    
@transform(starting_files, suffix("R1.fastq.gz"), "assembled.fastq")
def pair_reads(input_files, output_file):
    os.system("pear -j 4 -f {} -r {} -o {}".format(input_files[0], input_files[1], output_file))

cmdline.run(options)
    
