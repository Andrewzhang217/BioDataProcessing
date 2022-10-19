from ast import arg
import re
from time import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import edlib
from pathlib import Path
from collections import defaultdict, Counter
import time
import argparse
from concurrent.futures import ProcessPoolExecutor
import re
import sys

def get_reads(input_file, output_file1, output_file2):
    file_type = Path(input_file).suffix
    file_type = file_type[1:]
    if file_type in ['fq']:
        file_type = 'fastq'

    seq_lst1 = []
    seq_lst2 = []
    for record in SeqIO.parse(input_file, file_type):
        # print(record.id)
        # break
        if '_' in record.id :
            underscore_idx = record.id.index('_')
            hap_id = int(record.id[underscore_idx  + 1: ])
        
        else :
            desc = record.description
            desc = desc.split(" ")
            hap_id = int(desc[1][-1])        

        if hap_id == 1:
            seq_lst1.append(record) 
        else:
            seq_lst2.append(record)    

    SeqIO.write(seq_lst1, output_file1, 'fasta')
    SeqIO.write(seq_lst2, output_file2, 'fasta')
    return None


def parse_args():
    parser = argparse.ArgumentParser(description='Baseline model for error correction')
    parser.add_argument('-i', '--input1', type=str, help='path of the mixed FASTA/Q file.')
    parser.add_argument('-o1', '--output1', type=str, help='path1 of the separated FASTA/Q file.')
    parser.add_argument('-o2', '--output2', type=str, help='path2 of the separated FASTA/Q file.')
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    args = parse_args()
    get_reads(args.input1, args.output1, args.output2)
