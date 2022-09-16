from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import edlib
from pathlib import Path
from collections import defaultdict
import argparse


def correct_error():
    return None

def generate_cigar(align_map):

def get_reads(input_file):
    file_type = Path(input_file).suffix
    file_type = file_type[1:]
    if file_type in ['fq']:
        file_type = 'fastq'

    records = list(SeqIO.parse(input_file, file_type))

    return records

def parse_paf(paf_file):
    f = open(paf_file, 'r')
    lines = f.readlines()
    f.close()
    align_map = defaultdict(None)
    for line in lines:
        copy = line.strip()
        (query_seq_name, query_seq_len, query_start, query_end, strand, target_seq_name,
         target_seq_len, target_start, target_end, num_residue_matchm, align_block_len, map_quality) = copy.split('\t')
        align_map[query_seq_name] = (query_seq_len, query_start, query_end, strand, target_seq_name,
                                     target_seq_len, target_start, target_end, num_residue_matchm, align_block_len,
                                     map_quality)
        print(query_seq_name)

    return align_map

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Baseline model for error correction')
    parser.add_argument('-i', '--input', type=str, help='path of the input FASTA/Q file.')
    parser.add_argument('-p', '--paf', type = str, help = 'path of the input PAF file.')
    parser.add_argument('-o', '--output', type=str, help='path of error corrected reads file')
    args = parser.parse_args()

    reads = get_reads(input)
    align_map = parse_paf(args.paf)
    print('Finished')
