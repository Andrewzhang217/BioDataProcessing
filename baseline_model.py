from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import argparse


def correct_error():
    return None


def parse_paf(paf_file):
    f = open(paf_file, 'r')
    lines = f.readlines()
    f.close()
    align_map = defaultdict(None)
    for line in lines:
        copy = line.strip()
        (query_seq_name, query_seq_len, query_start, query_end, strand, target_seq_name,
         target_seq_len, target_start, target_end, num_residue_matchm, align_block_len, map_quality) = copy.split(' ')
        align_map[query_seq_name] = (query_seq_len, query_start, query_end, strand, target_seq_name,
                                     target_seq_len, target_start, target_end, num_residue_matchm, align_block_len,
                                     map_quality)
        print(query_seq_name)

    return align_map

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Baseline model for error correction')
    parser.add_argument('-i', '--input', type=str, help='path of the input PAF file.')
    parser.add_argument('-o', '--output', type=str, help='path of error corrected reads file')
    args = parser.parse_args()

    align_map = parse_paf(args.input)
    print('Finished')
