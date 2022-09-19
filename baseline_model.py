from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import edlib
from pathlib import Path
from collections import defaultdict
import argparse


def correct_error(cigar_map, out_path):

    for query_read in cigar_map:
        for cigar in query_read:
            print(cigar)
    # for record in SeqIO.parse(in_path, file_type):
    #     r = rle(str(record.seq))
    #     sr = SeqRecord(Seq(r))
    #     sr.id = record.id
    #     sr.name = record.name
    #     sr.description = record.description
    #     rles.append(sr)
    #
    # SeqIO.write(ec_reads, out_path, 'fasta')
    return None


def generate_cigar(align_map, reads):
    cigar_map = defaultdict(lambda: defaultdict(None))
    for query_seq in align_map:
        start = query_seq[0][1]
        end = query_seq[0][2]
        exact_query_seq = reads[query_seq][start:end]
        for target_seq in query_seq:
            target_name = target_seq[4]
            target_start = target_seq[6]
            target_end = target_seq[7]
            exact_target_seq = reads[target_name][target_start:target_end]
            cigar = edlib.align(exact_query_seq, exact_target_seq)
            cigar_map[query_seq] += cigar

    print(len(cigar_map))
    return cigar_map


def get_reads(input_file):
    file_type = Path(input_file).suffix
    file_type = file_type[1:]
    if file_type in ['fq']:
        file_type = 'fastq'

    records = list(SeqIO.parse(input_file, file_type))
    records_map = defaultdict(None)
    for record in records:
        records_map[record.name] = record.seq
    return records_map


def parse_paf(paf_file):
    f = open(paf_file, 'r')
    lines = f.readlines()
    f.close()
    dst = defaultdict(list)
    for line in lines:
        copy = line.strip()
        (query_seq_name, query_seq_len, query_start, query_end, strand, target_seq_name,
         target_seq_len, target_start, target_end, num_residue_match, align_block_len, map_quality) = copy.split('\t')
        dst[query_seq_name].append((query_seq_len, query_start, query_end, strand, target_seq_name,
                                target_seq_len, target_start, target_end, num_residue_match, align_block_len,
                                map_quality))
    print(len(dst))
    return dst


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Baseline model for error correction')
    parser.add_argument('-i', '--input', type=str, help='path of the input FASTA/Q file.')
    parser.add_argument('-p', '--paf', type=str, help='path of the input PAF file.')
    parser.add_argument('-o', '--output', type=str, help='path of error corrected reads file')
    args = parser.parse_args()

    reads = get_reads(args.input)
    align_map = parse_paf(args.paf)
    cigars = generate_cigar(align_map, reads)
    correct_error(cigars, args.output)
    print('Finished')
