from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import edlib
from pathlib import Path
from collections import defaultdict
import argparse


def correct_error(cigar_list, out_path):

    for query_read in cigar_list:
        print(query_read) 
        break
    #     for cigar in cigar_list[query_read]:
    #         print(cigar)
    # for record in SeqIO.parse(in_path, file_type):
    #     r = rle(str(record.seq))
    #     sr = SeqRecord(Seq(r))
    #     sr.id = record.id
    #     sr.name = record.name
    #     sr.description = record.description
    #     rles.append(sr)
    for cigar in cigar_list['12_96972_aligned_0_F_20_10782_26']:
        print(cigar)
    #SeqIO.write(ec_reads, out_path, 'fasta')
    return None


def generate_cigar(align_map, reads):
    cigar_list = defaultdict(list)
    cnt = 0
    for query_seq in align_map: 
        # print(cnt)
        if cnt == 100: break
        cnt += 1
        for overlap in align_map[query_seq]:
            query_start = int (overlap[1])
            query_end = int (overlap[2])
            query_overlap = reads[query_seq][query_start:query_end + 1]
            target_name = overlap[4]
            target_start = int (overlap[6])
            target_end = int (overlap[7])
            target_overlap = reads[target_name][target_start:target_end + 1]
            # print(len(target_overlap))
            # print(len(query_overlap))
            cigar = edlib.align(target_overlap, query_overlap, task = 'path')
            # cigar_list[query_seq].append(cigar)

    print(len(cigar_list))
    return cigar_list


def get_reads(input_file):
    file_type = Path(input_file).suffix
    file_type = file_type[1:]
    if file_type in ['fq']:
        file_type = 'fastq'

    records = list(SeqIO.parse(input_file, file_type))
    for i in range(50):
        print(len(records[i].seq))
        cigar = edlib.align(records[i].seq, records[i + 1].seq, task = "path")
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
        (query_seq_name, query_seq_len, query_start, query_end, strand, overlap_name,
         overlap_len, target_start, target_end, num_residue_match, align_block_len, map_quality) = copy.split('\t')
        dst[query_seq_name].append((query_seq_len, query_start, query_end, strand, overlap_name,
                                overlap_len, target_start, target_end, num_residue_match, align_block_len,
                                map_quality))
    return dst


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Baseline model for error correction')
    parser.add_argument('-i', '--input', type=str, help='path of the input FASTA/Q file.')
    parser.add_argument('-p', '--paf', type=str, help='path of the input PAF file.')
    parser.add_argument('-o', '--output', type=str, help='path of error corrected reads file')
    args = parser.parse_args()

    reads = get_reads(args.input)
    #align_map = parse_paf(args.paf)
    print("before")
    # cigars = generate_cigar(align_map, reads)
    # correct_error(cigars, args.output)
    print('Finished')
