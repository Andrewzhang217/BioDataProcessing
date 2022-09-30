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


def correct_error(cigar_list, reads, out_path):
    seq_lst = []
    for target_read in cigar_list:
        # uncorrected
        ec_read = reads[target_read].seq
        # reads before error correction
        # seq_lst.append(reads[target_read])
        print(len(ec_read))

        # Add original base into frequency
        freq = []
        for i, base in enumerate(ec_read):
            freq.append([])
            freq[i].append(defaultdict(int))
            freq[i][0][base] += 1
    

        # freq of A C G T and D at each base
        for target_start, target_end, query_name, query_start, query_end, path, strand in cigar_list[target_read]:
            # sr = SeqRecord(reads[query_name].seq)
            # sr.id = reads[query_name].id
            # sr.name = reads[query_name].name
            # sr.description = reads[query_name].description
            # seq_lst.append(sr)
            cnt = 0
            # reverse complement
            # recalculate the query start and end in case of RC
            # rc_start = len - end - 1
            # rc_end = len - start
            target_pos = target_start
            query_pos = query_start
            query_read = reads[query_name].seq

            if strand == '-':
                query_read = query_read.reverse_complement()
                query_pos = len(query_read) - query_end
                query_end = len(query_read) - query_start

            for operation, length in path:
                if operation == '=' or operation == 'X':
                    for i, b in enumerate(query_read[query_pos:query_pos + length]):
                        freq[target_pos+i][0][b] += 1
                    
                    target_pos += length
                    query_pos += length
                elif operation == 'D':
                    for i in range(target_pos, target_pos + length):
                        freq[i][0]['D'] += 1

                    target_pos += length    
                elif operation == 'I':
                    for i, b in enumerate(query_read[query_pos:query_pos + length]):
                        if i == len(freq[target_pos]):
                            freq[target_pos].append(defaultdict(int))
                        freq[target_pos][i][b] += 1

                    query_pos += length
                else:
                    raise ValueError(f'{operation} - Invalid CIGAR operation.')  
            # assert pointers are at the end
            assert target_pos == target_end and query_pos == query_end, f'{target_start}, {target_pos}, {target_end}, {query_start}, {query_pos}, {query_end}, {strand}'
            # assert target_pos == target_end and query_pos == query_end, f'{target_read}, {query_name}'
            # assert target_pos == target_end and query_pos == query_end, f'{path}'
        # generate consensus
        corrected = generate_consensus(ec_read, freq)
        print(len(corrected))
        corrected_seq = SeqRecord(Seq(corrected))
        corrected_seq.id = reads[target_read].id
        corrected_seq.name = reads[target_read].name
        corrected_seq.description = reads[target_read].description
        seq_lst.append(corrected_seq)
    print(len(seq_lst))    
    SeqIO.write(seq_lst, out_path, 'fasta')
    return None

def generate_consensus(ec_read, freq):
    corrected = ''
    # counter = 0
    # dels = 0
    # haps = 0
    insertions, dels = 0, 0
    for pos in range(len(freq)):
        n_support = sum(freq[pos][0].values())
        for i in range(len(freq[pos])):
            if i > 1:
                freq[pos][i]['D'] = n_support - sum(freq[pos][i].values())

            # max_occ = max(freq[pos][i].values(), default=0)
            mc = (Counter(freq[pos][i]).most_common(2))
            if len(mc) == 1:
                b1, c1 = mc[0]
                c2 = 0
            else:
                (b1, c1), (_, c2) = mc[0], mc[1]
            # print(max_occ)
            if c1 == 0:
                break
            if c1 + c2 < 5:
                # #print(cnt)
                # counter += 1
                if i == 0:
                    corrected += ec_read[pos]
            else :
                if c1 <= 1.5 * c2:
                    if i == 0:
                        corrected += ec_read[pos]
                    else:
                        if b1 != 'D':
                            corrected += b1
                else:
                    if b1 != 'D':
                        if i != 0:
                            insertions += 1
                        corrected += b1
                    else:
                        if i == 0:
                            dels += 1

    print("insert" , insertions, dels)                    
    return corrected

def generate_cigar(overlap_map, reads):
    cigar_list = defaultdict(list)
    cnt = 0
    for target_seq in overlap_map: 
        if cnt == 10: break
        # print(cnt, "#######################")
        cnt += 1
        # print(len(overlap_map[query_seq]))
        for overlap in overlap_map[target_seq]:
            (target_start, target_end, query_name, query_start, query_end, path, strand) = calculate_path(overlap, reads, target_seq)
            cigar_list[target_seq].append((target_start, target_end, query_name, query_start, query_end, path, strand))
            # reverse_cigar = (query_start, query_end, target_seq, target_start, target_end, reverse(path), strand)
            # cigar_list[query_name].append(reverse_cigar)

    return cigar_list

def calculate_path(overlap, reads, target_seq):
    target_start = overlap[0]
    target_end = overlap[1]
    target_overlap = reads[target_seq].seq[target_start:target_end]
    # check paf convention of end
    strand = overlap[2]
    query_name = overlap[3]
    query_start = overlap[4]
    query_end = overlap[5]
    if strand == '-':
        query_overlap = reads[query_name].seq[query_start:query_end].reverse_complement()
    else :
        query_overlap = reads[query_name].seq[query_start:query_end]
    path = edlib.align(query_overlap, target_overlap, task = 'path')['cigar']
    generator = gen(path)
    path = list(generator)

    return target_start, target_end, query_name, query_start, query_end, path, strand

REVERSED_OP = {'D': 'I', 'I': 'D', '=': '=', 'X': 'X'}
def reverse(path):
    return [(REVERSED_OP[op], l) for (op, l) in reversed(path)]

def gen(string):
    pattern = re.compile('(\d+)([=XID])')
    for match in pattern.finditer(string):
        yield match.group(2), int(match.group(1))


def get_reads(input_file):
    file_type = Path(input_file).suffix
    file_type = file_type[1:]
    if file_type in ['fq']:
        file_type = 'fastq'

    # records = list(SeqIO.parse(input_file, file_type))
    records_map = SeqIO.to_dict(SeqIO.parse(input_file, file_type))
    return records_map


def parse_paf(paf_file):
    overlap_map = defaultdict(list)
    with open(paf_file, 'r') as f:
        for line in f:
            copy = line.strip()
            (query_seq_name, _, query_start, query_end, strand, overlap_name,
            _, target_start, target_end, _, _, _) = copy.split('\t')
            #   parse to correct data type here
            overlap_map[query_seq_name].append((int(query_start), int(query_end), strand, overlap_name, int(target_start), int(target_end)))
    
    return overlap_map


def parse_args():
    parser = argparse.ArgumentParser(description='Baseline model for error correction')
    parser.add_argument('-i', '--input', type=str, help='path of the input FASTA/Q file.')
    parser.add_argument('-p', '--paf', type=str, help='path of the input PAF file.')
    parser.add_argument('-o', '--output', type=str, help='path of error corrected reads file')
    args = parser.parse_args()
    return args


def main(args):
    reads = get_reads(args.input)
    overlap_map = parse_paf(args.paf)
    cigars = generate_cigar(overlap_map, reads)
    correct_error(cigars, reads, args.output)


if __name__ == '__main__':
   
    args = parse_args()
    
    t1 = time.time()
    main(args)
    t2 = time.time()
    print('Time elapsed: ', t2 - t1)
