from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import edlib
from pathlib import Path
from collections import defaultdict
import argparse


def correct_error(cigar_list, reads, out_path):
    seq_lst = []
    for target_read in cigar_list:
        
        #print(len(cigar_list[query_read]))
        ec_read = reads[target_read].seq
        seq_lst.append(reads[target_read])
        print(len(ec_read))
        freq = []
        for i in range(len(ec_read)):
            freq.append([])
            freq[i].append(defaultdict(int))
        # freq of A C G T at each base
        for target_start, target_end, query_name, query_start, query_end, path in cigar_list[target_read]:
            # sr = SeqRecord(reads[query_name].seq)
            # sr.id = reads[query_name].id
            # sr.name = reads[query_name].name
            # sr.description = reads[query_name].description
            # seq_lst.append(sr)
            cnt = 0
            target_pos = target_start
            query_pos = query_start
            length = 0
            for i in range(len(path)):
                if path[i].isdigit():
                    length = length * 10 + int (path[i])
                else:
                    if path[i] == '=':
                        for r in range(target_pos, target_pos + length):
                            original_base = ec_read[r]
                            freq[r][0][original_base] += 1
                        target_pos += length
                        query_pos += length
                    elif path[i] == 'X':
                        for r in range(target_pos, target_pos + length):
                            mismatch_base = reads[query_name].seq[query_pos]
                            freq[r][0][mismatch_base] += 1
                        target_pos += length
                        query_pos += length
                    elif path[i] == 'D':
                        for r in range(target_pos, target_pos + length):
                            freq[r][0]['N'] += 1
                        target_pos += length    
                    elif path[i] == 'I':
                        insertion_bases = reads[query_name].seq[query_pos:query_pos + length]
                        for index in range(len(insertion_bases)):
                            if index < len(freq[target_pos]):
                                freq[target_pos][index][insertion_bases[index]] += 1
                            else :
                                freq[target_pos].append(defaultdict(int))
                                freq[target_pos][index][insertion_bases[index]] += 1

                        query_pos += length    
                    length = 0

        # generate consensus
        corrected = generate_consensus(ec_read, freq)
        print(len(corrected))
        #print(ec_read)
        # cc = edlib.align(ec_read, corrected)
        # print(cc)
        # corrected_seq = SeqRecord(Seq(corrected))
        # corrected_seq.id = reads[target_read].id
        # corrected_seq.name = reads[target_read].name
        # corrected_seq.description = reads[target_read].description
        # seq_lst.append(corrected_seq)
        #seq_lst.append(SeqRecord(Seq(ec_read)))
    print(len(seq_lst))    
    SeqIO.write(seq_lst, out_path, 'fasta')
    return None

def generate_consensus(ec_read, freq):
    corrected = ''
    cnt = 0
    counter = 0
    dels = 0
    haps = 0
    for pos in range(len(freq)):
        for i in range(len(freq[pos])):
            max_occ = max(freq[pos][i].values(), default=0)
            # print(max_occ)
            if max_occ < 5:
                # #print(cnt)
                # cnt += 1
                counter += 1
                corrected += ec_read[pos]
            else :
                max2_occ = 0
                for v in freq[pos][i].values():
                    if (v > max2_occ and v < max_occ):
                        max2_occ = v
                if max_occ <= 2 * max2_occ:
                    haps += 1
                    corrected += ec_read[pos]
                    continue        
                majority = max (freq[pos][i], key = freq[pos][i].get)
                if majority == 'N':
                    dels += 1
                    continue
                else:
                    # print(cnt)
                    # cnt += 1
                    corrected += majority
    # print("countr", counter)
    # print(dels)
    # print(haps)
    return corrected

def generate_cigar(align_map, reads):
    cigar_list = defaultdict(list)
    cnt = 0
    for target_seq in align_map: 
        if cnt == 1: break
        # print(cnt, "#######################")
        cnt += 1
        # print(len(align_map[query_seq]))
        for overlap in align_map[target_seq]:
            target_start = int (overlap[0])
            target_end = int (overlap[1])
            target_overlap = reads[target_seq].seq[target_start:target_end + 1]
            strand = reads[target_seq].seq[2]
            query_name = overlap[3]
            query_start = int (overlap[4])
            query_end = int (overlap[5])
            if strand == '-':
                query_overlap = reads[query_name].seq[query_start:query_end + 1].reverse_complement()
            else :
                query_overlap = reads[query_name].seq[query_start:query_end + 1]
            path = edlib.align(query_overlap, target_overlap, task = 'path')['cigar']
            cigar = (target_start, target_end, query_name, query_start, query_end, path)
            cigar_list[target_seq].append(cigar)
    return cigar_list


def get_reads(input_file):
    file_type = Path(input_file).suffix
    file_type = file_type[1:]
    if file_type in ['fq']:
        file_type = 'fastq'

    records = list(SeqIO.parse(input_file, file_type))
    records_map = defaultdict(None)
    for record in records:
        records_map[record.name] = record
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
        dst[query_seq_name].append((query_start, query_end, strand, overlap_name, target_start, target_end))
    return dst


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Baseline model for error correction')
    parser.add_argument('-i', '--input', type=str, help='path of the input FASTA/Q file.')
    parser.add_argument('-p', '--paf', type=str, help='path of the input PAF file.')
    parser.add_argument('-o', '--output', type=str, help='path of error corrected reads file')
    args = parser.parse_args()

    # print("before")
    reads = get_reads(args.input)
    align_map = parse_paf(args.paf)
    cigars = generate_cigar(align_map, reads)
    correct_error(cigars, reads, args.output)
    # print('Finished')
