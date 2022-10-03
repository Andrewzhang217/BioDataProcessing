#!/usr/bin/env python
import pysam
import argparse





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='clip')
    parser.add_argument('-i', '--input', type=str, help='path of the input bam file.')
    args = parser.parse_args()
    input_sam = args.input
    samfile = pysam.AlignmentFile(input_sam, "r")
    for record in samfile.fetch():
        if record.is_unmapped or record.is_secondary or record.is_supplementary:
            continue
        segment = record.query_alignment_sequence
        name = record.query_name
        print('>' + name)
        print(segment)