#!/usr/bin/env python

"""
Name : Sara Nicholson
Date : March 2024
Description: Counts kmers for each fasta file in a directory
"""
from pathlib import Path
import sys
import re
import os


def count_kmer(fastq_file, kmer_length, outfile):
    # initialize variable to place sequence
    seq = ''
    # Initialize kmer dictionary
    kmer_dictionary = {}

    # open file & add kmers from each sequence (seq) to kmer dictionary
    with open(fastq_file, 'r') as fastq_seq:
        for line in fastq_seq:
            line = line.rstrip()
            if re.match('^[ATGCN]+$', line):
                seq += line
                # set stop to length of sequence
                stop = len(seq) - kmer_length + 1
                for start in range(0, stop):
                    kmer = seq[start:start + kmer_length]

                    if kmer in kmer_dictionary:
                        # If kmer is already in dictionary, add 1 to count
                        kmer_dictionary[kmer] += 1
                    else:
                        # If kmer is not in dictionary, add the kmer to dictionary with count = 1
                        kmer_dictionary[kmer] = 1
            # clear the sequence variable to add next sequence
            seq = ''

    # Set variable for tab seperation
    sep = "\t"

    # Open output file for writing
    with open(outfile, 'w') as out:
        for kmer in kmer_dictionary:
            count = kmer_dictionary[kmer]
            output = (kmer, str(count))
            # format output to kmer + tab + count + newline
            out.write(sep.join(output) + "\n")


if __name__ == "__main__":
    # Initializations
    fasta = str(sys.argv[1])
    k_size = int(sys.argv[2])
    # outfile = str(sys.argv[3])

    for filename in os.listdir(fasta):
        f = os.path.join(fasta, filename)
        outfile = Path(f).stem
        count_kmer(f, k_size, outfile)
