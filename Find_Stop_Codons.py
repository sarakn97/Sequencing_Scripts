"""
Description: Find Sequences with Stop Codons in Open-Reading Frame in Fasta File,
            Run on command line: "python3 Find_Stop_Codons.py /path/to/fasta.fa"
             returns headers of fasta files with Stop Codons in ORF
Name: Sara Nicholson
Date: 27 April 2024
"""
import numpy
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import re
import os


def findStopCodons(orf):
    catch = numpy.arange(0, len(orf), 3)
    stopCodon = []
    for i in catch:
        codon = orf[i:i + 3]
        if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
            stopCodon.append(i + 1)
    return stopCodon


records = SeqIO.parse(str(sys.argv[1]), "fasta")

stops = []
for record in records:
    if findStopCodons(record.seq):
        stops.append(record.id)

print(stops)

#### update to remove stops from fasta
records = SeqIO.parse(str(sys.argv[2]), "fasta")

fasta = []
for record in records:
    RecID = record.id
    if RecID not in stops:
        record_final = SeqRecord(
            record.seq,
            id=RecID
        )
        print(record_final)
        fasta.append(record_final)

# SAVE FILE
SeqIO.write(fasta, str(sys.argv[3]), "fasta")
