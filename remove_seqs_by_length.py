'''
Name: Sara K Nicholson
Title: FASTA FILE PREP FOR PAL2NAL+PAML
Date: February 23, 2024
Description: Remove Sequences that do not meet a length threshold
'''

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

records = SeqIO.parse("/home/sara/PRIYA/mix_transfusiion_R2/D14/TransR2_D14_AA.fa", "fasta")

fasta = []
for record in records:
    if len(record.seq) >= 55:
        record_final = SeqRecord(
            record.seq,
            id=record.id,
            name="ASV",
            description="Mix Transfusion R2 Day 14",
        )
        fasta.append(record_final)

print(len(fasta))

# SAVE FILE
SeqIO.write(fasta, "/home/sara/PRIYA/mix_transfusiion_R2/D14/filt_by_length55.fa", "fasta")






