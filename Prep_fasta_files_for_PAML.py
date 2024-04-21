'''
Name: Sara K Nicholson
Title: FASTA FILE PREP FOR PAL2NAL+PAML
Description: Shorten Names of Sequence Names, Remove any Gaps & Translate Sequences
'''

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# CHANGE SEQ NAMES
seq_records = SeqIO.parse("/home/sara/Marilia/individuals/proteins/MSAs/hypermut_removed/nefW7.fa", "fasta")

fasta = []
for record in seq_records:
    rec = record
    rec = rec[2:]

    RecID = record.id
    record_final = SeqRecord(
        record.seq,
        id=RecID[0:8],
        name="NEF",
        description="WEEK7",
    )
    print(record_final)
    fasta.append(record_final)

# SAVE FILE
SeqIO.write(fasta, "/home/sara/Marilia/individuals/proteins/MSAs/hypermut_removed/Analyzed/nefW7_nc.fa", "fasta")

# REMOVE GAPS FROM SEQS
with open("/home/sara/Marilia/individuals/proteins/MSAs/hypermut_removed/Analyzed/nefW7_nc_ng.fa", "w") as o:
    for record in SeqIO.parse("/home/sara/Marilia/individuals/proteins/MSAs/hypermut_removed/Analyzed/nefW7_nc.fa", "fasta"):
        record.seq = record.seq.replace("-", "")
        SeqIO.write(record, o, "fasta")

# TRANSLATE SEQUENCES
def pad_seq(sequence):
    """ Pad sequence to multiple of 3 with N """

    remainder = len(sequence) % 3

    return sequence if remainder == 0 else sequence + Seq('N' * (3 - remainder))


seq_records = SeqIO.parse("/home/sara/Marilia/individuals/proteins/MSAs/hypermut_removed/Analyzed/nefW7_nc_ng.fa", "fasta")

fasta_aa = []
# ids = []

for record in seq_records:
    #rec = record.seq.reverse_complement()
    rec = record.seq
    rec = rec[0:]
    #rec1 = pad_seq(rec).translate()
    RecID = record.id
    # ids.append(RecID)
    record_final = SeqRecord(
        pad_seq(rec).translate(),
        id=RecID,
        name="NEF",
        description="w7",
    )
    print(record_final)
    fasta_aa.append(record_final)

# WRITE TRANSLATED SEQS TO FILE
SeqIO.write(fasta_aa, "/home/sara/Marilia/individuals/proteins/MSAs/hypermut_removed/Analyzed/nefW7_AA.fa", "fasta")