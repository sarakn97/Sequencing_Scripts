from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

removals = file = open('/home/sara/stops.txt','r')
content = (removals.read())
removals = content.split(",")
print(removals)

stops = []
for x in removals:
    rec = x.replace("'", "")
    rec = rec.replace("\n", "")
    rec = rec.replace(" ", "")
    rec = rec[0:8]
    stops.append(rec)

print(stops)

#### update to remove stops from fasta
records = SeqIO.parse(str(sys.argv[1]), "fasta")

fasta = []
for record in records:
    RecID = record.id[0:8]
    if RecID not in stops:
        record_final = SeqRecord(
            record.seq,
            id=RecID
        )
        # print(record_final)
        fasta.append(record_final)

# SAVE FILE
SeqIO.write(fasta, str(sys.argv[2]), "fasta")

