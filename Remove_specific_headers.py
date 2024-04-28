#from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord

#seq_records = SeqIO.parse("/home/sara/Marilia/individuals/proteins_24/gag_allseqs.fa", "fasta")
removals = file = open('/home/sara/gag_stops.txt','r')
content = (removals.read())
print(content)
l = content.split(",")
print(l)

