# This is a python script for changing headers in a fa file

# python 3

### import modules we need
from Bio import SeqIO
import sys

ID = []
SEQ = []

with open("hbl.fa") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        ID.append(str(record.id))
        SEQ.append(str(record.seq))

with open("hbl.fa", "w") as handle:
    for i in range(len(ID)):
        handle.write(">" + ID[i] + "_left_end" + "\n" + SEQ[i] + "\n")



