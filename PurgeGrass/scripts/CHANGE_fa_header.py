# This is a python script for changing headers in a fa file

# python 3

### import modules we need
from Bio import SeqIO
import sys

ID = []
SEQ = []

with open("gmap_protein.fasta") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        ID.append(str(record.id))
        SEQ.append(str(record.seq))

oldnames = []
newnames = []
with open("names") as handle:
    for name in handle:
        oldnames.append(name.split()[0])
        newnames.append(name.split()[1])

res = dict(zip(oldnames, newnames))

names = []
for name in ID:
    names.append(res[name])


with open("gmap_protein_A.fa", "w") as handle:
    for i in range(len(ID)):
        handle.write(">A_" + names[i] + "\n" + SEQ[i] + "\n")

with open("gmap_protein_B.fa", "w") as handle:
    for i in range(len(ID)):
        handle.write(">B_" + names[i] + "\n" + SEQ[i] + "\n")

