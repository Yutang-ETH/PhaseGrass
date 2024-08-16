#!/usr/bin/env python

# 25.12.2021
# I have H1, H2 and untagged reads, I want to find from all reads the remainning reads excluding H1, H2 and untagged reads
# let's call the remainning reads as unmapped (unmapped to chrs)

H1 = []
with open("whatshap_phased_read_list/H1.list", 'rb') as f:
    for read in f:
        H1.append(read)

H2 = []
with open("whatshap_phased_read_list/H2.list", 'rb') as f:
    for read in f:
        H2.append(read)
 
untagged = []
with open("whatshap_phased_read_list/none.list", 'rb') as f:
    for read in f:
        untagged.append(read)

phased = H1 + H2 + untagged

All = []
with open("whatshap_phased_read_list/all.list", 'rb') as f:
    for read in f:
        All.append(read)

unmapped = list(set(All) - set(phased))

with open("whatshap_phased_read_list/unmapped.list", "wb") as f:
    for read in unmapped:
        f.write(read)

    
