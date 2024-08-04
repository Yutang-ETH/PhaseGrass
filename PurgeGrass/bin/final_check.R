# check with busco table and output final primary and allelic contigs

args = commandArgs(trailingOnly=TRUE)

# args[1] = fai
# args[2] = secondary_contigs.txt
# args[3] = primary_contigs.txt
# args[4] = purge_haplotigs.txt
# args[5] = full_table.tsv


fai <- read.table(args[1], header = F, stringsAsFactors = F)
secondary_contig <- read.table(args[2], header = F, stringsAsFactors = F)
primary_contig <- read.table(args[3], header = F, stringsAsFactors = F)
purge_hap <- read.table(args[4], header = F, stringsAsFactors = F)
bus <- read.table(args[5], header = F, stringsAsFactors = F, sep = "\t", fill = T, quote = "\"")

# check how many duplicated busco genes in haplotig
bus_dup <- bus[bus$V2 == "Duplicated", ]

# allelic contigs to be removed
hp <- union(purge_hap$V1, secondary_contig$V1)

# check how many duplicated busco genes are left
bus_left <- bus[!bus$V3 %in% hp, ]

bus_dup_remain <- NULL
for(x in unique(bus_left$V1)){
  xx <- bus_left[bus_left$V1 == x, ]
  if (length(unique(xx$V3)) > 1){
    bus_dup_remain <- rbind(bus_dup_remain, xx)
  }
}

# number of remaining duplicated busco
print(paste("Number of remaining duplicated BUSCOs after removing allelic contigs",
            length(unique(bus_dup_remain$V1)),
            sep = ": "))

sprintf("Percent of remaining duplicated BUSCOs after removing allelic contigs: %0.2f%%", 
        length(unique(bus_dup_remain$V1))/length(unique(bus$V1)) * 100)

# number of remaining busco, excluding missing ones
print(paste("Number of total remaining BUSCOs after removing allelic contigs",
            length(unique(bus_left[bus_left$V2 != "Missing", 1])),
            sep = ": "))

sprintf("Percent of total remaining BUSCOs after removing allelic contigs: %0.2f%%", 
        length(unique(bus_left[bus_left$V2 != "Missing", 1]))/length(unique(bus$V1)) * 100)


# final haplotigs
final_haplotigs <- hp
print(paste("Size of final allelic contigs", sum(fai[fai$V1 %in% final_haplotigs, 2]), sep = ": ")) 

# final primary contigs
final_primary_contigs <- setdiff(fai$V1, final_haplotigs)
print(paste("Size of final primary contigs", sum(fai[fai$V1 %in% final_primary_contigs, 2]), sep = ": "))

# write out final haplotigs and primary contigs
write.table(final_haplotigs, "final_suspect_haplotigs.txt", quote = F, row.names = F, col.names = F, sep = "\t")
write.table(final_primary_contigs, "final_primary_contigs.txt", quote = F, row.names = F, col.names = F, sep = "\t")

