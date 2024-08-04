# get common haplotigs between purge haplotigs and canu

source("bin/cb_function.R")

args = commandArgs(trailingOnly=TRUE)

repeat_list <- read.table(args[1], header = F, stringsAsFactors = F)
haplotig_list <- read.table(args[2], header = F, stringsAsFactors = F)

find_common_haplotigs(repeat_list = repeat_list, haplotig = haplotig_list)


