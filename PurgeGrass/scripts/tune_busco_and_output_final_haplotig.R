

args = commandArgs(trailingOnly=TRUE)


# import data
myfai <- read.table(args[1], header = F, stringsAsFactors = F)
secondary <- read.table(args[2], header = F, stringsAsFactors = F)
primary <- read.table(args[3], header = F, stringsAsFactors = F)
haplotig <- read.table(args[4], header = F, stringsAsFactors = F)
bus <- read.table(args[5], header = F, stringsAsFactors = F, sep = "\t", fill = T, quote = "\"")
purge_haplotig_primary <- read.table(args[6], header = F, stringsAsFactors = F)

script <- args[7]
source(paste(script, "/", "cb_function.R", sep=""))

remove_redundant_busco_new(bus = bus, haplotig = haplotig, purge_haplotig_primary = purge_haplotig_primary, secondary = secondary, 
                       primary = primary, myfai = myfai)