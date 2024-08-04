# produce trim files

source("bin/cb_function.R")

args = commandArgs(trailingOnly=TRUE)

fai <- args[1]

myfai <- read.table(fai, header = F, stringsAsFactors = F)
mycb <- read.table("mymcscanx/mcscanx_cb_position.txt", header = F, stringsAsFactors = F)
mycoords <- read.table("ps.coords", header = F, stringsAsFactors = F)

Coords_and_haplotype_boundary_bed(cbp = mycb, myfai = myfai, coords = mycoords)