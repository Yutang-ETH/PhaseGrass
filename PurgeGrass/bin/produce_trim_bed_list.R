# produce trim files

source("bin/cb_function.R")

args = commandArgs(trailingOnly=TRUE)

mycoords <- read.table("ps.coords", header = F, stringsAsFactors = F)

Coords_and_haplotype_boundary_bed(coords = mycoords)