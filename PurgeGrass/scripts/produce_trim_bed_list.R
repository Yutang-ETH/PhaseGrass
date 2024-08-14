# produce trim files

args = commandArgs(trailingOnly=TRUE)

script <- args[1]

source(paste(script, "/", "cb_function.R", sep=""))

mycoords <- read.table("ps.coords", header = F, stringsAsFactors = F)

Coords_and_haplotype_boundary_bed(coords = mycoords)