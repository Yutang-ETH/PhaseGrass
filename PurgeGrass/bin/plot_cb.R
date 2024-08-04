# visualize cb

source("bin/cb_function.R")

args = commandArgs(trailingOnly=TRUE)

fai <- args[1]

myfai <- read.table(fai, header = F, stringsAsFactors = F)
mycb <- read.table("mymcscanx/mcscanx_cb_position.txt", header = F, stringsAsFactors = F)

mcscanx_linux_cb <- reshape_cb(cb = mycb)

# visualize the overlap
pdf(file = "mymcscanx/mcscanx_linux_cb.pdf", width = 7, height = 5)
for(i in 1:length(unique(mcscanx_linux_cb$contig1))){
  plot_collinear_block(cb = mcscanx_linux_cb[mcscanx_linux_cb$contig1 == unique(mcscanx_linux_cb$contig1)[i], ])
}

dev.off()