
args = commandArgs(trailingOnly=TRUE)

myinput <- args[1]
myoutput <- args[2]

ct <- read.table(myinput, header = F, stringsAsFactors = F)
primary <- rep(ct[1,1], length(ct[-1, 1]))
ctf <- cbind.data.frame(primary, ct[-1, ])
# str(ctf)
# names(ctf) <- c("primary", "associated", "tag")


write.table(ctf, myoutput, quote = F, row.names = F, col.names = F, sep = "\t")
