
################################################################################################################################

# reshape the data table, put longer contig first
reshape_cb <- function(cb = newcb){
  
  names(cb) <- c("contig1", "start1", "end1", "contig2", "start2", "end2")
  
  # reshpae the data table
  newcb_reshaped <- cb[NULL, ]
  for(i in 1:dim(cb)[1]){
    if(myfai[myfai$V1 == cb$contig1[i], 2] >= myfai[myfai$V1 == cb$contig2[i], 2]){
      xx <- cb[i, ]
      names(xx) <- colnames(cb)
      newcb_reshaped <- rbind.data.frame(newcb_reshaped, xx)
    }else{
      xx <- cbind.data.frame(cb[i, 4:6], cb[i, 1:3])
      names(xx) <- colnames(cb)
      newcb_reshaped <- rbind.data.frame(newcb_reshaped, xx)
    }
  }
  
  # remove dplicated pairs
  newcb_reshaped$ID <- paste(newcb_reshaped$contig1, newcb_reshaped$start1, newcb_reshaped$contig2, newcb_reshaped$start2,sep = "&")
  # print(newcb_reshaped[duplicated(newcb_reshaped$ID), ])
  Final_cb <- newcb_reshaped[!duplicated(newcb_reshaped$ID), ]
  
  # return(newcb_reshaped)
  Final_cb$ID <- paste(Final_cb$contig1, Final_cb$contig2, sep = "&")
  Final_cb <- Final_cb[Final_cb$contig1 != Final_cb$contig2, ]
  
  # write out primary contigs and secondary contigs
  primary_contigs <- unique(Final_cb$contig1)
  secondary_contigs <- unique(Final_cb$contig2)
  secondary_contigs <- setdiff(secondary_contigs, intersect(primary_contigs, secondary_contigs))
  write.table(primary_contigs, "primary_contigs.txt", quote = F, row.names = F, col.names = F)
  write.table(secondary_contigs, "secondary_contigs.txt", quote = F, row.names = F, col.names = F)
  
  # write out the ID too
  ID <- unique(Final_cb$ID)
  write.table(ID, "collinear_block_ID.txt", quote = F, row.names = F, col.names = F)
  
  return(Final_cb)
}

############################################################################################################################

plot_collinear_block <- function(cb = cb, fai = myfai, xbase = 500, yinter = 200){
  
  # sort cb based on start1
  cb <- cb[order(cb$start1), ]
  
  # what haplotigs are there
  contig2_uniq <- unique(cb$contig2)
  print(contig2_uniq)
  
  # find the length of the primary contig
  contig1_len <- fai[fai$V1 == unique(cb$contig1), 2]
  print(contig1_len)
  
  # find the length of every haplotigs
  contig2_len <- c()
  for(i in 1:length(contig2_uniq)){
    contig2_len <- c(contig2_len, fai[fai$V1 == contig2_uniq[i], 2])
  }
  print(contig2_len)
  
  ybase <- round((round(contig1_len/10^3, 1) + 2*xbase)/5, 1)
  
  # plot the primary contig first
  par(mar = c(10, 5, 4, 4))
  plot(x = 1:(round(contig1_len/10^3, 1) + 2*xbase), y = 1:(round(contig1_len/10^3, 1) + 2*xbase), type = "n",
      axes = F, xlab = unique(cb$contig1), ylab = " ")
  segments(x0 = xbase, y0 = ybase, x1 = xbase + round(contig1_len/10^3, 1), y1 = ybase)
  axis(1, at = seq(xbase, xbase + round(contig1_len/10^3, 1), 100), labels = seq(0, round(contig1_len/10^3, 1), 100))
  mtext("kb", side = 1, at = xbase, las = 1, line = -1.5)
  
  mycol <- c("lightblue3", "indianred2")
  mycolxx <- c()
  for(i in 1:length(contig2_uniq)){
    mycolxx <- c(mycolxx, mycol)
  }
  
  # plot every haplotigs
  for(i in 1:length(contig2_uniq)){
    xx <- cb[cb$contig2 == contig2_uniq[i], ]
    if(xx$start2[1] < xx$start2[2]){
      xx_start <- xbase + round(xx$start1[1]/10^3, 1) - round(xx$start2[1]/10^3, 1)
      xx_end <- xbase + round(tail(xx$start1, 1)/10^3, 1) + round(contig2_len[i]/10^3, 1) - round(tail(xx$start2, 1)/10^3, 1)
      segments(x0 = xx_start, y0 = ybase + yinter*i, x1 = xx_end, y1 = ybase + yinter*i)
      text(x = mean(xbase + round(head(xx$start1, 1)/10^3, 1), xbase + round(tail(xx$start1, 1)/10^3, 1)),
           y = 1.2*ybase + yinter*i, 
           labels = contig2_uniq[i], 
           cex = 0.7)
      for(j in 1:dim(xx)[1]){
        segments(x0 = xbase + round(xx$start1[j]/10^3, 1), y0 = ybase, 
                 x1 = xbase + round(xx$start1[j]/10^3, 1), y1 = ybase + yinter*i,
                 col = mycolxx[i])
      }
    }else{
      xx_start <- xbase + round(xx$start1[1]/10^3, 1) - (round(contig2_len[i]/10^3, 1) - round(head(xx$start2, 1)/10^3, 1))
      xx_end <- xbase + round(tail(xx$start1, 1)/10^3, 1) + round(tail(xx$start2, 1)/10^3, 1)
      segments(x0 = xx_start, y0 = ybase + yinter*i, x1 = xx_end, y1 = ybase + yinter*i)
      text(x = mean(xbase + round(head(xx$start1, 1)/10^3, 1), xbase + round(tail(xx$start1, 1)/10^3, 1)),
           y = 1.2*ybase + yinter*i, 
           labels = contig2_uniq[i], 
           cex = 0.7)
      for(j in 1:dim(xx)[1]){
        segments(x0 = xbase + round(xx$start1[j]/10^3, 1), y0 = ybase, 
                 x1 = xbase + round(xx$start1[j]/10^3, 1), y1 = ybase + yinter*i,
                 col = mycolxx[i])
      }
    }
  }
}


#############################################################################################################################

# process coords and define haplotype boundary, use this one to produce a bed file and then use bedtool getfasta to extract sequences

Coords_and_haplotype_boundary_bed <- function(coords = coords){
  
  # reshape the coords file, put large contig first and shorter contig next
  
  if(unique(coords[, 8] >= coords[, 9])){
    
    coords_reshape <- coords
    colnames(coords_reshape) <- c("ps", "pe", "hs", "he", "pal", "hal", "identity", "pl", "hl", "pac", "hac", "primary", "haplotig") # al means alignment length, ac means alignment coverage
    
    
  }else{
    
    xx <- coords[coords[, 8] >= coords[, 9], ]
    colnames(xx) <- c("ps", "pe", "hs", "he", "pal", "hal", "identity", "pl", "hl", "pac", "hac", "primary", "haplotig")
    
    yy <- coords[coords[, 8] < coords[, 9], ]
    yy <- cbind.data.frame(yy[, 3:4], yy[, 1:2], yy[, 6], yy[, 5], yy[, 7], yy[, 9], yy[, 8], yy[, 11], yy[, 10], yy[, 13], yy[, 12])
    colnames(yy) <- c("ps", "pe", "hs", "he", "pal", "hal", "identity", "pl", "hl", "pac", "hac", "primary", "haplotig")
    
    coords_reshape <- rbind.data.frame(xx, yy)
    
  }
  
  
  
  # calculate how long the collinear block spans the suspect haplotig
  # should calculate the alignment coverage based on the coords file
  
  coords_reshape$tag <- paste(coords_reshape$primary, coords_reshape$haplotig, sep = "&")
  hcov <- c()
  ascore <- c() # alingment score
  for(x in unique(coords_reshape$tag)){
    xx <- coords_reshape[coords_reshape$tag == x, ]
    y <- (max(c(xx$hs, xx$he)) - min(c(xx$hs, xx$he)))/xx$hl
    z <- rep(sum(xx$hac), length(y))
    hcov <- c(hcov, y)
    ascore <- c(ascore, z)
  }
  coords_reshape$hcov <- hcov
  coords_reshape$ascore <- ascore
  
  # calculate how long is left unaligned
  coords_reshape$hcovl <- coords_reshape$hl*(1-coords_reshape$hcov)
  
  # write out the reshaped coords table
  write.table(coords_reshape, "coords_reshape.txt", quote = F, sep = "\t", col.names = T, row.names = F)
  
  # only select haplotigs with hcov <= 80% and hcovl >= 500 Kb for trimming
  coords_trim <- coords_reshape[coords_reshape$hcov <= 0.8 & coords_reshape$hcovl >= 500000, ]
  
  # trim the larger end, the longer unaligned end of haplotig, add the longer unaligned end to the final assembly
  hbl_hbr <- NULL
  for(x in unique(coords_trim$tag)){
    xx <- coords_trim[coords_trim$tag == x, ]
    if(min(c(xx$hs, xx$he)) >= (unique(xx$hl) - max(c(xx$hs, xx$he)))){
      
      # in this case, left end is longer than the right end, trim the left end 
      y <- cbind.data.frame(unique(xx$haplotig), 0, (min(c(xx$hs, xx$he))-1), stringsAsFactors = FALSE)
      colnames(y) <- c("haplotig", "start", "end")
      hbl_hbr <- rbind.data.frame(hbl_hbr, y)
      
    }else{
      
      # in this case, right end is longer than the left end, trim the right end
      y <- cbind.data.frame(unique(xx$haplotig), (max(c(xx$hs, xx$he))-1), unique((xx$hl-1)), stringsAsFactors = FALSE)
      colnames(y) <- c("haplotig", "start", "end")
      hbl_hbr <- rbind.data.frame(hbl_hbr, y)
      
    }
    
  }
  
  # now check if the trimmed contig is in both primary and haplotig list
  # if so, exclude the contig, because it might be a chimera
  # if not, keep it and output the hbl_hbr file for trimming
  # also check if there are any duplications in hbl_hbr
  
  hbl_hbr$sed <- hbl_hbr$end - hbl_hbr$start # sed means difference between start and end
  hbl_hbr <- hbl_hbr[order(hbl_hbr$sed, decreasing = T), ]
  hbl_hbr <- hbl_hbr[!duplicated(hbl_hbr$haplotig), ]
  hbl_hbr <- hbl_hbr[hbl_hbr$haplotig %in% setdiff(hbl_hbr$haplotig, coords_reshape$primary), ]
  
  # write out the trim list
  write.table(hbl_hbr[, 1:3], "hbl_hbr.bed", quote = F, row.names = F, col.names = F, sep = "\t")
}


#############################################################################################################################

# process coords and define haplotype boundary

Coords_and_haplotype_boundary <- function(cbp = cbp, myfai = myfai, coords = coords){
  
  # put larger contig first and shorter contig next
  
  large <- NULL
  small <- NULL
  for(i in 1:dim(cbp)[1]){
    if(myfai[myfai$V1 == cbp[i, 1], 2] >= myfai[myfai$V1 == cbp[i, 4], 2]){
      xx <- cbind.data.frame(cbp[i, 1:6], myfai[myfai$V1 == cbp[i, 1], 2], myfai[myfai$V1 == cbp[i, 4], 2])
      # print(xx)
      large <- rbind.data.frame(large, xx)
    }else{
      xx <- cbind.data.frame(cbp[i, 4:6], cbp[i, 1:3], myfai[myfai$V1 == cbp[i, 4], 2], myfai[myfai$V1 == cbp[i, 1], 2])
      # print(xx)
      small <- rbind.data.frame(small, xx)
    }
  }
  
  
  names(small) <- c("primary", "ps", "pe", "haplotig", "hs", "he", "pl", "hl")
  names(large) <- c("primary", "ps", "pe", "haplotig", "hs", "he", "pl", "hl")
  
  unique(large$pl > large$hl)
  unique(small$pl > small$hl)
  
  cbp <- rbind.data.frame(large, small)
  unique(cbp$pl > cbp$hl)
  rm(large, small, xx)
  
  # calculate how long the collinear block spans the suspect haplotig
  cbp$tag <- paste(cbp$primary, cbp$haplotig, sep = "&")
  hcov <- c()
  for(x in unique(cbp$tag)){
    xx <- cbp[cbp$tag == x, ]
    y <- (max(xx$he) - min(xx$hs))/unique(xx$hl)
    hcov <- c(hcov, rep(y, dim(xx)[1]))
  }
  cbp$hcov <- hcov
  
  # export haplotigs with small overlaps and trim the overlap
  length(unique(cbp[cbp$hcov <= 1, 10]))  
  
  # find collinear blocks in coords
  coords$V14 <- paste(coords$V12, coords$V13, sep = "&")
  coords_cb <- coords[coords$V14 %in% cbp[cbp$hcov <= 1, 9], ]
  
  # calculate the alignment coverage of the short contig
  aba_cov <- c()
  for(x in unique(coords_cb$V14)){
    xx <- coords_cb[coords_cb$V14 == x, ]
    if(xx[1, 3] >= xx[1, 4]){
      y <- (max(xx$V3) - min(xx$V4))/unique(xx$V9)
    }else{
      y <- (max(xx$V4) - min(xx$V3))/unique(xx$V9)
    }
    aba_cov <- c(aba_cov, rep(y, dim(xx)[1]))
  }
  coords_cb$V15 <- aba_cov
  
  # define the haplotig boundary
  output <- NULL
  for(x in unique(coords_cb$V14)){
    xx <- coords_cb[coords_cb$V14 == x, ]
    if(xx[1, 3] >= xx[1, 4]){
      y <- (max(xx$V3) - min(xx$V4))/unique(xx$V9)
      hs <- min(xx$V4)
      he <- max(xx$V3)
      dhe <- unique(xx$V9) - he
    }else{
      y <- (max(xx$V4) - min(xx$V3))/unique(xx$V9)
      hs <- min(xx$V3)
      he <- max(xx$V4)
      dhe <- unique(xx$V9) - he
    }
    output <- rbind(output, c(unique(xx$V13), unique(xx$V12), hs, he, dhe, unique(xx$V9), y))
  }
  
  hb <- type.convert(as.data.frame(output), as.is = T)
  names(hb) <- c("haplotig", "primary", "hs", "he", "dhe", "hl", "cov")
  rm(output)
  
  
  # calculate overall alignment coverage as some haplotigs aligns to many primary contigs
  allcov <- c()
  for(i in 1:dim(hb)[1]){
    
    xx <- hb[hb$haplotig == hb$haplotig[i], ]
    allcov <- c(allcov, sum(xx$cov))
    
  }
  
  hb$allcov <- allcov
  
  # select haplotigs with more than 500kb not overlapping with primary contigs and the all coverage lower than 0.8
  # first lb >= 300kb. For these contigs, trim from right to keep the left end
  hbl <- hb[hb$hs >= 500000 & hb$allcov < 0.8, c(1, 3, 6)]
  hbl$trim <- hbl$hl - hbl$hs
  hbl <- hbl[, c(1, 4)]
  hbl <- hbl[order(hbl$trim, decreasing = T), ]
  if(length(hbl[duplicated(hbl$haplotig), 1]) > 0){
    hbl <- hbl[!duplicated(hbl$haplotig), ]
  }else{
    hbl <- hbl
  }
  
  # second rb >= 500kb and allcoverage < 0.8. For these contigs, trim from left to keep the right end
  hbr <- hb[hb$dhe >= 500000 & hb$allcov < 0.8, c(1, 4)]
  hbr <- hbr[order(hbr$he, decreasing = T), ]
  if(length(hbr[duplicated(hbr$haplotig), 1]) > 0){
    hbr <- hbr[!duplicated(hbr$haplotig), ]
  }else{
    hbr <- hbr
  }
  
  write.table(hbl, "hbl.txt", quote = F, row.names = F, col.names = F, sep = "\t")
  write.table(hbr, "hbr.txt", quote = F, row.names = F, col.names = F, sep = "\t")
}


#############################################################################################################################

find_common_haplotigs <- function(repeat_list = repeat_list, haplotig = haplotig){
  
  
  suspect_haplotig <- setdiff(haplotig[, 2], repeat_list[, 1])
  
  write.table(suspect_haplotig, "suspect_haplotigs.txt", quote = F, row.names = F, col.names = F, sep = "\t")
  
}


#############################################################################################################################

remove_redundant_busco <- function(bus = bus, haplotig = haplotig, suspect = suspect, secondary = secondary, 
                                   primary = primary, myfai = myfai){
  
  # check how many duplicated busco genes in haplotig
  length(unique(bus[bus$V2 == "Duplicated", 1])) # 1022 duplicated genes
  bus_dup <- bus[bus$V2 == "Duplicated", ]
  bus_single <- bus[bus$V2 == "Complete", ]
  
  # check whether duplicates are all having two copies
  test <- c()
  for (i in 1: length(unique(bus_dup$V1))){
    xx <- bus_dup[bus_dup$V1 == unique(bus_dup$V1)[i], ]
    test[i] <- dim(xx)[1]
  }
  
  unique(test) # 2 3 4 5, there are some busco genes having more than 2 copies
  
  # gene paris, some genes have more than two copies
  gp <- unique(bus_dup$V1)[which(test == 2)]
  bus_gp <- bus_dup[bus_dup$V1 %in% gp, ]
  
  # 
  y <- c()
  for(x in unique(haplotig$V2)){
    if(x %in% bus_dup$V3){
      y <- c(y, length(unique(bus_dup[bus_dup$V3 == x, 1])))
    }else{
      y <- c(y, 0)
    }
  }
  haplotig$V7 <- y
  
  # suspect haplotigs from haplotig and cbp
  hp <- unique(c(suspect$V1, secondary$V1))
  
  # check how many duplicated busco genes are left
  bus_left <- bus_gp[!bus_gp$V3 %in% hp, ]
  
  hap_bus <- c()
  for(x in unique(bus_left$V1)){
    xx <- bus_left[bus_left$V1 == x, ]
    if (length(unique(xx$V3)) > 1){
      yy <- myfai[myfai$V1 %in% unique(xx$V3), ]
      yy <- yy[order(yy$V2, decreasing = T), ]
      hap_bus <- c(hap_bus, yy[-1, 1])
    }
  }
  
  hap_bus <- myfai[myfai$V1 %in% hap_bus, ]
  hap_bus <- hap_bus[!hap_bus$V1 %in% unique(c(primary$V1, suspect$V1)), ]
  sum(hap_bus$V2) # 142257687
  
  sum(myfai[myfai$V1 %in% unique(c(hap_bus$V2, suspect$V1, secondary$V1)), 2])
  
  bus_single_final_haplotigs <- intersect(bus_single$V3, unique(c(hap_bus$V1, suspect$V1, secondary$V1)))
  bus_rescue <- intersect(bus_single_final_haplotigs, setdiff(bus_single_final_haplotigs, c(hap_bus$V1, secondary$V1)))
  
  final_suspect_haplotigs <- setdiff(c(hap_bus$V1, suspect$V1, secondary$V1), bus_rescue)
  
  write.table(final_suspect_haplotigs, "final_suspect_haplotigs.txt", quote = F, row.names = F, col.names = F, sep = "\t")
  write.table(setdiff(myfai$V1, final_suspect_haplotigs), "final_primary_contigs.txt", quote = F, row.names = F, col.names = F, sep = "\t")
}

#############################################################################################################################
parse_coords <- function(myfai = myfai, mycoords = mycoords){
  
  library(foreach)
  library(doParallel)
  
  # myfai <- read.table(args[1], header = F, stringsAsFactors = F)
  
  # mycoords <- read.table(args[2], header = F, stringsAsFactors = F)
  
  # remove self alignment, V12 is ref, V13 is query
  mycoords_f <- mycoords[mycoords$V12 != mycoords$V13, ]
  
  # add tags to every alignment
  mycoords_f$V14 <- paste(mycoords_f$V12, mycoords_f$V13, sep = "&")
  
  registerDoParallel(48)
  
  newcoords <- foreach(x=unique(mycoords_f$V14), .combine=rbind) %dopar% {
    
    xx <- mycoords_f[mycoords_f$V14 == x, ]
    
    col1 <- xx[1, 1]
    col2 <- xx[dim(xx)[1], 2]
    col3 <- xx[1, 3]
    col4 <- xx[dim(xx)[1], 4]
    col5 <- sum(xx[, 5])
    col6 <- sum(xx[, 6])
    col7 <- round(mean(xx[, 7]), 2)
    col8 <- unique(xx[, 8])
    col9 <- unique(xx[, 9])
    col10 <- sum(xx[, 10])
    col11 <- sum(xx[, 11])
    col12 <- unique(xx[, 12])
    col13 <- unique(xx[, 13])
    
    newxx <- cbind.data.frame(col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13)
    colnames(newxx) <- colnames(mycoords)
    
    print(which(unique(mycoords_f$V14) == x))
    
    newxx
    
  }
  
  # newcoords <- read.table("ru_new.coords", header = F, stringsAsFactors = F)
  
  newcoords_f <- newcoords[newcoords$V10 >= 60 | newcoords$V11 >= 60, ]
  
  # reshape newcoords_f, small contig first, large contig second
  coords_large <- newcoords_f[newcoords_f$V8 > newcoords_f$V9, ]
  coords_small <- newcoords_f[newcoords_f$V8 <= newcoords_f$V9, ]
  
  coords_large_small <- cbind.data.frame(coords_large[, 4:6],
                                         coords_large[, 1:3],
                                         coords_large[, 7],
                                         coords_large[, 9],
                                         coords_large[, 8],
                                         coords_large[, 11],
                                         coords_large[, 10],
                                         coords_large[, 13],
                                         coords_large[, 12])
  
  colnames(coords_large_small) <- colnames(coords_small)
  
  coords_final <- rbind.data.frame(coords_small, coords_large_small)
  coords_final$V14 <- paste(coords_final$V12, coords_final$V13, sep = "&")
  coords_final <- coords_final[!duplicated(coords_final$V14), ]
  
  write.table(unique(coords_final$V12), "haplotigs.txt", sep = "\t", row.names = F, col.names = F, quote = F)
  write.table(coords_final[, 12:13], "haplotig_pairs.txt", sep = "\t", row.names = F, col.names = F, quote = F)
  
}


#############################################################################################################################
parse_busco <- function(bus=bus, myfai=myfai){
  
  # check how many duplicated busco genes in haplotig
  print(paste("Number of duplicated BUSCOs:", length(unique(bus[bus$V2 == "Duplicated", 1])), sep = " ")) # 1164 duplicated genes
  bus_dup <- bus[bus$V2 == "Duplicated", ]
  
  # check whether duplicates are all having two copies
  test <- c()
  for (i in 1: length(unique(bus_dup$V1))){
    xx <- bus_dup[bus_dup$V1 == unique(bus_dup$V1)[i], ]
    test[i] <- dim(xx)[1]
  }
  
  print(paste("duplication level:", sort(unique(test)), sep = " "))
  
  # gene paris, some genes have more than two copies
  gp <- unique(bus_dup$V1)[which(test == 2)]
  bus_gp <- bus_dup[bus_dup$V1 %in% gp, ]
  
  # make the pair
  mypair <- NULL
  for (i in 1: length(gp)){
    yy <- bus_gp[bus_gp$V1 == gp[i], ]
    zz <- cbind.data.frame(yy[1, c(1, 3:7)], yy[2, c(3:7)])
    mypair <- rbind.data.frame(mypair, zz)
  }
  
  Len_1 <- c()
  for (i in 1:length(mypair[, 2])){
    Len_1[i] <- myfai[myfai$V1 == mypair[, 2][i], 2]
  }
  
  Len_2 <- c()
  for (i in 1:length(mypair[, 7])){
    Len_2[i] <- myfai[myfai$V1 == mypair[, 7][i], 2]
  }
  
  colnames(mypair) <- c("gene", "an", "as", "ae", "ascore", "aglen",
                        "bn", "bs", "be", "bscore", "bglen")
  mypair$aslen <- Len_1
  mypair$bslen <- Len_2
  
  # rebuild the pair, first B and then A, A is maller than B, large contig first
  busco_redundant <- NULL
  for (i in 1:dim(mypair)[1]){
    pp <- mypair[i, ]
    gn <- pp$gene
    if(pp$aslen <= pp$bslen){
      an <- pp$an
      as <- pp$as
      ae <- pp$ae
      aglen <- pp$aglen
      aslen <- pp$aslen
      bn <- pp$bn
      bs <- pp$bs
      be <- pp$be
      bglen <- pp$bglen
      bslen <- pp$bslen
    } else{
      an <- pp$bn
      as <- pp$bs
      ae <- pp$be
      aglen <- pp$bglen
      aslen <- pp$bslen
      bn <- pp$an
      bs <- pp$as
      be <- pp$ae
      bglen <- pp$aglen
      bslen <- pp$aslen
    }
    aabb <- c(gn, bn, bs, be, bglen, bslen, an, as, ae, aglen, aslen)
    busco_redundant <- rbind(busco_redundant, aabb)
  }
  
  # convert matrix to dataframe
  busco_redundant <- as.data.frame(busco_redundant, stringsAsFactors = F)
  busco_redundant <- type.convert(busco_redundant)
  busco_redundant <- busco_redundant[, c(2, 3, 4, 7, 8, 9)]
  busco_redundant <- busco_redundant[order(busco_redundant$V2, busco_redundant$V3), ]
  busco_redundant <- busco_redundant[busco_redundant$V2 != busco_redundant$V7, ]
  # busco_redundant$V10 <- paste(busco_redundant$V2, busco_redundant$V7, sep = "&")
  # busco_redundant <- busco_redundant[busco_redundant$V10 %in% busco_redundant$V10[duplicated(busco_redundant$V10)], ]
  # busco_redundant <- busco_redundant[, 1:6]
  
  write.table(busco_redundant, "busco_pairs.txt", quote = F, sep = "\t", row.names = F, col.names = F)
  
  return(busco_redundant)
  
}

#############################################################################################################################
# a new function to remove allelic contigs based on busco genes
remove_redundant_busco_new <- function(bus = bus, haplotig = haplotig, purge_haplotig_primary = purge_haplotig_primary, secondary = secondary, 
                                       primary = primary, myfai = myfai){
  # print(str(bus))
  # print(str(haplotig))
  # print(str(purge_haplotig_primary))
  # print(str(secondary))
  # print(str(primary))
  # print(str(myfai))
  

  # check how many duplicated busco genes in haplotig
  length(unique(bus[bus$V2 == "Duplicated", 1])) # 1330 duplicated genes
  bus_dup <- bus[bus$V2 == "Duplicated", ]
  bus_single <- bus[bus$V2 == "Complete", ]
  
  # check whether duplicates are all having two copies
  test <- c()
  for (i in 1: length(unique(bus_dup$V1))){
    xx <- bus_dup[bus_dup$V1 == unique(bus_dup$V1)[i], ]
    test[i] <- dim(xx)[1]
  }
  
  # print(unique(test)) # 2 3 4 5 6 8, there are some busco genes having more than 2 copies
  
  # gene paris, some genes have more than two copies
  gp <- unique(bus_dup$V1)[which(test == 2)]
  bus_gp <- bus_dup[bus_dup$V1 %in% gp, ]
  # print(str(bus_gp))
  
  # suspect haplotigs from purge haplotig and cbp
  hp <- as.data.frame(unique(c(haplotig$V1, secondary$V1)), stringsAsFactors = FALSE)
  colnames(hp) <- "haplotigs"
  # print(str(hp))
  
  # primary contigs from both purge haplotigs and collinear block
  total_primary <- as.data.frame(unique(c(primary$V1, purge_haplotig_primary$V1)), stringsAsFactors = FALSE)
  colnames(total_primary) <- 'primary'
  # print(str(total_primary))
  
  # check how many duplicated busco genes are left
  bus_left <- bus_gp[!bus_gp$V3 %in% hp$haplotigs, ]
  
  # based on busco genes, some left undetected allelic contigs by purge_haplotigs and collinear block can be found
  # for a pair of allelic contigs detected by busco genes, the one with shorter length is assigned to haplotig
  hap_bus <- c()
  for(x in unique(bus_left$V1)){
    xx <- bus_left[bus_left$V1 == x, ]
    if (length(unique(xx$V3)) > 1){
      yy <- myfai[myfai$V1 %in% unique(xx$V3), ]
      yy <- yy[order(yy$V2, decreasing = T), ]
      hap_bus <- c(hap_bus, yy[-1, 1])
    }
  }
  
  # check the total length of busco detected haplotigs
  hap_bus <- myfai[myfai$V1 %in% hap_bus, ]
  
  hap_bus <- hap_bus[!hap_bus$V1 %in% c(total_primary$primary, hp$haplotigs), ]
  # print(sum(hap_bus$V2)) # 281541568
  
  # total length of all haplotigs
  # print(sum(myfai[myfai$V1 %in% unique(c(hap_bus$V1, hp$haplotigs)), 2])) # 2598808938
  
  # well, some haplotigs contains single copy complete busco genes, maybe rescue them if they don't add too much duplication to the final primay assembly
  bus_rescue <- intersect(bus_single$V3, hap_bus$V1)
  # if number of complete busco gene is greater of equal than duplicated busco gene, then save this contig
  # bus_rescue_final <- c()
  # for(x in bus_rescue){
  #   
  #   if(length(bus_single[bus_single$V3 == x, 1]) >= length(bus_dup[bus_dup$V3 == x, 1])){
  #     bus_rescue_final <- c(bus_rescue_final, x)
  #   }
  #   
  # }
  
  # make a final haplotig list
  final_haplotigs <- setdiff(c(hap_bus$V1, hp$haplotigs), bus_rescue)
  print(paste("Size of final allelic contigs", sum(myfai[myfai$V1 %in% final_haplotigs, 2]), sep = ": "))
  
  
  # final primary contigs
  final_primary_contigs <- setdiff(myfai$V1, final_haplotigs)
  print(paste("Size of final primary contigs", sum(myfai[myfai$V1 %in% final_primary_contigs, 2]), sep = ": "))
  
  # check how many duplicated busco genes are left
  bus_left <- bus[!bus$V3 %in% final_haplotigs, ]
  
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
  
  # write out final haplotigs and primary contigs
  write.table(final_haplotigs, "final_suspect_haplotigs.txt", quote = F, row.names = F, col.names = F, sep = "\t")
  write.table(final_primary_contigs, "final_primary_contigs.txt", quote = F, row.names = F, col.names = F, sep = "\t")
  
  # add a column to show which method the haplotig was found
  final_haplotigs <- as.data.frame(final_haplotigs, stringsAsFactors = FALSE)
  colnames(final_haplotigs) <- "haplotigs"
  # print(final_haplotigs[1:5, ])
  
  xx <- c()
  for(x in final_haplotigs$haplotigs){
    
    if(x %in% setdiff(haplotig$V1, secondary$V1)){
      xx <- c(xx, "purge_haplotig")
    }else if(x %in% setdiff(secondary$V1, haplotig$V1)){
      xx <- c(xx, "mcscanx")
    }else if(x %in% intersect(haplotig$V1, secondary$V1)){
      xx <- c(xx, "phmc") # phmc means puge_haplotig + mcscanx
    }else if(x %in% hap_bus$V1){
      xx <- c(xx, "busco")
    }
    
  }
  
  final_haplotigs$methods <- xx
  
  # now add another column to show how many duplicated busco genes in each haplotig
  xx <- c()
  yy <- c()
  for(x in final_haplotigs$haplotigs){
    xx <- c(xx, length(bus[bus$V3 == x & bus$V2 == "Duplicated", 3]))
    yy <- c(yy, length(bus[bus$V3 == x & bus$V2 == "Complete", 3]))
  }
  
  final_haplotigs$complete <- yy
  final_haplotigs$duplicated <- xx
  
  # write out the final haplotigs table to show which method found which haplotig
  write.table(final_haplotigs, "final_suspect_haplotigs_table.txt", quote = F, row.names = F, col.names = T, sep = "\t")
}