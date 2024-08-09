if (!require('scales', quietly = TRUE))
    install.packages('scales', lib = '/home/aharder/.Rlibs')
library('scales')

args <- commandArgs(trailingOnly = TRUE)
taxon <- args[1]

`%notin%` <- Negate(`%in%`)

by.samp.fn <- paste0(taxon,'_prop_retained_by_sample.pdf')
by.chrom.fn <- paste0(taxon,'_prop_retained_by_chrom.pdf')

##### Get sample names #####
incl.beds <- list.files(path = '.', pattern = '_subpath_included_coords.sorted.bed')
samps <- gsub('_subpath_included_coords.sorted.bed', '', incl.beds)

OUT <- NULL
OUT1 <- NULL
##### 1. Get output of this loop and process on pants with bedtools intersect for now #####
for(samp in samps){
  ## read in coords for retained sequences
  kept <- read.table(paste0(samp,'_subpath_included_coords.sorted.bed'))
  genom <- read.table(paste0(samp,'.genome'))

  tot.kept <- sum(as.numeric(kept[,3]) - as.numeric(kept[,2]) + 1)/sum(as.numeric(genom$V2))
  save <- c(samp, tot.kept)
  OUT <- rbind(OUT, save)
  
  for(c in unique(kept[,1])){
    tmp <- kept[kept[,1] == c,]
    c.kept <- sum(as.numeric(tmp[,3]) - as.numeric(tmp[,2]) + 1)/sum(as.numeric(genom[genom$V1 == c, 2]))
    save <- c(samp, c, c.kept)
    OUT1 <- rbind(OUT1, save)
  }
}

## plot proportions retained by sample
OUT <- OUT[order(OUT[,2], decreasing = TRUE),]
colour <- 'springgreen4'
offset = 0.2
if(length(samps) < 5){
  plot.ht <- length(samps)*2
} else{
  plot.ht <- length(samps)*0.9
}
plot.ht <- min(plot.ht, 12)
pdf(by.samp.fn, width = 5.5, height = plot.ht)
bmar <- max(nchar(samps))/1.5
par(mar = c(5, bmar, 1, 1)+0.1)
plot(0,0, col = 'transparent', 
     xlim = c(min(as.numeric(OUT1[,3])), max(as.numeric(OUT1[,3]))),
     ylim = c(0.8, length(samps)+0.2), yaxt = 'n', xlab = 'Chromosome proportion retained',
     ylab = '')
  axis(2, at = c(length(samps):1), labels = OUT[,1], las = 2)
  y <- length(samps)
  c <- 1
  for(r in 1:nrow(OUT)){
    xvals <- as.numeric(OUT1[OUT1[,1] == OUT[r, 1], 3])
    points(xvals, rep(y, length(xvals)), col = alpha(colour, 0.6), pch = 16, cex = 1.25)
    lines(c(mean(xvals), mean(xvals)), c(y-offset, y+offset), lwd = 4, col = colour)
    y <- y-1
    c <- c+1
  }
dev.off()

## plot proportions retained by chromosome
colour <- 'springgreen4'
offset = 0.2
if(length(genom$V1) < 5){
  plot.ht <- length(genom$V1)*2
} else{
  plot.ht <- length(genom$V1)*0.7
}
plot.ht <- min(plot.ht, 12)
bmar <- max(nchar(genom$V1))/1.5
pdf(by.chrom.fn, width = 5, height = plot.ht)
par(mar = c(5, bmar, 1, 1)+0.1)
plot(0,0, col = 'transparent', 
     xlim = c(min(as.numeric(OUT1[,3])), max(as.numeric(OUT1[,3]))),
     ylim = c(0.8, length(unique(OUT1[,2]))+0.2), yaxt = 'n', xlab = 'Chromosome proportion retained',
     ylab = '')
  axis(2, at = c(length(unique(OUT1[,2])):1), labels = c(unique(OUT1[,2])), las = 2)
  y <- length(unique(OUT1[,2]))
  c <- 1
  for(c in unique(OUT1[,2])){
    print(paste0(c,' - ',y))
    xvals <- as.numeric(OUT1[OUT1[,2] == c, 3])
    points(xvals, rep(y, length(xvals)), col = alpha(colour, 0.6), pch = 16, cex = 1.25)
    lines(c(mean(xvals), mean(xvals)), c(y-offset, y+offset), lwd = 4, col = colour)
    y <- y-1
  }
dev.off()
