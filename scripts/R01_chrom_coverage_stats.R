args <- commandArgs(trailingOnly = TRUE)

dat <- read.table(args[1])
out.pdf <- args[2]

ref <- unlist(strsplit(gsub('id=', '', dat[grep('ref-contig', dat[,2]),][1,5]), '|', fixed = TRUE))[1]
dat <- dat[,c(5,6,7)]
dat[,1] <- gsub('id=', '', dat[,1])
dat[,2] <- gsub('len=', '', dat[,2])
dat[,3] <- gsub('cov=', '', dat[,3])
dat <- as.data.frame(cbind(do.call(rbind, strsplit(dat[,1], split = '|', fixed = TRUE)), dat[,2], dat[,3]))
colnames(dat) <- c('sample','chrom','length','covg')
dat$length <- as.numeric(dat$length)
dat$covg <- as.numeric(dat$covg)
dat <- dat[order(dat$chrom, dat$sample),]

n.samps <- length(unique(dat$sample))
y.dim <- max(5, n.samps/2.5)
x.dim <- max(6, max(nchar(dat$sample))/3)
l.mar <- max(9, max(nchar(dat$sample))/1.8)

pdf(out.pdf, width = x.dim, height = y.dim)
par(mar=c(4, l.mar, 4, 1)+0.1)
for(c in unique(dat$chrom)){
  sub <- dat[dat$chrom == c,]
  plot(0, 0, ylim = c(1, length(unique(dat$sample))), xlim = c(min(dat$covg), max(dat$covg)), 
       col = 'transparent', ylab = '', yaxt = 'n', xlab = 'Chromosome coverage in graph',
       main = c)
    mtext(unique(sub$sample), side = 2, at = c(1:length(unique(sub$sample))), las = 2, line = 0.25)
    for(r in 1:nrow(sub)){
      if(sub$sample[r] == ref){
        lines(c(0, sub$covg[r]), c(r, r), lwd = 3, col = 'darkorchid4')
        points(sub$covg[r], r, pch = 16, cex = 1.5, col = 'darkorchid4')
      } else{
        lines(c(0, sub$covg[r]), c(r, r), lwd = 3, col = 'darkorange3')
        points(sub$covg[r], r, pch = 16, cex = 1.5, col = 'darkorange3')
      }
    }
}
dev.off()