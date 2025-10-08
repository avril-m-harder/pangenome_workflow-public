args <- commandArgs(trailingOnly = TRUE)

dat <- read.table(args[1], comment.char = '?')
fais <- list.files(path = args[2], pattern = '.fai')

dat <- do.call(rbind, strsplit(dat$V1, split = '#'))
dat <- cbind(dat[,1], do.call(rbind, strsplit(dat[,3], split = ':')))
dat <- cbind(dat[,c(1:2)], do.call(rbind, strsplit(dat[,3], split = '-')))

## make this an arg
h <- max(6, length(unique(dat[,2]))/2.4)

pdf(args[3], width = 7, height = h)
for(s in unique(dat[,1])){
  ## this will break if sample name in fasta name != path name in 'dat';
  ## be sure to do this perfectly going forward -- base names only
  idx <- read.table(paste0(args[2],'/',fais[grep(s, fais)]))
  idx <- idx[grep('Chr', idx[,1]),]
  sub <- dat[dat[,1] == s,]
  
  plot(0, 0, col = 'transparent', xlim = c(0, max(idx[,2])), ylim = c(0.5, nrow(idx)),
       xlab = 'Chromosome position (bp)', xaxt = 'n', ylab = '', bty = 'n', yaxt = 'n',
       main = s)
    if(max(idx[,2]) > 10e6){
      axis(1, at = seq(from = 0, to = max(idx[,2]), by = 1e7))
    } else{
      axis(1, at = seq(from = 0, to = max(idx[,2]), by = 1e6))
    }
    axis(2, at = c(1:nrow(idx)), labels = rev(idx[,1]), las = 1, lwd = 0)
    y <- 1
    for(r in rev(seq(nrow(idx):1))){
      chr <- idx[r,1]
      real.len <- idx[idx[,1] == chr, 2]
      polygon(c(0, real.len, real.len, 0), c(y-0.4, y-0.4, y+0.4, y+0.4))
      temp <- sub[sub[,2] == chr,]
      if(is.null(nrow(temp))){
        polygon(c(as.numeric(temp[3]), as.numeric(temp[4]), as.numeric(temp[4]), as.numeric(temp[3])),
                c(y-0.4, y-0.4, y+0.4, y+0.4),
                border = NA, col = 'dodgerblue3')
      } else{
        for(p in 1:nrow(temp)){
          polygon(c(as.numeric(temp[p,3]), as.numeric(temp[p,4]), as.numeric(temp[p,4]), as.numeric(temp[p,3])),
                  c(y-0.4, y-0.4, y+0.4, y+0.4),
                  border = NA, col = 'dodgerblue3')
        }
      }
      y <- y+1  
    }
}
dev.off()
