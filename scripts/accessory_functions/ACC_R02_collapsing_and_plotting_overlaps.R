if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager', lib = '/home/aharder/.Rlibs')
if (!require('GenomicRanges'))
	BiocManager::install('GenomicRanges', lib = '/home/aharder/.Rlibs')
library('GenomicRanges')

args <- commandArgs(trailingOnly = TRUE)
taxon <- args[1]

`%notin%` <- Negate(`%in%`)

##### Summarize results of bedtools intersect #####
samps <- read.table('clipped_analysis_files.list')
samps <- samps$V1

R.OUT <- NULL
G.OUT <- NULL

for(s in samps){
  samp <- s
  
  idx <- read.table(paste0(samp,'.genome'))
  idx <- idx[grep('Chr', idx[,1]),]
  tot.genom.len <- sum(as.numeric(idx$V2))
  
  for(t in c('repeats','genes')){
    dat.type <- t
    intersects <- read.table(paste0(samp,'_clips_',t,'_intersects.txt'), sep = '\t')
    
    ##### process repeats #####
    if(dat.type == 'repeats'){
      
      colnames(intersects) <- c('clip.chrom','clip.start','clip.end',
      'rep.chrom','rep.method','basis','rep.start','rep.end',
      'rep.score','rep.strand','rep.frame','rep.info','olap.len')

      ## keep features that overlap with clipped sequences (i.e., overlap length > 0)
      intersects <- intersects[intersects$olap.len > 0,]
      
      ##### >>> collapse overlapping annotations and calc summary stats #####
      ## get repeat and clipped ranges and intersect them (bedtools isec doesn't provide collapsed overlap coords)
      feature.ranges <- GRanges(seqnames = intersects$rep.chrom,
                                ranges = IRanges(start = intersects$rep.start, end = intersects$rep.end))
      feature.ranges <- reduce(feature.ranges)
      
      clipped.ranges <- GRanges(seqnames = intersects$clip.chrom,
                                ranges = IRanges(start = intersects$clip.start, end = intersects$clip.end))
      clipped.ranges <- reduce(clipped.ranges)
      
      internal.ints <- intersect(feature.ranges, clipped.ranges)
      ## total length of intersections between features (repeats or genes) and clipped sequences
      non.overlap.int.len <- sum(internal.ints@ranges@width)
      
      clips <- read.table(paste0(samp,'_subpath_excluded_coords.bed'), sep = '\t')
      ## total length of clipped sequences
      tot.clip.len <- sum(as.numeric(clips$V3) - as.numeric(clips$V2))
      
      ## proportion of the input genome that is clipped out
      prop.genom.clip <- tot.clip.len/tot.genom.len
      
      ## proportion of clipped sequence that overlap annotated repeats
      prop.clip.repeats <- non.overlap.int.len/tot.clip.len
      
      ## check proportions of clipped repeats belonging to different categories (not accounting for overlap here,
      ## so really describing how all annotations overlapping clipped repeats are categorized)
      if(nrow(intersects) > 100e3){
        SAVE <- NULL
        for(r in 1:nrow(intersects)){
          SAVE <- c(SAVE, gsub('class=', '', unlist(strsplit(intersects$rep.info[r], split = ';'))[4]))
        }
        ints <- cbind(SAVE, intersects$olap.len)
        tot.ints <- sum(as.numeric(ints[,2]))
        rm(SAVE)
      } else{
        ints <- cbind(gsub('class=','',do.call(rbind, strsplit(intersects$rep.info, split = ';'))[,4]), intersects$olap.len)
        tot.ints <- sum(as.numeric(ints[,2]))
      }
      
      ltrs <-   ints[(grep('LTR', ints[,1])),]
      dnas <- ints[grep('DNA', ints[,1]),]
      prop.dnas <- sum(as.numeric(dnas[,2]))/tot.ints
      prop.ty3 <- sum(as.numeric(ltrs[ltrs[,1] == 'LTR/Gypsy', 2]))/tot.ints
      prop.copia <- sum(as.numeric(ltrs[ltrs[,1] == 'LTR/Copia', 2]))/tot.ints
      ltrs <- ltrs[which(ltrs[,1] != 'LTR/Gypsy' & ltrs[,1] != 'LTR/Copia'),]
      prop.other.ltrs <- sum(as.numeric(ltrs[,2]))/tot.ints
      others <- ints[which(ints[,1] %notin% c('LTR/Gypsy', 'LTR/Copia', dnas[,1], ltrs[,1])),]
      prop.others <- sum(as.numeric(others[,2]))/tot.ints
      if(sum(prop.dnas, prop.ty3, prop.copia, prop.other.ltrs, prop.others) != 1){
        print(paste0(samp,' - PROPORTION CALCULATION ERROR'))
      }
      
      save <- c(samp, tot.genom.len, tot.clip.len, prop.genom.clip, prop.clip.repeats, 
                prop.ty3, prop.copia, prop.other.ltrs, prop.dnas, prop.others)
      R.OUT <- rbind(R.OUT, save)
    
    ##### process genes #####
    } else{
      ## need to narrow down to 1 entry per locus, two ways:
      ##  i) genes
      ## ii) CDS
      
      colnames(intersects) <- c('clip.chrom','clip.start','clip.end',
                                'gene.chrom','gene.method','feature.type','gene.start','gene.end',
                                'gene.score','gene.strand','gene.frame','gene.info','olap.len')
      
      ## keep features that overlap with clipped sequences (i.e., overlap length > 0)
      intersects <- intersects[intersects$olap.len > 0,]
      
      genes <- intersects[intersects$feature.type == 'gene',]
      gene.ranges <- GRanges(seqnames = genes$gene.chrom,
                             ranges = IRanges(start = genes$gene.start, end = genes$gene.end))
      gene.ranges <- reduce(gene.ranges)
      clipped.ranges <- GRanges(seqnames = intersects$clip.chrom,
                                ranges = IRanges(start = intersects$clip.start, end = intersects$clip.end))
      clipped.ranges <- reduce(clipped.ranges)
      internal.ints <- intersect(gene.ranges, clipped.ranges)
      ## total length of intersections between features (repeats or genes) and clipped sequences
      non.overlap.int.len <- sum(internal.ints@ranges@width)
      prop.clip.genes <- non.overlap.int.len/tot.clip.len
      
      cds <- intersects[intersects$feature.type == 'CDS',]
      cds.ranges <- GRanges(seqnames = cds$gene.chrom,
                            ranges = IRanges(start = cds$gene.start, end = cds$gene.end))
      cds.ranges <- reduce(cds.ranges)
      clipped.ranges <- GRanges(seqnames = intersects$clip.chrom,
                                ranges = IRanges(start = intersects$clip.start, end = intersects$clip.end))
      clipped.ranges <- reduce(clipped.ranges)
      internal.ints <- intersect(cds.ranges, clipped.ranges)
      ## total length of intersections between features (repeats or genes) and clipped sequences
      non.overlap.int.len <- sum(internal.ints@ranges@width)
      prop.clip.cds <- non.overlap.int.len/tot.clip.len
      save <- c(samp, tot.genom.len, tot.clip.len, prop.genom.clip, prop.clip.genes, prop.clip.cds)
      G.OUT <- rbind(G.OUT, save)
    }
  }
}

R.sum.stats <- as.data.frame(R.OUT)
colnames(R.sum.stats) <- c('samp', 'tot.genom.len', 'tot.clip.len', 'prop.genom.clip', 'prop.clip.repeats', 
                           'prop.ty3', 'prop.copia', 'prop.other.ltrs', 'prop.dnas', 'prop.others')
for(c in 2:ncol(R.sum.stats)){
  R.sum.stats[,c] <- as.numeric(R.sum.stats[,c])
}

G.sum.stats <- as.data.frame(G.OUT)
colnames(G.sum.stats) <- c('samp', 'tot.genom.len', 'tot.clip.len', 'prop.genom.clip', 'prop.clip.genes','prop.clip.cds')
for(c in 2:ncol(G.sum.stats)){
  G.sum.stats[,c] <- as.numeric(G.sum.stats[,c])
}
R.sum.stats <- R.sum.stats[order(R.sum.stats$tot.clip.len, decreasing = FALSE),]
G.sum.stats <- G.sum.stats[order(G.sum.stats$tot.clip.len, decreasing = FALSE),]
if(nrow(R.sum.stats) < 5){
	wid <- nrow(R.sum.stats)*2
} else{
	wid <- nrow(R.sum.stats)*1
}
wid <- min(wid, 12)
pdf(paste0(taxon,'_clipped_sequence_summary.pdf'), width = 6, height = wid)
bmar <- max(nchar(R.sum.stats$samp))/1.5
par(mar = c(5, bmar, 6, 1)+0.1)
plot(0, 0, col = 'transparent', ylim = c(0.75, nrow(R.sum.stats)+0.25), xlim = c(0, max(R.sum.stats$prop.genom.clip)),
     xlab = 'Proportion of input genome',
     yaxt = 'n', ylab = '', main = '')
  for(r in 1:nrow(R.sum.stats)){
    prop.clip <- R.sum.stats$prop.genom.clip[r]
    prop.repeat <- R.sum.stats$prop.genom.clip[r]*R.sum.stats$prop.clip.repeats[r]
    prop.gene <- G.sum.stats$prop.genom.clip[r]*G.sum.stats$prop.clip.genes[r]
    prop.cds <- G.sum.stats$prop.genom.clip[r]*G.sum.stats$prop.clip.cds[r]
    points(prop.clip, r-.15, pch = 16, cex = 1.5, col = 'grey30')
    lines(c(0, prop.clip), c(r-.15, r-.15), lwd = 4, col = 'grey30')
    points(prop.repeat, r-.05, pch = 16, cex = 1.5, col = 'chocolate1')
    lines(c(0, prop.repeat), c(r-.05, r-.05), lwd = 4, col = 'chocolate1')
    points(prop.gene, r+.05, pch = 16, cex = 1.5, col = 'dodgerblue2')
    lines(c(0, prop.gene), c(r+.05, r+.05), lwd = 4, col = 'dodgerblue2')
    points(prop.gene, r+.15, pch = 16, cex = 1.5, col = 'springgreen3')
    lines(c(0, prop.gene), c(r+.15, r+.15), lwd = 4, col = 'springgreen3')
    axis(2, at = r, labels = R.sum.stats$samp[r], las = 2)
  }

  par(xpd = TRUE)
  legend('top', legend = c('Clipped sequence','Repetitive clipped sequence','Genic clipped sequence (\'gene\')',
                                'Genic clipped sequence (\'CDS\')'),
         col = c('grey30','chocolate1','dodgerblue2','springgreen3'), pch = 16, pt.cex = 1.5, inset = -0.11)
dev.off()

temp <- merge(R.sum.stats, G.sum.stats[,c(1,5,6)], by = 'samp')
write.table(temp, paste0(taxon,'_clipped_sequence_summary.txt'), sep = '\t', quote = FALSE, row.names = FALSE)