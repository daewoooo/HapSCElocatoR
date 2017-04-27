#' This function localize seqments of continous integer values
#' 
#' @param data A list of integer vectors of zeros and ones
#' @param minSeg A minimal length of continuous segments of a certain value
#' @param smooth A number of subsequent integers considered as errors
#' @param read.len
#' @import GenomicRanges
#' @importFrom fastseg fastseg
#' @author David Porubsky
#' @export
 
findSegmentsCBC <- function(data=NULL, minSeg=0, minSize=0, smooth=3, read.len=100) {
  
  switchValue <- function(x) {
    if (x == 1) {
      x <- 0
    } else {
      x <- 1  
    } 
  }
  
  collapseBins <- function(gr, id.field=0) {
    ind.last <- cumsum(runLength(Rle(mcols(gr)[,id.field]))) ##get indices of last range in a consecutive(RLE) run of the same value
    ind.first <- c(1,cumsum(runLength(Rle(mcols(gr)[,id.field]))) + 1) ##get indices of first range in a consecutive(RLE) run of the same value
    ind.first <- ind.first[-length(ind.first)]  ##erase last index from first range indices 
    collapsed.gr <- GenomicRanges::GRanges(seqnames=seqnames(gr[ind.first]), ranges=IRanges(start=start(gr[ind.first]), end=end(gr[ind.last])), mcols=mcols(gr[ind.first]))
    names(mcols(collapsed.gr)) <- names(mcols(gr[ind.first]))
    return(collapsed.gr)
  }
  
  #localize segments in 0/1 vector
  segments <- GenomicRanges::GRangesList()
  recombs <- GenomicRanges::GRangesList()
  for (i in 1:length(data)) {
    chrom.gr <- data[[i]]
    index <- names(data[i])
    
    if (length(chrom.gr) <= 2*minSeg) next
    
    #Translate read directionality into a binary vector (0|1)
    frag.df <- as(chrom.gr, "data.frame")
    frag.df$strand <- as.character(frag.df$strand)
    frag.df$strand[frag.df$strand == "+"] <- 0
    frag.df$strand[frag.df$strand == "-"] <- 1
    comp.vector <- as.numeric(frag.df[,5])
    
    if (minSize != 0 & minSeg==0) {
      tiles <- unlist(tileGenome(seqlengths(chrom.gr)[i], tilewidth = minSize))
      #suppressWarnings( tiles.shift <- shift(tiles, shift =  minSeg/2) )
      #suppressWarnings( tiles <- c(tiles, tiles.shift) )
      counts <- countOverlaps(tiles, chrom.gr)
      minSeg <- max(10, round(mean(counts[counts>0], trim=0.05)))
    }
    
    if (length(comp.vector) > 2*minSeg) {
      segs <- fastseg(comp.vector, minSeg=minSeg)
      
      while (any(segs$num.mark <= smooth)) {
        toSwitch <- which(segs$num.mark <= smooth)
        switch.segs <- segs[toSwitch]
        switch.pos <- mapply(function(x,y) {x:y}, x=switch.segs$startRow, y=switch.segs$endRow)
        switch.pos <- unlist(switch.pos)
        
        switched.vals <- sapply(comp.vector[switch.pos], switchValue) #SWITCH
        #comp.vector <- comp.vector[-switch.pos]  #DELETE
        comp.vector[switch.pos] <- switched.vals
        segs <- fastseg(comp.vector, minSeg=minSeg)
      }
      
      chrom.gr$breakpointBefore <- 0
      chrom.gr$breakpointAfter <- 0
      chrom.gr[strand(chrom.gr) == "+"]$breakpointBefore <- end(chrom.gr[strand(chrom.gr) == "+"])
      chrom.gr[strand(chrom.gr) == "+"]$breakpointAfter <- start(chrom.gr[strand(chrom.gr) == "+"])
      chrom.gr[strand(chrom.gr) == "-"]$breakpointBefore <- start(chrom.gr[strand(chrom.gr) == "-"])
      chrom.gr[strand(chrom.gr) == "-"]$breakpointAfter <- end(chrom.gr[strand(chrom.gr) == "-"])
      
      #gen.ranges <- IRanges(start=end(chrom.gr)[segs$startRow], end=start(chrom.gr)[segs$endRow])
      gen.ranges <- IRanges(start=chrom.gr[segs$startRow]$breakpointBefore, end=chrom.gr[segs$endRow]$breakpointAfter)
      #gen.ranges <- IRanges(start=start(chrom.gr)[segs$startRow], end=end(chrom.gr)[segs$endRow])
      #gen.ranges <- disjoin(gen.ranges)
      #mask <- which(width(gen.ranges)>read.len)
      #gen.ranges <- gen.ranges[mask]
      ranges(segs) <- gen.ranges
      segs$index <- index
      
      segs$direction[segs$seg.mean >= 0.75] <- 'w'
      segs$direction[segs$seg.mean <= 0.25] <- 'c'
      segs$direction[segs$seg.mean > 0.25 & segs$seg.mean < 0.75] <- 'wc'  
      
      segments[[index]] <- segs
    } else {
      #TODO
    }
    
    #recombs <- GenomicRanges::GRangesList()
    seqlevels(recombs) <- seqlevels(segments)
    for (i in 1:length(segments)) {
      segm <- segments[[i]]
      index <- names(segments[i])
      if (length(segm) > 1) {
        segm <- collapseBins(gr = segm, id.field = 7)
        segments[[i]] <- segm
      } 
      
      if (length(segm) > 1) {
        #suppressWarnings( recomb <- gaps(segm) )
        recomb <- GRanges(seqnames=seqnames(segm), ranges=IRanges(start=c(1, end(segm)[-length(segm)]), end=start(segm)))
        start(recomb) <- start(recomb)-1
        end(recomb) <- end(recomb)+1
        recomb$index <- index
        recomb <- recomb[-1]
        recomb$genoT <- paste(segm$direction[1:length(segm$direction)-1], segm$direction[2:length(segm$direction)], sep='-')
        recombs[[index]] <- recomb
      }  
    }
  }
  
  #seqments.gr <- GRanges()
  segments.gr <- unlist(segments)
  segments.gr <- GRanges(seqnames=segments.gr$index, ranges=ranges(segments.gr), num.read=segments.gr$num.mark, seg.mean=segments.gr$seg.mean, direction=segments.gr$direction)
  names(segments.gr) <- NULL
  #segments.df <- data.frame()
  #segments.df <- as(segments.gr, "data.frame")
  #segments.df <- data.frame(chromosome=segments.df$index, start=segments.df$start, end=segments.df$end, num.read=segments.df$num.mark, seg.mean=segments.df$seg.mean, strand=segments.df$direction)
  
  #recombs.gr <- GRanges()
  recombs.gr <- unlist(recombs)
  recombs.gr <- GRanges(seqnames=recombs.gr$index, ranges=ranges(recombs.gr), genoT=recombs.gr$genoT) 
  names(recombs.gr) <- NULL
  #recombs.df <- data.frame()
  #recombs.df <- as(recombs.gr, "data.frame")
  #recombs.df <- data.frame(chromosome=recombs.df$index, start=recombs.df$start, end=recombs.df$end, range=recombs.df$width)
  
  results <- list()
  seqlengths(segments.gr) <- seqlengths(data)[seqlevels(segments.gr)]
  results[['segments']] <- segments.gr
  seqlengths(recombs.gr) <- seqlengths(data)[seqlevels(recombs.gr)]
  results[['recombs']] <- recombs.gr
  return(results)
# end
}

