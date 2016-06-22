#' This function localize seqments of continous integer values
#' 
#' @param data A list of integer vectors of zeros and ones
#' @param minSeg A minimal length of continuous segments of a certain value
#' @param smooth A number of subsequent integers considered as errors
#' @import GenomicRanges
#' @importFrom fastseg fastseg
#' @author David Porubsky
#' @export
 
findSegmentsCBC <- function(data=NULL, minSeg=0, smooth=3) {
  
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
  for (i in 1:length(data)) {
    chrom.binary <- data[[i]]
    index <- names(data[i])
    comp.vector <- as.numeric(chrom.binary[,5])
    
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
      gen.ranges <- IRanges(start=chrom.binary$start[segs$startRow], end=chrom.binary$end[segs$endRow])
      ranges(segs) <- gen.ranges
      segs$index <- index
    
      segs$direction[segs$seg.mean >= 0.75] <- 'w'
      segs$direction[segs$seg.mean <= 0.25] <- 'c'
      segs$direction[segs$seg.mean > 0.25 & segs$seg.mean < 0.75] <- 'wc'  
    
      segments[[index]] <- segs
    } else {
      #TODO
    }
  }
    
  recombs <- GenomicRanges::GRangesList()
  for (i in 1:length(segments)) {
    segm <- segments[[i]]
    index <- names(segments[i])
    if (length(segm) > 1) {
      segm <- collapseBins(gr = segm, id.field = 7)
      segments[[i]] <- segm
    } 
      
    if (length(segm) > 1) {
      suppressWarnings( recomb <- gaps(segm) )
      start(recomb) <- start(recomb)-1
      end(recomb) <- end(recomb)+1
      recomb$index <- index
      recombs[[index]] <- recomb[-1]
    }  
  }
  
  #seqments.gr <- GRanges()
  seqments.gr <- unlist(segments)
  names(seqments.gr) <- NULL
  #segments.df <- data.frame()
  segments.df <- as(seqments.gr, "data.frame")
  segments.df <- data.frame(chromosome=segments.df$index, start=segments.df$start, end=segments.df$end, num.read=segments.df$num.mark, seg.mean=segments.df$seg.mean, strand=segments.df$direction)
  
  #recombs.gr <- GRanges()
  recombs.gr <- unlist(recombs)
  names(recombs.gr) <- NULL
  #recombs.df <- data.frame()
  recombs.df <- as(recombs.gr, "data.frame")
  recombs.df <- data.frame(chromosome=recombs.df$index, start=recombs.df$start, end=recombs.df$end, range=recombs.df$width)
  
  results <- list()
  results[['segments']] <- segments.df
  results[['recombs']] <- recombs.df
  return(results)
# end
}

