#' This function localize seqments of continous integer values
#' 
#' @param data A list of integer vectors of zeros and ones
#' @import GenomicRanges
#' @import HMM
#' @author David Porubsky
#' @export

findSegmentsHMM <- function(data=NULL) {
  
  switchValue <- function(x) {
    if (x == '1') {
      x <- '0'
    } else {
      x <- '1'  
    }
    return(x)
  }
  
  collapseBins <- function(gr, id.field=0) {
    ind.last <- cumsum(runLength(Rle(mcols(gr)[,id.field]))) ##get indices of last range in a consecutive(RLE) run of the same value
    ind.first <- c(1,cumsum(runLength(Rle(mcols(gr)[,id.field]))) + 1) ##get indices of first range in a consecutive(RLE) run of the same value
    ind.first <- ind.first[-length(ind.first)]  ##erase last index from first range indices 
    collapsed.gr <- GenomicRanges::GRanges(seqnames=seqnames(gr[ind.first]), ranges=IRanges(start=start(gr[ind.first]), end=end(gr[ind.last])), mcols=mcols(gr[ind.first]))
    names(mcols(collapsed.gr)) <- names(mcols(gr[ind.first]))
    return(collapsed.gr)
  }
  
  States <- c("0","1","2")
  Symbols <- c("0","1")
  startProbs <- c(0.3,0.3,0.3)
  emissionProbs <- matrix(c(0.95,0.05,0.5,0.05,0.95,0.5),3,2)
  transProbs <- matrix(0.05,3,3)
  diag(transProbs) <- 0.9
  hmm <- initHMM(States, Symbols, startProbs, transProbs, emissionProbs)
  
  segments <- GenomicRanges::GRangesList()
  for (i in 1:length(data)) {
    chrom.binary <- data[[i]]
    index <- names(data[i])
    comp.vector <- as.numeric(chrom.binary[,5])  
    
    if (length(comp.vector) > 1) {
      observed <- as.character(comp.vector)
      postProbs <- posterior(hmm, observed)
    
      maxProbs <- apply(postProbs, 2, max)
      toswitch <- which(maxProbs < 0.5)
    
      maxiter <- 1
      while (length(toswitch)>0 & maxiter<=5) {
        switched.vals <- sapply(observed[toswitch], switchValue)
        observed[toswitch] <- switched.vals
        postProbs <- posterior(hmm, observed)
        maxProbs <- apply(postProbs, 2, max)
        toswitch <- which(maxProbs < 0.5)
        maxiter <- maxiter+1
      } 
    
      obsStates <- apply(postProbs, 2, which.max)
    
      dirs <- obsStates
      dirs[dirs == 1] <- "c"
      dirs[dirs == 2] <- "w"
      dirs[dirs == 3] <- "wc" 
    
      strand <- chrom.binary$strand
      strand[strand == 0] <- "+"
      strand[strand == 1] <- "-"
    
      frag.gr <- GenomicRanges::GRanges(seqnames=chrom.binary$seqnames, IRanges(start=chrom.binary$start, end=chrom.binary$end), strand=strand, state=obsStates, direction=dirs)
      segs <- collapseBins(frag.gr, id.field = 2)
      segments[[index]] <- segs
    }  
  }
  
  recombs <- GenomicRanges::GRangesList()
  seqlevels(recombs) <- seqlevels(segments)
  for (i in 1:length(segments)) {
    segm <- segments[[i]]
    index <- names(segments[i])
    
    if (length(segm) > 1) {
      suppressWarnings( recomb <- gaps(segm) )
      start(recomb) <- start(recomb)-1
      end(recomb) <- end(recomb)+1
      recomb$index <- index
      recombs[[index]] <- recomb[-1]
    }  
  }
      
  seqments.gr <- unlist(segments)
  names(seqments.gr) <- NULL
  segments.df <- as(seqments.gr, "data.frame")
  segments.df <- data.frame(chromosome=segments.df$seqnames, start=segments.df$start, end=segments.df$end, strand=segments.df$direction)
  
  recombs.gr <- unlist(recombs)
  names(recombs.gr) <- NULL
  recombs.df <- as(recombs.gr, "data.frame")
  recombs.df <- data.frame(chromosome=recombs.df$seqnames, start=recombs.df$start, end=recombs.df$end, range=recombs.df$width)
    
  results <- list()
  results[['segments']] <- segments.df
  results[['recombs']] <- recombs.df
  return(results)
}
  
    