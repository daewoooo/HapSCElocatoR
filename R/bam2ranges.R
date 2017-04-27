#' Import BAM file into GRanges
#'
#' Import aligned reads from a BAM file into a \code{\link{GRanges}} object.
#'
#' @param file Bamfile with aligned reads.
#' @param bamindex Bam-index file with or without the .bai ending. If this file does not exist it will be created and a warning is issued.
#' @param pairedEndReads Set to \code{TRUE} if you have paired-end reads in your file.
#' @param min.mapq Minimum mapping quality when importing from BAM files.
#' @importFrom Rsamtools indexBam scanBamHeader ScanBamParam scanBamFlag
#' @importFrom GenomicAlignments readGAlignmentPairs readGAlignments first last
#' @author David Porubsky
#' @export

bam2ranges <- function(file, bamindex=file, min.mapq=10, pairedEndReads=FALSE, chromosomes=NULL) {
  ## Check if bamindex exists
  bamindex.raw <- sub('\\.bai$', '', bamindex)
  bamindex <- paste0(bamindex.raw,'.bai')
  if (!file.exists(bamindex)) {
    bamindex.own <- Rsamtools::indexBam(file)
    warning("Couldn't find BAM index-file ",bamindex,". Creating our own file ",bamindex.own," instead.")
    bamindex <- bamindex.own
  }
  file.header <- Rsamtools::scanBamHeader(file)[[1]]
  chrom.lengths <- file.header$targets
  chroms.in.data <- names(chrom.lengths)
  if (is.null(chromosomes)) {
    chromosomes <- chroms.in.data
  }
  chroms2use <- intersect(chromosomes, chroms.in.data)
  if (length(chroms2use)==0) {
    chrstring <- paste0(chromosomes, collapse=', ')
    stop('The specified chromosomes ', chrstring, ' do not exist in the data. Please try ', paste(paste0('chr',chromosomes), collapse=', '), ' instead.')
  }
  ## Issue warning for non-existent chromosomes
  diff <- setdiff(chromosomes, chroms.in.data)
  if (length(diff)>0) {
    diffs <- paste0(diff, collapse=', ')
    warning(paste0('Not using chromosomes ', diffs, ' because they are not in the data.'))
  }
  ## Import the file into GRanges
  gr <- GenomicRanges::GRanges(seqnames=chroms2use, ranges=IRanges(start=rep(1, length(chroms2use)), end=chrom.lengths[chroms2use]))

  if (pairedEndReads) {
      data.raw <- GenomicAlignments::readGAlignmentPairs(file, index=bamindex, param=Rsamtools::ScanBamParam(tag="XA", which=range(gr), what='mapq', flag=scanBamFlag(isDuplicate=F)))
  } else {
      data.raw <- GenomicAlignments::readGAlignments(file, index=bamindex, param=Rsamtools::ScanBamParam(tag="XA", which=range(gr), what='mapq', flag=scanBamFlag(isDuplicate=F)))
  }
  
  
  ## Second mate of the pair will inherit directionality from the first mate of the pair
  if (pairedEndReads) {
    data.first <- as(GenomicAlignments::first(data.raw), 'GRanges')
    data.last <- as(GenomicAlignments::last(data.raw), 'GRanges')
    strand(data.last) <- strand(data.first)
    data <- sort(c(data.first, data.last))
  } else {
    data <- as(data.raw, 'GRanges')
  }
  
  ## Filter by mapping quality
  if (!is.null(min.mapq)) {
    if (any(is.na(mcols(data)$mapq))) {
      warning(paste0(file,": Reads with mapping quality NA (=255 in BAM file) found and removed. Set 'min.mapq=NULL' to keep all reads."))
      mcols(data)$mapq[is.na(mcols(data)$mapq)] <- -1
    }
    data <- data[mcols(data)$mapq >= min.mapq]
  }
  seqlevels(data) <- seqlevels(gr)
  return(data)
}

