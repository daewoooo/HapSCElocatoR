#' Wrapper function for the YeastSCElocatoR package
#'
#' @param bamfilepath Folder with BAM files to analyze
#' @param dataDirectory Folder to output the results. If it does not exist it will be created.
#' @param method

#' @inheritParams bam2ranges
#' @inheritParams findSegmentsCBC
#' @importFrom gtools mixedsort

#' @author David Porubsky
#' @export

localizeSCE <- function(bamfilepath, dataDirectory="./SCE_analysis", chromosomes=NULL, min.mapq=10, minSeg=0, minSize=0, smooth=1, method="CBC", hotspots=FALSE, read.len=50) {
  
  if (!file.exists(dataDirectory)) {
    dir.create(dataDirectory)  
  }
  
  results <- file.path(dataDirectory, 'Results')
  if (!file.exists(results)) {
    dir.create(results)  
  }
  
  browserpath <- file.path(dataDirectory, 'browserFiles')
  if (!file.exists(browserpath)) {
    dir.create(browserpath)  
  }
  
  #plots <- file.path(dataDirectory, 'Plots')
  #if (!file.exists(plots)) {
  #  dir.create(plots)  
  #}
  
  bamfiles <- list.files(bamfilepath, full.names=TRUE, pattern = paste0('.bam$'))
  
  plots <- list()
  for (bam in bamfiles) {
      message("Working on ", bam)
      filename <- basename(bam)
      prefix <- gsub("\\.bam", "", filename)
      fragments <- bam2ranges(file=bam, bamindex=bam, min.mapq=min.mapq, chromosomes=chromosomes)
      #seqlevels(fragments) <- sub(pattern='Mito', replacement="M", seqlevels(fragments))
      
      exportUCSC(index=filename, outputDirectory=browserpath, fragments=fragments)
      
      #Translate read directionality into a binary vector (0|1)
      #frag.df <- as(fragments, "data.frame")
      #frag.df$strand <- as.character(frag.df$strand)
      #frag.df$strand[frag.df$strand == "+"] <- 0
      #frag.df$strand[frag.df$strand == "-"] <- 1
      
      #frag.binary <- split(frag.df, frag.df$seqnames)
      
      frag.grl <- split(fragments, seqnames(fragments))
      
      bam.segm <- findSegmentsCBC(data=frag.grl, minSeg=minSeg, minSize=minSize, smooth=smooth, read.len=read.len)

      if (length(bam.segm$segments) > 0) {
        segments.df <- as(bam.segm$segments, "data.frame")
        names(segments.df)[1] <- "chromosome"
        #segments.df <- data.frame(chromosome=segments.df$index, start=segments.df$start, end=segments.df$end, num.read=segments.df$num.mark, seg.mean=segments.df$seg.mean, strand=segments.df$direction)
        recombs.df <- as(bam.segm$recombs, "data.frame")
        names(recombs.df)[1] <- "chromosome"
        #recombs.df <- data.frame(chromosome=recombs.df$index, start=recombs.df$start, end=recombs.df$end, range=recombs.df$width, genoT=recombs.df$genoT)
        exportUCSC(index=filename, outputDirectory=browserpath, segments=segments.df)
        exportUCSC(index=filename, outputDirectory=browserpath, recombs=recombs.df)
      
        #destination <- file.path(results, paste0(prefix, "_segments.txt"))
        #write.table(bam.segm$segments, file=destination, row.names=F, quote=F)
      
        #destination <- file.path(results, paste0(prefix, "_breaks.txt"))
        #write.table(bam.segm$recombs, file=destination, row.names=F, quote=F)
      
        results.data <- list(reads=fragments, segments=bam.segm$segments, recombs=bam.segm$recombs)
        destination <- file.path(results, paste0(prefix, ".RData"))
        save(results.data, file=destination)

        # Plotting resulting segments	
        #segments <- bam.segm$segments
        segments.df$chromosome <- factor(segments.df$chromosome, levels=gtools::mixedsort(as.character(unique(segments.df$chromosome))))
        my_theme <- theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())  
      
        plt <- ggplot(segments.df) + geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=1, fill=direction)) + facet_grid(chromosome ~ .) + scale_fill_manual(values = c("paleturquoise4", "sandybrown", "olivedrab")) + ggtitle(filename) + theme_bw() + my_theme
        plots[[length(plots)+1]] <- plt
      } else {
        message("Insufficient read number, Skipping!!! ")
      }
      
  }
  
  if (hotspots) {
    message("Calling SCE hotspots ...")
    files <- list.files(results, pattern=".RData$", full.names=TRUE)
    breaks.all.files <- GenomicRanges::GRangesList()
    for (file in files) {
      breakpoints <- get(load(file))$recombs
      if (length(breakpoints)) {
        suppressWarnings( breaks.all.files[[file]] <- breakpoints ) #TODO check if this can be done without warnings
      }  
    }
    ## write a bedgraph of data (overlapping breaks)
    
    ## Search for SCE hotspots
    hotspots <- hotspotter(gr.list=breaks.all.files, bw = 10000, pval = 0)
    hotspots <- hotspots[hotspots$num.events > floor(length(breaks.all.files)*0.1)]
  
    breakpoints <- unlist(breaks.all.files)
    hits <- findOverlaps(breakpoints, hotspots)
    breaks.overlap <- breakpoints[queryHits(hits)]
    breaks.grl <- split(breaks.overlap, subjectHits(hits))
  
    hotspot.ranges <- GenomicRanges::GRangesList()
    for (i in 1:length(breaks.grl)) {
      gr <- breaks.grl[[i]]
      outer.range <- reduce(gr)
    
      if (length(outer.range) > 1) {
        outer.range <- GRanges(seqnames=seqnames(outer.range)[1], ranges=IRanges(start=start(outer.range)[1], end=end(outer.range)[length(outer.range)]))
      }
    
      outer.range$range <- "outer"
      outer.range$ID <- i
      ranges.disjoin <- disjoin(gr) # redefine ranges in df
      hits <- countOverlaps(ranges.disjoin, gr) # counts number of breaks overlapping at each range
      mcols(ranges.disjoin)$hits <- hits 
      max.cov <- ranges.disjoin[ranges.disjoin$hits == max(hits)]
    
      if (length(max.cov) > 1) {
        inner.range <- GRanges(seqnames=seqnames(max.cov)[1], ranges=IRanges(start=start(max.cov)[1], end=end(max.cov)[length(max.cov)]))
      } else {
        inner.range <- max.cov[,0]
      }
      
      inner.range$range <- "inner"
      inner.range$ID <- i
      refined.ranges <- c(inner.range, outer.range) 
      hotspot.ranges[[i]] <- refined.ranges
    }
    hotspot.ranges <- unlist(hotspot.ranges)
    hotspot.ranges.df <- as(hotspot.ranges, "data.frame")
    destination <- file.path(dataDirectory, "hotspot_ranges.txt")
    write.table(hotspot.ranges.df, file = destination, quote = F, row.names = F)
  
    all.breaks <- unlist(breaks.all.files)
    filenames <- basename(names(all.breaks))
    names(all.breaks) <- NULL
    all.breaks.df <- as(all.breaks, 'data.frame')
    all.breaks.df$filenames <- filenames
    destination <- file.path(dataDirectory, "breakpoint_ranges.txt")
    write.table(all.breaks.df, file = destination, quote = F, row.names = F)
  
    message("Exporting heatmap ...")
    destination <- file.path(dataDirectory, "heatmap.pdf")
    plotHeatmap(datapath=results, hotspots=hotspots, plotLibraries = NULL, file=destination)
  }
  
  message("Exporting ideograms ...")
  destination <- file.path(dataDirectory, "results_plots.pdf")
  pdf(destination, width = 10, height = 6)
  bquiet = lapply(plots, print)
  d <- dev.off()
}
