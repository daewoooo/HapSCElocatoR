#' Wrapper function for the YeastSCElocatoR package
#'
#' @param bamfilepath Folder with BAM files to analyze
#' @param dataDirectory Folder to output the results. If it does not exist it will be created.
#' @param method

#' @inheritParams bam2ranges
#' @inheritParams findSegmentsCBC

#' @author David Porubsky
#' @export

localizeSCE <- function(bamfilepath, dataDirectory="./SCE_analysis", min.mapq=10, minSeg=5, smooth=1, method="CBC") {
  
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
      fragments <- bam2ranges(file=bam, bamindex=bam, min.mapq=min.mapq)
      seqlevels(fragments) <- sub(pattern='Mito', replacement="M", seqlevels(fragments))
      
      exportUCSC(index=filename, outputDirectory=browserpath, fragments=fragments)
      
      #Translate read directionality into a binary vector (0|1)
      frag.df <- as(fragments, "data.frame")
      frag.df$strand <- as.character(frag.df$strand)
      frag.df$strand[frag.df$strand == "+"] <- 0
      frag.df$strand[frag.df$strand == "-"] <- 1
      
      frag.binary <- split(frag.df, frag.df$seqnames)
      
      if (method == "CBC") {
        bam.segm <- findSegmentsCBC(data=frag.binary, minSeg=minSeg, smooth=smooth)
      } else if (method == "HMM") {
        bam.segm <- findSegmentsHMM(data=frag.binary)
      } else {
        warning("No method submitted!!!")
      }  
      
      exportUCSC(index=filename, outputDirectory=browserpath, segments=bam.segm$segments)
      exportUCSC(index=filename, outputDirectory=browserpath, recombs=bam.segm$recombs)
      
      destination <- file.path(results, paste0(prefix, "_segments.txt"))
      write.table(bam.segm$segments, file=destination, row.names=F, quote=F)
      
      destination <- file.path(results, paste0(prefix, "_breaks.txt"))
      write.table(bam.segm$recombs, file=destination, row.names=F, quote=F)

      # Plotting resulting segments	
      segments <- bam.segm$segments
      my_theme <- theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())  
      
      plt <- ggplot(segments) + geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=1, fill=strand)) + facet_grid(chromosome ~ .) + scale_fill_manual(values = c("paleturquoise4", "sandybrown", "olivedrab")) + ggtitle(filename) + theme_bw() + my_theme
      plots[[length(plots)+1]] <- plt
      
  }
  message("Exporting plot ...")
  destination <- file.path(dataDirectory, "results_plots.pdf")
  pdf(destination, width = 10, height = 6)
  bquiet = lapply(plots, print)
  d <- dev.off()
}
