#' Genome wide heatmap of template inheritance states
#'
#' Plot a genome wide heatmap of template inheritance states.
#'
#' @param datapath location of RData files storing reads and breakpoint coordinates
#' @param plotLibraries subset of libaries to be plotted, otherwise all libraries will be plotted
#' @param file name of the file to store libraries in
#' 
#' @author David Porubsky, Aaron Taudt, Ashley Sanders
#' @export


plotHeatmap <- function(datapath, hotspots=NULL, plotLibraries = NULL, file=NULL) {

	files <- list.files(datapath, pattern=".RData$", full=T)

	if (!is.null(plotLibraries)) {
		libs2plot <- plotLibraries
	} else {
		libs2plot <- c(1:length(files))
	}

	message("Preparing heatmap from ",length(libs2plot), " libraries")

	IDs <- list()
	grl <- GRangesList()
	breaks <- GRangesList()
	for (i in libs2plot) {	
		data <- get(load(files[[i]]))
		IDs[[i]] <- basename(files[[i]])
		grl[[i]] <- data$segments
		suppressWarnings( breaks[[i]] <- data$recombs )
	}

	#transform bin coordinates of each chromosome into genomewide coordinates (cumulative sum of bin coordintes)
	cum.seqlengths <- cumsum(as.numeric(seqlengths(grl[[1]])))
	cum.seqlengths.0 <- c(0,cum.seqlengths[-length(cum.seqlengths)])
	names(cum.seqlengths.0) <- seqlevels(grl[[1]])

	#get positions of ends of each chromosome to plot lones between the chromosomes
	chr.lines <- data.frame( y=cum.seqlengths[-length(cum.seqlengths)] )	

	transCoord <- function(gr) {
		gr$start.genome <- start(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
		gr$end.genome <- end(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
		return(gr)
	}
	grl <- endoapply(grl, transCoord)
	
#	seqlevels(breaks) <- seqlevels(grl)
#	seqlengths(breaks) <- seqlengths(grl)

	# disjoin overlaping breaks
#	breaks <- unlist(breaks, use.names=F)
#	ranges <- disjoin(breaks) # redefine ranges in df
#	hits <- countOverlaps(ranges, breaks) # counts number of breaks overlapping at each range
# 	mcols(ranges)$hits <- hits

#	disjoin.breaks <- transCoord(ranges)
#	dfplot.disjoin.breaks <- as.data.frame(disjoin.breaks)

	# Chromosome lines for heatmap
	label.pos <- round( cum.seqlengths.0 + 0.5 * seqlengths(grl[[1]]) )
	names(label.pos) <- gsub("chr", "", names(label.pos)) 
	df.chroms <- data.frame(y=c(0,cum.seqlengths))

	# Plot breaks summary
#	ggplt1 <- ggplot(dfplot.disjoin.breaks) + geom_rect(aes(ymin=0, ymax=hits, xmin=start.genome, xmax=end.genome), fill="red", color="red")
#	ggplt1 <- ggplt1 + geom_vline(aes_string(xintercept='y'), data=chr.lines, col='black') + scale_y_continuous(expand = c(0,0))

	# Data
	df <- list()
	for (i1 in 1:length(grl)) {
		df[[length(df)+1]] <- data.frame(start=grl[[i1]]$start.genome, end=grl[[i1]]$end.genome, seqnames=seqnames(grl[[i1]]), sample=IDs[[i1]], state=grl[[i1]]$direction)
	}
	df <- do.call(rbind, df)

	## PLOT
	ggplt2 <- ggplot(df) + geom_linerange(aes_string(ymin='start', ymax='end', x='sample', col='state'), size=5) + scale_y_continuous(breaks=label.pos, labels=names(label.pos)) + coord_flip() + scale_color_manual(values=c('c'="paleturquoise4", 'wc'="olivedrab",'w'="sandybrown",'?'="red")) + theme(panel.background=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_text(size=10))
	ggplt2 <- ggplt2 + geom_hline(aes_string(yintercept='y'), data=df.chroms, col='black')
	if (!is.null(hotspots)) { 
	  hot <- transCoord(hotspots)
	  hot$midpoint <- hot$start.genome + ((hot$end.genome - hot$start.genome)/2) 
	  ggplt2 <- ggplt2 + geom_hline(yintercept = hot$midpoint, color="red", alpha=0.5) 
	}
	
	## PRINT TO FILE
	## printing to a file or returning plot object
	if (!is.null(file)) {
		message("Printing to PDF ",file)
		height.cm <- length(libs2plot) * 0.5
		width.cm <- 200
		pdf(file, width=width.cm/2.54, height=height.cm/2.54)
		print(ggplt2)
		d <- dev.off()
	} else {
		return(ggplt2)
	}
}


