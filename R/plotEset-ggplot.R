#' visualize and \link{esetPlot} with the 'ggplot2' package
#' @param object object of class \link{esetPlot}
#' @return \code{ggplot} object
#' @importFrom grid unit
#' @importFrom utils packageVersion
#' @author Laure Cougnaud
ggPlotEset <- function(object){
		
	if(!requireNamespace("ggplot2", quietly = TRUE))
		stop(paste("The package 'ggplot2' need to be loaded to create static plots."))
	
	
	plotGeneSets <- if(length(object@geneSets) > 0 & object@topGeneSets > 0){
		if(length(object@geneSetsVar) == 0){
			warning(paste("No gene sets are plotted because the variable",
				"describing the mapping of genes IDs between the 'geneSets' and the 'eset'",
				"objects is not specified in the 'geneSetsVar' parameter."
			))
			FALSE
		}else TRUE
	}else FALSE

	# bind data samples with annotation
	dataPlotSamplesWithAnnotation <- getDataPlotSamplesWithAnnotation(object)
	
	# base plot
	g <- ggplot2::ggplot()
	
	## gene plot first
	if(object@cloudGenes & nrow(object@dataPlotGenes) > 0){
				
		# implementation for fill
		baseFillColor <- do.call("rgb", c(as.list(
			colorRamp(c("white", object@cloudGenesColor))(0.2)), list(maxColorValue = 255))
		)
		
		# because parameter/function names changed in ggplot2
		if (utils::packageVersion("ggplot2") < "2.0.0"){
			fillAes <- '..count..'
			statBinHexFct <- ggplot2::stat_binhex
		}else{
			# ggplot2 2.0.0: stat_binhex() has been renamed to stat_bin_hex()
			statBinHexFct <- ggplot2::stat_bin_hex
			# ggplot2 2.1.0: output of stat_bin_hex() is now value instead of count.
			fillAes <- if (utils::packageVersion("ggplot2") == "2.1.0")	'..value..'	else	'..count..'
		}
		
		g <- g + statBinHexFct(
			ggplot2::aes_string(x = 'X', y = 'Y', fill = fillAes),
				data = object@dataPlotGenes, bins = object@cloudGenesNBins) +
			ggplot2::scale_fill_gradientn(colours = c(baseFillColor, object@cloudGenesColor))#c(0.5, 0.8)
		
		#if(includeLegend)	nLegendsSamples <- length(c(colorVar, shapeVar, sizeVar, alphaVar))
		# set legend if specified, and in the first position
		fillGuideArgs <- if(!object@cloudGenesIncludeLegend){
			ifelse(packageVersion("ggplot2") >= "3.3.4", "none", FALSE)
		}else{
			ggplot2::guide_legend(
				order = 1, 
				# change legend title, by default 'count'
				title = object@cloudGenesTitleLegend
			)
		}
		g <- g + ggplot2::guides(fill = fillGuideArgs)
			
		
	}
	
	## sample plot
	
	# fix for ggplot when variable names contained space
	formatVariableSpace <- function(var)
		if(grepl(" ", var))	paste0("`", var, "`")	else var
	
	aesArgSamplePlot <- c(list(x = 'X', y = 'Y'),
		if(length(object@colorVar) > 0)	list(color = formatVariableSpace(object@colorVar)),
		if(length(object@shapeVar) > 0)	list(shape = formatVariableSpace(object@shapeVar)),
		if(length(object@sizeVar) > 0)	list(size = formatVariableSpace(object@sizeVar)),
		if(length(object@alphaVar) > 0)	list(alpha = formatVariableSpace(object@alphaVar))
	)
	
	geomPointArgs <- c(
		list(data = dataPlotSamplesWithAnnotation,
			mapping = do.call(ggplot2::aes_string, aesArgSamplePlot)), 
		if(setFixElement(object@colorVar, object@color))	list(color = object@color),
		if(setFixElement(object@shapeVar, object@shape))	list(shape = object@shape),
		if(setFixElement(object@sizeVar, object@size))	list(size = object@size),
		if(setFixElement(object@alphaVar, object@alpha))	list(alpha = object@alpha)
	)
	
	g <- g + do.call(ggplot2::geom_point, geomPointArgs)
	
	#horizontal/vertical lines
	if(object@includeLineOrigin){
		g <- g + ggplot2::geom_vline(xintercept = 0, linetype = 'dashed')			  # X axis
		g <- g + ggplot2::geom_hline(yintercept = 0, linetype = 'dashed')			  # Y axis
	}
	
	#manual specifications: custom scales
	#only if variable, values are specified and if the variable is not numeric or integer (doesn't work with ggplot2)
	setManualScaleStatic <- function(typeVar, nameVar, valVar){
		values <- formatManualScale(x = dataPlotSamplesWithAnnotation, valVar, nameVar)
		do.call(getFromNamespace(paste("scale", typeVar, "manual", sep = "_"), ns = "ggplot2"),
			list(values = values) )#, name = nameVar, 
	}
	if (setManualScale(dataPlotSamplesWithAnnotation, object@colorVar, object@color))	
		g <- g + setManualScaleStatic(typeVar = "color", nameVar = object@colorVar,
			valVar = object@color)
	
	if (setManualScale(dataPlotSamplesWithAnnotation, object@shapeVar, object@shape))	
		g <- g + setManualScaleStatic("shape", object@shapeVar, object@shape)
	
	if (setManualScale(dataPlotSamplesWithAnnotation, object@sizeVar, object@size))	
		g <- g + setManualScaleStatic("size", object@sizeVar, object@size)
	
	if (setManualScale(dataPlotSamplesWithAnnotation, object@alphaVar, object@alpha))	
		g <- g + setManualScaleStatic("alpha", object@alphaVar, object@alpha)
	
	## axes 
	
	# plot axes labels
	if(object@xlab != "")	g <- g + ggplot2::xlab(paste0("\n", object@xlab))
	if(object@ylab != "")	g <- g + ggplot2::ylab(paste(object@ylab, "\n"))
	
	# increase margin between ticks and axes labels
	getElementText <- function(...)
		if (utils::packageVersion("ggplot2") < "2.0.0")
			ggplot2::element_text(...)	else
			ggplot2::element_text(..., margin = unit(0.5, "lines"))
	
	argsTheme <- c(
		list(
			# labels y axes centered on ticks
			axis.text.y = getElementText(vjust = 5),
			# title x axis further away from the axes labels
			axis.title.x = ggplot2::element_text(vjust = -1)
		),
		# code for previous ggplot2 version
		if (utils::packageVersion("ggplot2") < "2.0.0"){
			list(axis.ticks.margin = unit(0.5, "lines"))
			# code for more recent version
			}else{
				list(axis.text.x = ggplot2::element_text(
					margin = unit(0.5, "lines")))
			}
	)
	g <- g + do.call(ggplot2::theme, argsTheme)
	
	# symmetry
	if(object@symmetryAxes != "none"){
		
		# define axes limits
		axesLimits <- getAxesLimits(object)
		# set axes limits
		g <- g + 
			ggplot2::scale_x_continuous(limits = axesLimits[, "x"]) +
			ggplot2::scale_y_continuous(limits = axesLimits[, "y"])
		
	}
	
	# custom size range, works only if size variable is numeric, or integer
	if(length(object@sizeVar) > 0 &&
		class(dataPlotSamplesWithAnnotation[, object@sizeVar]) %in% c("numeric", "integer") & 
		length(object@sizeRange) > 0)	
		g <- g + ggplot2::scale_size(range = object@sizeRange)
	
	# custom transparency range, works only if alpha variable is numeric, or integer
	if(length(object@alphaVar) > 0 &&
		class(dataPlotSamplesWithAnnotation[, object@alphaVar]) %in% c("numeric", "integer")
			& length(object@alphaRange) > 0)	
		g <- g + ggplot2::scale_alpha(range = object@alphaRange)
	
	# add title
	if((is.character(object@title) && object@title != "") | (is.expression(object@title) & length(object@title) > 0))	
		g <- g + ggplot2::ggtitle(object@title) +
			# in version ggplot2 == 2.2.0 title is left adjusted by default
			ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
	
	
	# theme
	g <- g + ggplot2::theme_bw()
	
	## annotation top genes/samples, at the end to avoid overlapping with plot
	
	topElements <- NULL
	
	# top genes
	if(object@topGenes > 0 & nrow(object@dataPlotGenes) > 0){
		outputTopGenes <- plotTopElements(
			top = object@topGenes, type = "gene", 
			var = object@topGenesVar, cex = object@topGenesCex, 
			just = object@topGenesJust,
			color = object@topGenesColor,
			dataPlotGenes = object@dataPlotGenes, 
			esetUsed =  object@eset,
			returnTopElements = object@returnTopElements,
			packageTextLabel = object@packageTextLabel
		)
		if(object@returnTopElements){
			topElements <- c(topElements, list(topGenes = outputTopGenes$topElements))
			if(!is.null(outputTopGenes$geomText)){
				g <- g + outputTopGenes$geomText
			}else message("No labels are available for the top genes to annotate so they are not represented in the plot.")
		}else{
			if(!is.null(outputTopGenes))	g <- g + outputTopGenes	else
				message("No labels are available for the top genes to annotate so they are not represented in the plot.")
		}
	}
	
	# top samples
	if(object@topSamples > 0){
		outputTopSamples <- plotTopElements(top = object@topSamples, type = "sample", 
			var = object@topSamplesVar, cex = object@topSamplesCex, just = object@topSamplesJust,
			color = object@topSamplesColor,
			dataPlotSamples = dataPlotSamplesWithAnnotation, esetUsed = object@eset,
			returnTopElements = object@returnTopElements,
			packageTextLabel = object@packageTextLabel)
		if(object@returnTopElements){
			topElements <- c(topElements, list(topSamples = outputTopSamples$topElements))
			g <- g + outputTopSamples$geomText
		}else	g <- g + outputTopSamples
	}
	
	#top gene sets
	if(plotGeneSets){
		outputTopGeneSets <- plotTopElements(
			top = object@topGeneSets, type = "geneSets", 
			cex = object@topGeneSetsCex, just = object@topGeneSetsJust,
			color = object@topGeneSetsColor,
			geneSets = object@geneSets, esetUsed = object@eset, 
			dataPlotGenes = object@dataPlotGenes,
			geneSetsVar = object@geneSetsVar,
			geneSetsMaxNChar = object@geneSetsMaxNChar,
			returnTopElements = object@returnTopElements,
			packageTextLabel = object@packageTextLabel)
		if(object@returnTopElements){
			topElements <- c(topElements, list(topGeneSets = outputTopGeneSets$topElements))
			g <- g + outputTopGeneSets$geomText
		}else	g <- g + outputTopGeneSets
	}
	
	if(!object@includeLegend)	g <- g + ggplot2::theme(legend.position = "none")
	
	#if(!is.null(nColLegend))	guide_legend(ncol = nColLegend)
	
	res <- if(object@returnTopElements & !is.null(topElements)){
		list(topElements = topElements, plot = g)
	}else g
	
	return(res)
	
}
