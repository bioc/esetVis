#' visualize and \link{esetPlot} with the the 'plotly' package
#' @param object object of class \link{esetPlot}
#' @return \code{plotly} plot
#' @author Laure Cougnaud
#' @importFrom stats setNames
#' @keywords internal
plotlyPlotEset <- function(object){

	esetMethods <- getMethodsInputObjectEsetVis(object@eset)

	if(!requireNamespace("plotly", quietly = TRUE))
		stop(paste("The package 'plotly' need to be loaded to create",
			"interactive plots with plotly."))

	## base plot

	figInteractiveSize <- if(length(object@sizePlot) != 2){
		warning(paste("The size of the window should contain two",
			"elements: width and height. The default values c(600, 400) are used."))
			c(600, 400)
	}else object@sizePlot
	argsPlotly <- c(
		list(width = figInteractiveSize[1], height = figInteractiveSize[2]),
		# colors/shapes should be set in plot_ly call, not add_trace
		if(length(object@colorVar) > 0 && length(object@color) > 0)
			list(colors = object@color),
		if(length(object@shapeVar) > 0 && length(object@shape) > 0)
			list(symbols = object@shape),
		if(length(object@sizeVar) > 0 && length(object@size) > 0)
			list(sizes = object@size)
	)
	fig <- do.call(plotly::plot_ly, argsPlotly)

	## add gene cloud
	if(object@cloudGenes & nrow(object@dataPlotGenes) > 0){

		# create data.frame with x/y position of hexbin + counts
		dataBin <- hexbin::hexbin(
			x = object@dataPlotGenes$X, y = object@dataPlotGenes$Y,
			xbins = object@cloudGenesNBins
		)
		dataBinDf <- do.call(cbind.data.frame, hexbin::hcell2xy(dataBin))
		cBin <- dataBin@count
		dataBinDf[, "count"] <- cBin/max(cBin)

		# colors
		baseFillColor <- do.call("rgb", c(as.list(
			colorRamp(c("white", object@cloudGenesColor))(0.2)),
			list(maxColorValue = 255))
		)
		colorscale <- list(list(0, baseFillColor), list(1, object@cloudGenesColor))

		# create hexbin plot
		fig <- plotly::add_markers(
			p = fig, data = dataBinDf, x = ~ x, y = ~ y,
			type = "scatter", showlegend = FALSE,
			hoverinfo = "none",
			marker = list(color = ~ count, colorscale = colorscale)
		)

	}

	# add samples
	dataSamples <- getDataPlotSamplesWithAnnotation(object)
	setAes <- function(var, value, type){
		isVar <- !is.null(var) && length(var) > 0
		x <- list(
			if(isVar){varToFm(var)}else{I(value)}
		)
		names(x) <- type
		x <- x[!sapply(x, is.null)]
		return(x)
	}
	argsSamplePlot <- c(
		list(
			p = fig, x = varToFm("X"), y = varToFm("Y"),
			type = "scatter", mode = "markers", showlegend = TRUE
		),
		setAes(var = object@colorVar, value = object@color, type = "color"),
		setAes(var = object@sizeVar, value = object@size, type = "size"),
		setAes(var = object@shapeVar, value = object@shape, type = "symbol"),
		setAes(var = object@alphaVar, value = object@alpha, type = "opacity")
	)
	if(object@includeTooltip & length(object@tooltipVars)){
		hovertext <- lapply(object@tooltipVars, function(var)
			paste0(var, ": ", dataSamples[, var])
		)
		dataSamples[, "hovertext"] <- do.call(what = paste,
			args = c(hovertext, list(sep = "<br>")))
		argsSamplePlot[["hovertemplate"]] <- varToFm("hovertext")
	}
	argsSamplePlot[["data"]] <- dataSamples
	fig <- do.call(plotly::add_trace, argsSamplePlot)

	topElements <- NULL

	# add top genes
	if(object@topGenes > 0 & nrow(object@dataPlotGenes) > 0){

		# get coordinates of top elements
		dataTopGenes <- getTopElements(
			top = object@topGenes, type = "gene",
			var = object@topGenesVar,
			dataPlotGenes = object@dataPlotGenes,
			esetUsed =  object@eset
		)

		# plot top elements
		fig <- plotly::add_text(
			p = fig, x = ~X, y = ~Y, text = ~labels,
			data = dataTopGenes,
			textfont = list(
				color = object@topGenesColor,
				size = object@topGenesCex
			),
			showlegend = FALSE
#			textposition = object@topGenesJust
		)

		if(object@returnTopElements){
			topGenes <- setNames(dataTopGenes$labels, rownames(dataTopGenes))
			topElements <- c(topElements, list(topGenes = topGenes))
		}

	}

	if(object@topSamples > 0 & nrow(object@dataPlotSamples) > 0){

		# get coordinates of top elements
		dataTopSamples <- getTopElements(
			top = object@topSamples, type = "sample",
			var = object@topSamplesVar,
			dataPlotSamples = object@dataPlotSamples,
			esetUsed =  object@eset
		)

		# plot top elements
		fig <- plotly::add_text(
			p = fig, x = ~X, y = ~Y, text = ~labels,
			data = dataTopSamples,
			textfont = list(
				color = object@topSamplesColor,
				size = object@topSamplesCex
			),
			showlegend = FALSE
#			textposition = object@topSamplesJust
		)

		if(object@returnTopElements){
			topSamples <- setNames(dataTopSamples$labels, rownames(dataTopSamples))
			topElements <- c(topElements, list(topSamples = topGenes))
		}

	}

	# top gene sets
	plotGeneSets <- if(length(object@geneSets) > 0 & object@topGeneSets > 0){
		if(length(object@geneSetsVar) == 0){
			warning(paste("No gene sets are plotted because the variable",
				"describing the mapping of genes IDs between the 'geneSets' and the 'eset'",
				"objects is not specified in the 'geneSetsVar' parameter."
				))
			FALSE
		}else TRUE
	}else FALSE

	if(plotGeneSets){

		# get coordinates of top elements
		dataTopGeneSets <- getTopElements(
			top = object@topGeneSets, type = "geneSets",
			dataPlotGenes = object@dataPlotGenes,
			geneSets = object@geneSets,
			esetUsed =  object@eset,
			geneSetsVar = object@geneSetsVar,
			geneSetsMaxNChar = object@geneSetsMaxNChar
		)

		# plot top elements
		fig <- plotly::add_text(
			p = fig, x = ~X, y = ~Y, text = ~labels,
			data = dataTopGeneSets,
			textfont = list(
				color = object@topGeneSetsColor,
				size = object@topGeneSetsCex
			),
			showlegend = FALSE
#			textposition = object@topSamplesJust
		)

		if(object@returnTopElements){
			topGeneSets <- setNames(dataTopGeneSets$labels, rownames(dataTopGeneSets))
			topElements <- c(topElements, list(topGeneSets = topGeneSets))
		}

	}

	## axes
	argsLayout <- list()

	# title
	if(object@title != "")	argsLayout$title$text <- object@title

	# plot axes labels
	if(object@xlab != "")	argsLayout$xaxis$title$text <- object@xlab
	if(object@ylab != "")	argsLayout$yaxis$title$text <- object@ylab

	# symmetry
	if(object@symmetryAxes != "none"){

		# define axes limits
		axesLimits <- getAxesLimits(object)
		# set axes limits
		argsLayout$xaxis$range <- axesLimits[, "x"]
		argsLayout$yaxis$range <- axesLimits[, "y"]

	}

	if(length(argsLayout) > 0){
		argsLayout[["p"]] <- fig
		fig <- do.call(plotly::layout, argsLayout)
	}

	return(fig)

}

#' Get formula for a specific variable,
#' to be used in aesthetic specification in \code{\link[plotly]{plot_ly}}.
#' @param var Character vector with variable to combine.
#' Otherwise with the '+' operator.
#' @return \code{\link[stats]{as.formula}}
#' @author Laure Cougnaud
#' @importFrom stats as.formula
#' @keywords internal
varToFm <- function(var){
	fm <- as.formula(paste0("~", paste(paste0("`", var, "`"), collapse = "+")))
	return(fm)
}