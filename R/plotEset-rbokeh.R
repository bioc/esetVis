#' visualize and \link{esetPlot} with the the 'rbokeh' package
#' @param object object of class \link{esetPlot}
#' @return \code{rbokeh} plot
#' @author Laure Cougnaud
rbokehPlotEset <- function(object){

	esetMethods <- getMethodsInputObjectEsetVis(object@eset)
	
	if(!requireNamespace("rbokeh", quietly = TRUE))
		stop(paste("The package 'rbokeh' need to be loaded to create",
			"interactive plots with rbokeh."))
	
	if(length(object@alphaVar) > 0 && is.factor(esetMethods$pData(object@eset)[, object@alphaVar])){
		warning("The transparency aesthetic (alpha) is not ",
			"yet implemented for factors in rbokeh interactive plot, ",
			"so no transparency is used.")
		alphaVar <- NULL; alpha <- 1
	}else{alphaVar <- object@alphaVar; alpha <- object@alpha}
	
	# bind data samples with annotation
	dataPlotSamplesWithAnnotation <- getDataPlotSamplesWithAnnotation(object)
	
	#define axes limits
	axesLimits <- getAxesLimits(object)
	
	## create empty figure
	argsFigure <- c(
		if(is.numeric(object@sizePlot)) 
			list(width = object@sizePlot[1], 
				height = object@sizePlot[2]
			),
		if(object@title != "") list(title = object@title),
		if(object@xlab != "")	list(xlab = object@xlab),
		if(object@ylab != "")	list(ylab = object@ylab),
		if(length(axesLimits) > 0)	
			list(xlim = axesLimits[, "x"], ylim = axesLimits[, "y"]),
		list(xgrid = FALSE, ygrid = FALSE)
	)
	g <- do.call(rbokeh::figure, argsFigure)
	
	# add axes limits
	
	## gene plot first
	if(object@cloudGenes & nrow(object@dataPlotGenes) > 0){
		
		# can only specify color as palette for color ramp: don't use white at lower palette
		paletteFctUsed <- colorRampPalette(
			c(colorRampPalette(c("white", object@cloudGenesColor))(10)[2], object@cloudGenesColor))
		hoverGenes <- object@dataPlotGenes
		hoverGenes[, c("X", "Y")] <- round(hoverGenes[, c("X", "Y")], digits = 2)
		g <- rbokeh::ly_hexbin(
			fig = g, 
			x = object@dataPlotGenes$X, 
			y = object@dataPlotGenes$Y,
			xbins = object@cloudGenesNBins, 
			palette = paletteFctUsed, trans = sqrt,
			hover = hoverGenes
		)#, alpha = 0.8
		
	}
	
	## samples plot
	
	# need to remove NA values for glyph?
	varsPlot <- c(object@colorVar, object@shapeVar, object@sizeVar)
	if(length(varsPlot) > 0){
		idxRowsKept <- rowSums(is.na(dataPlotSamplesWithAnnotation[, varsPlot, drop = FALSE])) == 0
		dataPlotWithAnnotationWthtNA <- dataPlotSamplesWithAnnotation[idxRowsKept, ]
	}else	dataPlotSamplesWithAnnotation
		
	# possible to specify a data.frame for the hoover, 
	# so add additional variables if requested
	
	# remove variable if already present in the data
	# (so should have the same name several times)
	tooltipVars <- if(length(object@tooltipVars) > 0){
		object@tooltipVars[
			!object@tooltipVars %in%  colnames(dataPlotWithAnnotationWthtNA)]
	}else character()
	
	hoverDf <- if(object@includeTooltip){
		hoverDf <- dataPlotWithAnnotationWthtNA
		hoverDf[, c("X", "Y")] <- round(hoverDf[, c("X", "Y")], digits = 2)
		if(length(tooltipVars) > 0){
			hoverDf <- cbind.data.frame(hoverDf, 
				esetMethods$pData(object@eset)[as.character(hoverDf$sampleName), tooltipVars]
			)
			colnames(hoverDf) <- c(colnames(dataPlotWithAnnotationWthtNA), tooltipVars)
		}
		# remove redundant column
		hoverDf <- as.data.frame(t(unique(t(hoverDf))))
	}else	NULL
	
	# samples plot
	
	color <- if(length(object@colorVar) > 0)	object@colorVar	else	object@color
	glyph <- if(length(object@shapeVar) > 0)	object@shapeVar	else	object@shape
	size <- if(length(object@sizeVar) > 0)	object@sizeVar	else	object@size
	alpha <- if(length(alphaVar) > 0)	alphaVar	else	alpha
	g <- rbokeh::ly_points(
		fig = g,
		data = dataPlotWithAnnotationWthtNA,
		x = "X", y = "Y",
		hover = hoverDf,
		#lname = "sample",
		legend = object@includeLegend,
		color = color,
		glyph = glyph,
		size = size,
		alpha = alpha
	)
	
	## add horizontal/vertical lines
	
	if(object@includeLineOrigin)
		g <- rbokeh::ly_abline(fig = g, v = 0, type = 2)
	g <- rbokeh::ly_abline(fig = g, h = 0, type = 2)
	
	# TODO: add custom palettes when will be available in rbokeh
	
	return(g)
	
}