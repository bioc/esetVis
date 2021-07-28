#' visualize and \link{esetPlot} with the the 'ggvis' package
#' @param object object of class \link{esetPlot}
#' @return \code{ggvis} plot object
#' @author Laure Cougnaud
ggvisPlotEset <- function(object){
	
	esetMethods <- getMethodsInputObjectEsetVis(object@eset)
	
	if(!requireNamespace("ggvis", quietly = TRUE))
		stop(paste("The package 'ggvis' need to be loaded to create",
			"interactive plots with ggvis."))
	
	# bind data samples with annotation
	dataPlotSamplesWithAnnotation <- getDataPlotSamplesWithAnnotation(object)

	includeCloudGenes <- object@cloudGenes & length(object@dataPlotGenes) > 0
	
	# ggvis doesn't allow different variable type for same aes, here transparency of cloud of genes is numeric and sample can be a factor
	alphaVar <- if(includeCloudGenes && length(object@alphaVar) > 0 && !is.numeric(dataPlotSamplesWithAnnotation[, object@alphaVar])){
		warning("A factor variable for transparency cannot be used while cloud of genes is required so 'shapeVar' is ignored.")
		character(0)
	}else	object@alphaVar
	
	## sample plot
	getProp <- function(type, typeVar)
		if(length(typeVar) > 0)	list(ggvis::prop(type, as.name(typeVar)))
	
	ggvisArgsSamplePlot <- c(
		list(data = dataPlotSamplesWithAnnotation),
		getProp("x", "X"),
		getProp("y", "Y"),
		getProp("fill", object@colorVar),
		getProp("shape", object@shapeVar),
		getProp("size", object@sizeVar),
		getProp("opacity", alphaVar)
		#if(interactiveTooltip)	list(props(key := ~"keyggvis"))
	)
	
	## gene plot first
	if(includeCloudGenes){
		
		hexbinGeneData <- hexbin(
			object@dataPlotGenes$X, 
			object@dataPlotGenes$Y, 
			xbins = object@cloudGenesNBins
		)
		hexbinGeneDataDf <- data.frame(
			xGene = hexbinGeneData@xcm, 
			yGene = hexbinGeneData@ycm, 
			geneCol = hexbinGeneData@count
		)
		
		g <- ggvis::ggvis(x = ~xGene, y= ~yGene, data = hexbinGeneDataDf)
		g <- ggvis::layer_points(vis = g, fillOpacity = ~geneCol, fill = object@cloudGenesColor)
		g <- ggvis::hide_legend(vis = g, c("fillOpacity", "fill")) 
		g <- ggvis::scale_numeric(vis = g, property = "fillOpacity", trans = "sqrt")
		
		# add sample plot to the gene plot
		ggvisArgsSamplePlotWithGenePlot <- c(list(vis = g), ggvisArgsSamplePlot)
		g <- do.call(ggvis::layer_points, ggvisArgsSamplePlotWithGenePlot)
		
	}else	g <- ggvis::layer_points(vis = do.call(ggvis::ggvis, ggvisArgsSamplePlot))
	
	#TODO: horizontal/vertical lines
	#			g <- g + geom_vline(xintercept = 0, linetype = 'dashed')			  # X axis
	#			g <- g + geom_hline(yintercept = 0, linetype = 'dashed')			  # Y axis
	
	# add fixed elements
	if(setFixElement(object@colorVar, object@color))
		g <- ggvis::layer_points(g, fill := object@color)
	if(setFixElement(object@shapeVar, object@shape))
		g <- ggvis::layer_points(g, shape := object@shape)
	if(setFixElement(object@sizeVar, object@size))
		g <- ggvis::layer_points(g, size := object@size)
	if(setFixElement(alphaVar, object@alpha))
		g <- ggvis::layer_points(g, opacity := object@alpha)

	setManualScaleGgvis <- function(g, typeVar, nameVar, valVar, typeScale = "logical"){
		values <- if(length(valVar) > 0)	
			formatManualScale(x = dataPlotSamplesWithAnnotation, valVar, nameVar)
		do.call(
			getFromNamespace(paste0("scale_", typeScale), ns = "ggvis"), 
			list(vis = g, property = typeVar, range = values)
		)
	}	
	

	
	# manual specifications: custom scales
	if (setManualScale(dataPlotSamplesWithAnnotation, object@colorVar, object@color))	
		g <- setManualScaleGgvis(g, typeVar = "fill", nameVar = object@colorVar, valVar = object@color)
	if (setManualScale(dataPlotSamplesWithAnnotation, object@shapeVar, object@shape))	
		g <- setManualScaleGgvis(g, typeVar = "shape", nameVar = object@shapeVar, valVar = object@shape)
	if (setManualScale(dataPlotSamplesWithAnnotation, object@sizeVar, object@size))	
		g <- setManualScaleGgvis(g, typeVar = "size", nameVar = object@sizeVar, valVar = object@size)
	
	# need to specify range for transparence
	alphaRange <- if(length(object@alphaRange) == 0){
		c(0.1, 1)
	}else	object@alphaRange
	if (length(alphaVar) > 0)
		g <- ggvis::scale_ordinal(vis = g, property = "opacity", range = alphaRange)
	
	#plot axes labels and title
	g <- ggvis::add_axis(vis = g, "x", title = object@xlab)
	g <- ggvis::add_axis(vis = g, "y", title = object@ylab)
	
	#symmetry
	if(object@symmetryAxes != "none"){
		
		#define axes limits
		axesLimits <- getAxesLimits(object)
		#set axes limits
		g <- ggvis::scale_numeric(vis = g, "x", trans = "linear", domain = axesLimits[, "x"])
		g <- ggvis::scale_numeric(vis = g, "y", trans = "linear", domain = axesLimits[, "y"])
		
	}
	
	#TODO: custom size range, works only if size variable is numeric, or integer
	#			if(class(dataPlotWithAnnotation[, sizeVar]) %in% c("numeric", "integer")
	#				& !is.null(sizeRange))	g <- g + scale_size(range = sizeRange)
	
	#add title
	if (object@title != ""){
		g <- ggvis::add_axis(
			vis = g, "x", orient = "top", ticks = 0, 
			title = object@title,
			properties = ggvis::axis_props(
				axis = list(stroke = "white"),
				labels = list(fontSize = 2)
			)
		)
	}
	
	figInteractiveSize <- if(length(object@sizePlot) != 2){
		warning(paste("The size of the ggvis window should contain two",
			"elements: width and height. The default values c(600, 400) are used."))
		c(600, 400)
		# default options
		# str(default_options())
	}else object@sizePlot
	
	g <- ggvis::set_options(vis = g, 
		width = figInteractiveSize[1], 
		height = figInteractiveSize[2]
	)
	
#	orderLegendLog <- c(!is.null(colorVar), !is.null(shapeVar), !is.null(sizeVar))
#	names(orderLegendLog) <- c("fill", "shape", "size")
	#			if(sum(orderLegendLog) > 0){
	#				legendToSet <- names(orderLegendLog)[orderLegendLog]
	#				argsAddLegend <- c(as.list(legendToSet), list(orient = "right"))
	#				g <- g %>% do.call("add_legend", argsAddLegend)
	#			}
	
	# adjust legend manually because overlap if several present
	# only if the window size is numeric (not auto)
	# this will be fixed in future version of ggvis? 
	
	orderLegendLog <- c(
		length(object@colorVar) > 0, 
		length(object@shapeVar) > 0, 
		length(object@sizeVar) > 0
	)
	names(orderLegendLog) <- c("fill", "shape", "size")
	legendToSet <- names(orderLegendLog)[orderLegendLog]
	
	if(is.numeric(figInteractiveSize) & 
		object@includeLegend & 
		object@adjustLegend & 
		length(legendToSet) > 0){
		
		# for legend side by side in y direction
		setLegendPos <- function(vis, typeVar, y)
			ggvis::add_legend(vis = vis, scales = typeVar, 
				properties = ggvis::legend_props(legend = list(y = y)))
		
		for(i in 1:length(legendToSet))
			g <- setLegendPos(vis = g, typeVar = legendToSet[i],
				y = (i-1) * figInteractiveSize[2]/length(legendToSet))
		
		# for legend side by side in x direction
		#					setLegendPos <- function(vis, typeVar, x)
		#						add_legend(vis = vis, scales = typeVar, 
		#								properties = legend_props(legend = list(x = x)))
		#					
		#					for(i in 1:length(legendToSet))
		#						g <- setLegendPos(g, typeVar = legendToSet[i],
		#							x = (0.1*(i-1) + 1) * figInteractiveSize[1])
		#				
	}
	
	
	# hide legend
	if((!object@includeLegend) & length(legendToSet) > 0)	
		g <- ggvis::hide_legend(vis = g, scales = legendToSet)
	
	if(object@includeTooltip){
		
		tooltipVars <- object@tooltipVars[
			object@tooltipVars %in% esetMethods$varLabels(object@eset)
		]
		
		if(length(tooltipVars) == 0)
			tooltipVars <- NULL
		
		tooltipLabelsFct <- function(x)
			return(
					
				if(!is.null(x) && !any(c("xGene", "yGene") %in% names(x))){
						
					formatCoor <- function(x)
						paste0("(", paste(formatC(x, width = 3), collapse = ", "), ")")
					
					xy <- as.numeric(x[, c("X", "Y")])
					xyST <- formatCoor(xy)
					
					# use bracket rather than subset ta avoid having notes in R CMD check
					dataPoints <- dataPlotSamplesWithAnnotation[
						round(dataPlotSamplesWithAnnotation$X, 3) == round(xy[1], 3) & 
							round(dataPlotSamplesWithAnnotation$Y, 3) == round(xy[2], 3), 
					]
					
					if(nrow(dataPoints) > 0 & length(tooltipVars) > 0){
						addedAnnot <- esetMethods$pData(object@eset)[dataPoints[1, "sampleName"], 
							tooltipVars, drop = FALSE]
						colnames(addedAnnot) <- tooltipVars
						x <- data.frame(x, addedAnnot)
					}
					
					extraCols <- colnames(x)[!colnames(x) %in% c("X", "Y")]
					
					if(length(extraCols) > 0){
						sampleAnnot <- paste0(names(x[, extraCols]), ": ", x[, extraCols], collapse = "<br />")
						paste0(xyST, "<br />", sampleAnnot)
					}else	xyST
				}
			)
		
			g <- ggvis::add_tooltip(vis = g, tooltipLabelsFct, "hover") # can use 'click' too
		
	}
		
	## TODO: annotation top genes/samples, at the end to avoid overlapping with plot
	
	# to keep legend separated
	g <- ggvis::set_options(vis = g, duration = 0)
	
	return(g)
	
}
