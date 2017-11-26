#' plot an \link{plotEset} object
#' @param object object of class \link{esetPlot}
#' @param returnEsetPlot logical, 
#' if TRUE return also the \link{esetPlot} object,
#' such as can be re-use for future call to \link{plotEset}
#' @return the plot object if \code{returnEsetPlot} is FALSE,
#' otherwise a list with 'plot': the plot object
#'  and 'esetPlot': the \link{esetPlot} object
#' @rdname plotEset-methods
#' @docType methods
#' @author Laure Cougnaud
#' @export
setGeneric(name = "plotEset",
	def = function(object, returnEsetPlot = FALSE) {
		standardGeneric("plotEset")
	}
)

#' @aliases plotEset,ggplotEsetPlot-method
#' @rdname plotEset-methods
setMethod("plotEset", "ggplotEsetPlot", function(object, returnEsetPlot){		
	plot <- ggPlotEset(object)
	res <- formatOutput(plot, object, "ggplot", returnEsetPlot)
})

#' @aliases plotEset,ggvisEsetPlot-method
#' @rdname plotEset-methods
setMethod("plotEset", "ggvisEsetPlot", function(object, returnEsetPlot){		
	plot <- ggvisPlotEset(object)
	res <- formatOutput(plot, object, "ggvis", returnEsetPlot)
})

#' @aliases plotEset,rbokehEsetPlot-method
#' @rdname plotEset-methods
setMethod("plotEset", "rbokehEsetPlot", function(object, returnEsetPlot){		
	plot <- rbokehPlotEset(object)
	res <- formatOutput(plot, object, "rbokeh", returnEsetPlot)
})

#' generic for get axes limits
#' @param object \link{plotEset} object
#' @return matrix with limits for axes: columns \code{x} and \code{y}
#' @author Laure Cougnaud
#' @docType methods
#' @rdname getAxesLimits-methods
setGeneric(
	name = "getAxesLimits", 
	def = function(object)	standardGeneric("getAxesLimits")
)

#' @aliases getAxesLimits,esetPlot-method
#' @rdname getAxesLimits-methods
setMethod("getAxesLimits",
	signature = "esetPlot",
	definition = function(object){
			
	# get axes limits, depending on the 'symmetryAxes' parameter
	
	# if no genes plotted, only samples, otherwise use genes
	dataXY <- if((object@cloudGenes | object@topGenes > 0) & nrow(object@dataPlotGenes) > 0)	
		rbind(object@dataPlotSamples[, c('X', 'Y')], object@dataPlotGenes[, c('X', 'Y')])	else	object@dataPlotSamples[, c('X', 'Y')]
	
	# get maximum by axis
	maxCoordByAxis <- apply(abs(dataXY), 2, max)
	
	getAxisLimit <- function(x) c(-x, x)
	
	# define axes limits
	axesLimits <- switch(object@symmetryAxes, 
		'separate' = {
			res <- sapply(maxCoordByAxis, getAxisLimit); 
			colnames(res) <- c("x", "y"); 
			res
		},
		'combine' = sapply(c("x", "y"), function(x) 
			getAxisLimit(max(maxCoordByAxis))),
		
		'none' = {
			res <- apply(dataXY, 2, range);
			colnames(res) <- c("x", "y");
			res
		}
	)
	
	return(axesLimits)
	
})

#' get sample data for plot
#' @param object \link{plotEset} object
#' @return data.frame with 'dataPlotSamples' binded
#' with variables displayed in the plot
#' @author Laure Cougnaud
#' @docType methods
#' @rdname getDataPlotSamplesWithAnnotation-methods
setGeneric(
	name = "getDataPlotSamplesWithAnnotation", 
	def = function(object)	
		standardGeneric("getDataPlotSamplesWithAnnotation")
)

#' @aliases getDataPlotSamplesWithAnnotation,esetPlot-method
#' @rdname getDataPlotSamplesWithAnnotation-methods
setMethod("getDataPlotSamplesWithAnnotation",
	signature = "esetPlot",
	definition = function(object){
			
	esetMethods <- getMethodsInputObjectEsetVis(object@eset)		
		
	sampleAnnotation <- esetMethods$pData(object@eset)[, 
		c(object@colorVar, object@shapeVar, object@sizeVar, object@alphaVar), 
		drop = FALSE
	]

	dataPlotWithAnnotation <- cbind(object@dataPlotSamples, sampleAnnotation)
	
	colnames(dataPlotWithAnnotation) <- c(
		colnames(object@dataPlotSamples), 
		c(object@colorVar, object@shapeVar, 
			object@sizeVar, object@alphaVar)
	)
	
	return(dataPlotWithAnnotation)
	
})

#' @aliases getDataPlotSamplesWithAnnotation,ggvisEsetPlot-method
#' @rdname getDataPlotSamplesWithAnnotation-methods
#' @importFrom methods callNextMethod
setMethod("getDataPlotSamplesWithAnnotation",
	signature = "ggvisEsetPlot",
	definition = function(object){
		
	sampleAnnotation <- callNextMethod(object)	
	
	esetMethods <- getMethodsInputObjectEsetVis(object@eset)
	
	if(object@includeTooltip)
		sampleAnnotation <- cbind(sampleAnnotation, 
			keyggvis = esetMethods$sampleNames(object@eset))

	return(sampleAnnotation)
			
})

