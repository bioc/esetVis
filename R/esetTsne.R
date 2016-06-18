#' @title plot a t-SNE of an \linkS4class{eSet} object
#' @description \code{esetTsne} reduces the dimension of the data contained in the \linkS4class{eSet} via t-Distributed Stochastic Neighbor Embedding
#'  with the \code{Rtsne} function and plot the subsequent biplot, possibly with sample annotation contained in the eSet.
#' @param eset expressionSet (or SummarizedExperiment) object with data
#' @param psids featureNames of genes to include in the plot, all by default
#' @param trace logical, if TRUE (by default), print some messages during tsne is running
#' @param Rtsne.args arguments for the Rtsne function, by default:
#' perplexite parameter = optimal number of neighbours, 
#' theta = speed/accuracy trade-off (increase for less accuracy), set to 0.0 for exact TSNE
#' @param fctTransformDataForInputTsne function which transform the data in the eSet object before
#' calling the \code{RTnse} function. 
#' This should be a function which takes a matrix as input and return a matrix, e.g. the dist function.
#' @param returnAnalysis logical, if TRUE (FALSE by default), return also the output of the analysis,
#' and the outlying samples in the topElements element if any, otherwise only the plot object
#' @inheritParams esetPlotWrapper
#' @references L.J.P. van der Maaten and G.E. Hinton (2008). Visualizing
#' High-Dimensional Data Using t-SNE. Journal of Machine Learning
#' Research, 2579--2605
#' @seealso the function used internally: \code{\link[Rtsne]{Rtsne}} or \url{http://homepage.tudelft.nl/19j49/t-SNE.html}
#'  for further explanations about this technique.
#' @return if \code{returnAnalysis} is TRUE, return a list:
#' \itemize{
#'  \item{analysis: }{output of the spectral map analysis, can be given as input to the \code{esetPlotWrapper} function}
#'    \itemize{
#' 		\item{dataPlotSamples: }{coordinates of the samples}
#' 		\item{esetUsed: }{expressionSet used in the plot}
#'      \item{axisLabels: }{axes labels indicating percentage of variance explained by the selected axes}
#' 	  }
#'   \item{topElements: }{list with top outlying elements if any, possibly genes, samples and gene sets}
#'   \item{plot: }{the plot output}
#' }
#' otherwise return only the plot
#' @examples 
#' library(ALL)
#' data(ALL)
#' 
#' ## complete example (most of the parameters are optional)
#' 
#' # create custom color palette
#' colorPalette <- c("dodgerblue", colorRampPalette(c("white","dodgerblue2", "darkblue"))(5)[-1], 
#'		"red", colorRampPalette(c("white", "red3", "darkred"))(5)[-1])
#' 
#' # create tsne
#' print(esetTsne(eset = ALL, 
#' 	title = "Acute lymphoblastic leukemia dataset \n Tsne complete",
#' 	colorVar = "BT", color = colorPalette,
#' 	shapeVar = "sex", shape = 15:16,
#' 	sizeVar = "age", sizeRange = c(2, 6),
#' 	symmetryAxes = "separate",
#' 	topSamples = 15, topSamplesVar = "cod", topSamplesColor = "black",
#' 	topSamplesJust = c(1, 0), topSamplesCex = 3)
#' )
#' @author Laure Cougnaud
#' @import Biobase
#' @importFrom utils capture.output
#' @importFrom Rtsne Rtsne
#' @export
esetTsne <- function(eset, 
	psids = 1:nrow(eset),
	trace = TRUE,
	colorVar = NULL, #color = "black", colorVarValues = NULL, 
	color = if(is.null(colorVar))	"black"	else	NULL,
	shapeVar = NULL, 
	shape = if(is.null(shapeVar))	15	else	NULL,
	sizeVar = NULL, 
	size = if(is.null(sizeVar))	2.5	else	NULL, 
	sizeRange = NULL, #c(1, 6),
	alphaVar = NULL, 
	alpha = if(is.null(alphaVar))	1	else	NULL, alphaRange = NULL,
	title = "", 
	#parameters for tsne
	Rtsne.args = list(perplexity = floor((ncol(eset)-1)/3),
		theta = 0.5, dims = 2, initial_dims = 50),
	fctTransformDataForInputTsne = NULL,
	symmetryAxes = c("combine", "separate", "none"),
	packageTextLabel = c("ggrepel", "ggplot2"),
	#cloudGenes = TRUE, cloudGenesColor = "black",
	#topGenes = 10, topGenesCex = 2.5, topGenesVar = NULL, topGenesJust = c(0.5, 0.5),
	topSamples = 10, topSamplesCex = 2.5, topSamplesVar = NULL, topSamplesJust = c(0.5, 0.5), topSamplesColor = "black",
	#geneSets = list(), geneSetsVar = NULL, geneSetsMaxNChar = NULL, topGeneSets = 10, topGeneSetsCex = 2.5, topGeneSetsJust = c(0.5, 0.5),
	includeLegend = TRUE, includeLineOrigin = TRUE,
	typePlot = c("static", "interactive"), packageInteractivity = c("rbokeh", "ggvis"),
	figInteractiveSize  = c(600, 400), ggvisAdjustLegend = TRUE, interactiveTooltip = TRUE, interactiveTooltipExtraVars = NULL,
	returnAnalysis = FALSE){

	symmetryAxes <- match.arg(symmetryAxes)
	packageInteractivity <- match.arg(packageInteractivity)
	packageTextLabel <- match.arg(packageTextLabel)
	
	# get methods depending on the class of the object
	esetMethods <- getMethodsInputObjectEsetVis(eset)

	# to have reproducible results
	set.seed(123)
	
	# criteria on perplexity in Rtsne:
	#nrow(dat) - 1 < 3 * perplexity
	
	symmetryAxes <- match.arg(symmetryAxes)
	
	if(length(psids) <= 1)
		stop(paste("At least two genes should be used for the visualization.",
			"Please change the 'psids' argument accordingly."))
	
	## Extract exprsMat with specified psids
	esetUsed <- eset[psids, ]
	
	## Prepare data for tsne
	inputTsne <- if(!is.null(fctTransformDataForInputTsne))	
		fctTransformDataForInputTsne(esetMethods$exprs(esetUsed))	else	t(esetMethods$exprs(esetUsed))
	#as.matrix((1 - cor(exprs(esetUsed)))/2)# or as.dist((1 - cor(dat))/2)
	argsRtsne <- c(Rtsne.args, list(X = inputTsne))
	
	## Run tsne
	wrapperTnse <- function()	do.call("Rtsne", argsRtsne)
	if(trace)	tsneOut <- wrapperTnse()	else
		outputMessages <- capture.output(tsneOut <- wrapperTnse())
	
	## Extract data for final plot
	dataPlotSamples <- data.frame(tsneOut$Y, esetMethods$sampleNames(esetUsed))
	colnames(dataPlotSamples) <- c("X", "Y", "sampleName")
	rownames(dataPlotSamples) <- as.character(dataPlotSamples$sampleName)
	
	plot <- esetVis::esetPlotWrapper(dataPlotSamples = dataPlotSamples,
		dataPlotGenes = NULL, 
		esetUsed = esetUsed, 
		colorVar = colorVar, color = color,
		shapeVar = shapeVar, shape = shape,
		sizeVar = sizeVar, size = size, sizeRange = sizeRange,
		alphaVar = alphaVar, alpha = alpha, alphaRange = alphaRange,
		title = title, symmetryAxes = symmetryAxes,
		topSamples = topSamples, topSamplesCex = topSamplesCex, topSamplesVar = topSamplesVar, topSamplesJust = topSamplesJust,
		includeLegend = includeLegend, includeLineOrigin = includeLineOrigin,
		typePlot = typePlot, 
		figInteractiveSize = figInteractiveSize, ggvisAdjustLegend = ggvisAdjustLegend, interactiveTooltip = interactiveTooltip, interactiveTooltipExtraVars = interactiveTooltipExtraVars,
		xlab = "X", ylab = "Y",
		returnTopElements = returnAnalysis,
		packageInteractivity = packageInteractivity,
		packageTextLabel = packageTextLabel)

	someTopElementsReturned <- any(!class(plot) == "ggplot")
	
	res <- if(!returnAnalysis)	plot else
		c(
			list(analysis = list(
				dataPlotSamples = dataPlotSamples,
				esetUsed = esetUsed
			)),
			if(someTopElementsReturned)	list(topElements = plot$topElements)	else NULL,
			list(plot = plot$plot)
		)

	return(res)

}
