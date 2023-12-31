#' @title plot a spectral map biplot of an \linkS4class{eSet}.
#' @description \code{esetSpectralMap} reduces the dimension of the data contained in the \linkS4class{eSet} with the \code{\link[mpm]{mpm}} function
#'  and plot the subsequent biplot of the specified dimensions, possibly with gene and sample annotation contained in the \linkS4class{eSet}.
#' A spectral map with the default parameters is equivalent to a 
#' principal component analysis on the log-transformed, double centered and 
#' global normalized data (from documentation of the \code{\link[mpm]{mpm}} function).
#' @param eset expressionSet (or SummarizedExperiment) object with data
#' @param psids featureNames of genes to include in the plot, all by default
#' @param dim dimensions of the analysis to represent, first two dimensions by default
#' @param mpm.args list with input parameters for the \code{\link[mpm]{mpm}} function.
#' The default value is:
#' \code{list(closure = 'none', center = 'double', normal = 'global', 'row.weight' = 'mean',
#' col.weight = 'constant', logtrans = FALSE)}.
#' This assumes that the data are already in a log scale.
#' @param plot.mpm.args list with input parameters for the \code{\link[mpm]{plot.mpm}} function.
#' The default value is: list(scale = "uvc").
#' @param returnAnalysis logical, if TRUE (FALSE by default), return also the output of the analysis,
#' and the outlying samples in the topElements element if any, otherwise only the plot object
#' @inheritParams esetPlotWrapper
#' @references Lewi, P.J. (1976). Spectral mapping, a technique for
#' classifying biological activity profiles of chemical compounds.
#' Arzneimittel Forschung (Drug Research), 26, 1295--1300
#' @seealso the function used internally: \code{\link[mpm]{mpm}} and 
#' \code{\link[a4Base]{spectralMap}} for spectral map in base R graphics
#' @examples
#' library(ALL)
#' data(ALL)
#' 
#' ## complete example (most of the parameters are optional)
#' # create custom color palette
#' colorPalette <- c("dodgerblue", colorRampPalette(c("white","dodgerblue2", "darkblue"))(5)[-1], 
#'		"red", colorRampPalette(c("white", "red3", "darkred"))(5)[-1])
#' # plot the spectral map
#' print(esetSpectralMap(eset = ALL, 
#' 	title = "Acute lymphoblastic leukemia dataset \n Spectral map complete",
#'	colorVar = "BT", color = colorPalette,
#'	shapeVar = "sex", shape = 15:16,
#'	sizeVar = "age", sizeRange = c(2, 6),
#'	symmetryAxes = "separate",
#'	topGenes = 10, topGenesJust = c(1, 0), topGenesCex = 2, topGenesColor = "darkgrey",
#'	topSamples = 15, topSamplesVar = "cod", topSamplesColor = "black",
#'	topSamplesJust = c(1, 0), topSamplesCex = 3)
#')
#' 
#' # see vignette for other examples, especially one with gene sets specification
#' 
#' @return if \code{returnAnalysis} is TRUE, return a list:
#' \itemize{
#'  \item{analysis: }{output of the spectral map analysis, can be given as input to the \code{esetPlotWrapper} function}
#'    \itemize{
#' 		\item{dataPlotSamples: }{coordinates of the samples}
#' 		\item{dataPlotGenes: }{coordinates of the genes}
#' 		\item{esetUsed: }{expressionSet used in the plot}
#'      \item{axisLabels: }{axes labels indicating percentage of variance explained by the selected axes}
#'      \item{axesContributionsPercentages: }percentages of variance explained by each axis (not only the ones specified in \code{dim})
#' 	  }
#'   \item{topElements: }{list with top outlying elements if any, possibly genes, samples and gene sets}
#'   \item{plot: }{the plot output}
#' }
#' otherwise return only the plot
#' @author Laure Cougnaud
#' @import Biobase mpm
#' @export
esetSpectralMap <- function(eset, 
	psids = 1:nrow(eset),
	dim = c(1, 2),
	colorVar = character(),
	color = if(length(colorVar) == 0)	"black"	else	character(),
	shapeVar = character(), 
	shape = if(length(shapeVar) == 0)	15	else	numeric(),
	sizeVar = character(),
	size = if(length(sizeVar) == 0){
		ifelse(typePlot[1] == "interactive" && packageInteractivity[1] == "plotly",
			20, 2.5
		)
	}else{numeric()},
	sizeRange = numeric(), #c(1, 6),
	alphaVar = character(), 
	alpha = if(length(alphaVar) == 0)	1	else	numeric(),  
	alphaRange = numeric(),
	title = "", 
	# assume that data are already log-transformed
	mpm.args = list(closure = "none", center = "double", 
		normal = "global", row.weight = "mean", col.weight = "constant",
		logtrans = FALSE),
	plot.mpm.args = list(scale = "uvc"),#except (x, dim, do.plot)
	symmetryAxes = c("combine", "separate", "none"),
	packageTextLabel = c("ggrepel", "ggplot2"),
	cloudGenes = TRUE, cloudGenesColor = "black",
	cloudGenesNBins = sqrt(length(psids)),
	cloudGenesIncludeLegend = FALSE, cloudGenesTitleLegend = "nGenes",
	topGenes = 10,
	topGenesCex = ifelse(
		typePlot[1] == "interactive" && packageInteractivity[1] == "plotly",
		10, 2.5),
	topGenesVar = character(), topGenesJust = c(0.5, 0.5), topGenesColor = "black",
	topSamples = 10,
	topSamplesCex = ifelse(
		typePlot[1] == "interactive" && packageInteractivity[1] == "plotly",
		10, 2.5),
	topSamplesVar = character(),
	topSamplesJust = c(0.5, 0.5), topSamplesColor = "black",
	geneSets = list(), geneSetsVar = character(), geneSetsMaxNChar = numeric(),
	topGeneSets = 10,
	topGeneSetsCex = ifelse(
		typePlot[1] == "interactive" && packageInteractivity[1] == "plotly",
		10, 2.5),
	topGeneSetsJust = c(0.5, 0.5),
	topGeneSetsColor = "black",
	includeLegend = TRUE, includeLineOrigin = TRUE,
	typePlot = c("static", "interactive"),
	packageInteractivity = c("plotly", "ggvis"),
	figInteractiveSize  = c(600, 400), ggvisAdjustLegend = TRUE,
	interactiveTooltip = TRUE, interactiveTooltipExtraVars = character(),
	returnAnalysis = FALSE, returnEsetPlot = FALSE){

	if(identical(packageInteractivity, "rbokeh")){
		.Deprecated(msg = paste("The 'rbokeh' interactive plot is deprecated",
			"(as the rbokeh package is archived), a 'plotly' interactive plot",
			"is created instead."
		))
		packageInteractivity <- "plotly"
	}

	symmetryAxes <- match.arg(symmetryAxes)
	packageTextLabel <- match.arg(packageTextLabel)
	packageInteractivity <- match.arg(packageInteractivity)
	
	if(length(psids) <= 1)
		stop(paste("At least two genes should be used for the visualization.",
			"Please change the 'psids' argument accordingly."))

	# get methods depending on the class of the object
	esetMethods <- getMethodsInputObjectEsetVis(eset)
	
	## Extract exprsMat with specified psids
	esetUsed <- eset[psids, ]

	## Create dataframe for plotting
	exprsMat <- esetMethods$exprs(esetUsed) #exprs
	dataMPM <- data.frame(rownames(exprsMat), exprsMat, stringsAsFactors = FALSE)
	colnames(dataMPM) <- c("featureID", colnames(exprsMat))
	dataMPMWthtNA <- dataMPM[rowSums(is.na(dataMPM)) == 0, ]
	
	## Run mpm with default parameters 
	argsMpm <- c(mpm.args, list("data" = dataMPMWthtNA))
	mpmResults <- do.call("mpm", argsMpm)
	
	## Convert spm results into plotting coordinates
	argsPlotMpm <- c(plot.mpm.args, 
		list("x" = mpmResults, "do.plot" = FALSE, "dim" = dim))
		#labels = fData(eset)$SYMBOL
	mpmCoords <- do.call("plot.mpm", argsPlotMpm)

	## Extract data for final plot
	dataPlotSamples <- data.frame(mpmCoords$Columns[, c("X", "Y")],
		sampleName = rownames(mpmCoords$Columns))

	dataPlotGenes <- data.frame(mpmCoords$Rows[, c("X", "Y")],
		geneName = rownames(mpmCoords$Rows))

	## label axes
	#extract percentages of variance in each axis
	pctVar <- round(mpmResults$contrib*100)
	axisLabels <- paste0("PC", dim, ": ", pctVar[dim], "%")

	## Plot
	plot <- esetVis::esetPlotWrapper(
		dataPlotSamples = dataPlotSamples,
		dataPlotGenes = dataPlotGenes, 
		esetUsed = esetUsed, 
		xlab = axisLabels[1], ylab = axisLabels[2],
		colorVar = colorVar, color = color,
		shapeVar = shapeVar, shape = shape,
		sizeVar = sizeVar, size = size, sizeRange = sizeRange,
		alphaVar = alphaVar, alpha = alpha, alphaRange = alphaRange,
		title = title, symmetryAxes = symmetryAxes,
		cloudGenes = cloudGenes, cloudGenesColor = cloudGenesColor, 
		cloudGenesNBins = cloudGenesNBins, 
		cloudGenesIncludeLegend = cloudGenesIncludeLegend, cloudGenesTitleLegend = cloudGenesTitleLegend,
		topGenes = topGenes, topGenesCex = topGenesCex, topGenesVar = topGenesVar, 
		topGenesJust = topGenesJust, topGenesColor = topGenesColor,
		topSamples = topSamples, topSamplesCex = topSamplesCex, 
		topSamplesVar = topSamplesVar, topSamplesJust = topSamplesJust, topSamplesColor = topSamplesColor,
		geneSets = geneSets, geneSetsVar = geneSetsVar, geneSetsMaxNChar = geneSetsMaxNChar, 
		topGeneSets = topGeneSets, topGeneSetsCex = topGeneSetsCex, topGeneSetsJust = topGeneSetsJust,
		topGeneSetsColor = topGeneSetsColor, 
		includeLegend = includeLegend, includeLineOrigin = includeLineOrigin,
		#nColLegend = nColLegend,
		typePlot = typePlot, 
		figInteractiveSize = figInteractiveSize, ggvisAdjustLegend = ggvisAdjustLegend, interactiveTooltip = interactiveTooltip, 
		interactiveTooltipExtraVars = interactiveTooltipExtraVars,
		returnTopElements = returnAnalysis,
		packageInteractivity = packageInteractivity,
		packageTextLabel = packageTextLabel,
		returnEsetPlot = returnEsetPlot)

	res <- if(!returnAnalysis)	plot else
		c(
			list(
				analysis = list(
					dataPlotSamples = dataPlotSamples,
					dataPlotGenes = dataPlotGenes, 
					esetUsed = esetUsed,
					axisLabels = axisLabels,
					axesContributionsPercentages = mpmResults$contrib*100
				)
			),
			if(!is.null(plot$topElements))	list(topElements = plot$topElements),
			list(plot = if(inherits(plot, "ggplot"))	plot	else	plot$plot)
		)

	return(res)
	
}
