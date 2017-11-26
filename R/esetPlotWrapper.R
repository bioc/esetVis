#' wrapper for biplot of features/samples contained in a \linkS4class{eSet} object
#' 
#' Wrapper function used for all plots of the visualizations contained in the package.
#' @param dataPlotSamples data.frame with columns 'X', 'Y' with coordinates 
#' for the samples and with rownames which 
#' should correspond and be in the same order as the sampleNames of
#' \code{esetUsed}
#' @param dataPlotGenes data.frame with two columns 'X' and 'Y' with 
#' coordinates for the genes
#' @param esetUsed expressionSet (or SummarizedExperiment) object with data
#' @param xlab label for the x axis
#' @param ylab label for the y axis
#' @param colorVar name of variable (in varLabels of the \code{eset}) used 
#' for coloring, empty by default
#' @param color specified color(s) for the points, replicated if needed, 
#' used only if \code{colorVar} is empty, a factor or character
#' by default: 'black' if \code{colorVar} is not specified and default 
#' \code{ggplot} palette otherwise
#' @param shapeVar name of variable (in varLabels of the \code{eset}) 
#' used for the shape, empty by default
#' @param shape specified shape(s) (pch) for the points, replicated if 
#' needed, used only if \code{shapeVar} is empty, a factor or character
#' by default: '15' (filled square) if \code{shapeVar} is not specified 
#' and default \code{ggplot} shape(s) otherwise
#' @param sizeVar name of variable (in varLabels of the \code{eset}) 
#' used for the size, empty by default
#' @param size specified size(s) (cex) for the points, replicated if 
#' needed, used only if \code{sizeVar} is empty, a factor or character
#' by default: '2.5' if \code{sizeVar} is not specified and default 
#' \code{ggplot} size(s) otherwise
#' @param sizeRange, size (cex) range used in the plot, possible only 
#' if the \code{sizeVar} is 'numeric' or 'integer'
#' @param alphaVar name of variable (in varLabels of the \code{eset}) 
#' used for the transparency, empty by default.
#' This parameter is currently only available for static plot.
#' @param alpha specified transparency(s) for the points, replicated if 
#' needed, used only if \code{shapeVar} is empty, a factor or character
#' by default: '1' if \code{alphaVar} is not specified and default 
#' \code{ggplot} alpha otherwise
#' This parameter is currently only available for static plot.
#' @param alphaRange transparency (alpha) range used in the plot, 
#' possible only if the \code{alphaVar} is 'numeric' or 'integer'
#' This parameter is currently only available for static plot.
#' @param title plot title, '' by default
#' @param symmetryAxes set symmetry for axes, either:
#' \itemize{
#'  \item{'combine' (by default): }{both axes are symmetric and with the same limits}
#'  \item{'separate': }{each axis is symmetric and has its own limits}
#'  \item{'none': }{axes by default (plot limits)}
#' }
#' @param cloudGenes logical, if TRUE (by default), include the 
#' cloud of genes in the spectral map
#' @param cloudGenesColor if \code{cloudGenes} is TRUE, 
#' color for the cloud of genes, black by default
#' @param cloudGenesNBins number of bins to used for the clouds of genes,
#' by default the square root of the number of genes
#' @param cloudGenesIncludeLegend logical, if TRUE (FALSE by default) 
#' include the legend for the cloud of genes (in the top position if multiple legends)
#' @param cloudGenesTitleLegend string with title for the legend for the cloud of genes
#' 'nGenes' by default
#' @param packageTextLabel package used to label the outlying genes/samples/gene sets,
#' either \code{ggrepel} (by default, only used if package \code{ggrepel} is available),
#' or \code{ggplot2}
#' @param topGenes numeric indicating which percentile (if <1) or number (if >=1) of genes
#' most distant to the origin of the plot to annotate, by default: 10 genes are selected
#' If no genes should be annotated, set this parameter to 0
#' Currently only available for static plot.
#' @param topGenesCex cex for gene annotation (used when \code{topGenes} > 0)
#' @param topGenesVar variable of the featureData used to label the genes, 
#' by default: empty, the featureNames are used for labelling (used when \code{topGenes} > 0)
#' @param topGenesJust text justification for the genes 
#' (used when \code{topGenes} > 0 and if \code{packageTextLabel} is \code{ggplot2}),
#' by default: c(0.5, 0.5) so centered
#' @param topGenesColor text color for the genes 
#' (used when \code{topGenes} > 0), black by default
#' @param topSamples numeric indicating which percentile (if <1) or number (if >=1) 
#' of samples most distant to the origin of the plot to annotate, 
#' by default: 10 samples are selected
#' If no samples should be annotated, set this parameter to 0.
#' Currently available for static plot.
#' @param topSamplesCex cex for sample annotation (used when \code{topSamples} > 0)
#' @param topSamplesVar variable of the phenoData used to label the samples, 
#' by default: empty, the sampleNames are used for labelling 
#' (used when \code{topSample}s > 0)
#' @param topSamplesJust text justification for the samples 
#' (used when \code{topSamples} > 0 and if \code{packageTextLabel} is \code{ggplot2}),
#' by default: c(0.5, 0.5) so centered
#' @param topSamplesColor text color for the samples 
#' (used when \code{topSamples} > 0), black by default
#' @param geneSets list of gene sets/pathways, each containing 
#' identifiers of genes contained in the set.
#' E.g. pathways from Gene Ontology databases output from the 
#' \code{\link{getGeneSetsForPlot}} function or any custom list of pathways.
#' The genes identifiers should correspond to the variable 
#' \code{geneSetsVar} contained in the phenoData, if not specified
#' the featureNames are used.
#' If several gene sets have the same name, they will be 
#' combine to extract the top gene sets.
#' @param geneSetsVar variable of the featureData used to 
#' match the genes contained in geneSets,
#' most probably ENTREZID, if not specified the featureNames
#'  of the eSet are used
#' Only used when \code{topGeneSets} > 0 and the parameter 
#' geneSets is specified.
#' @param geneSetsMaxNChar maximum number of characters for pathway names, 
#' by default keep entire names
#' Only used when \code{topGeneSets} > 0 and the parameter \code{geneSets} is specified.
#' If \code{returnAnalysis} is set to TRUE and \code{geneSetsMaxNChar} specified, 
#' the top pathways will be returned in the output object, named with the identifiers used in the plot 
#' (so with maximum \code{geneSetsMaxNChar} number of characters)
#' @param topGeneSets numeric indicating which percentile (if <=1) or number (if >1) of gene sets
#' most distant to the origin of the plot to annotate, by default: 10 gene sets are selected
#' If no gene sets should be annotated, set this parameter to 0.
#' Currently available for static plot.
#' Only used when \code{topGeneSets} > 0 and the parameter geneSets is specified.
#' @param topGeneSetsCex cex for gene sets annotation
#' Only used when \code{topGeneSets} > 0 and the parameter geneSets is specified.
#' @param topGeneSetsJust text justification for the gene sets 
#' by default: c(0.5, 0.5) so centered
#' Only used when \code{topGeneSets} > 0, the parameter \code{geneSets} 
#' is specified and if \code{packageTextLabel} is \code{ggplot2}.
#' @param topGeneSetsColor color for the gene sets 
#' (used when \code{topGeneSets} > 0 and \code{geneSets} is specified), black by default
#' Only used when \code{topGeneSets} > 0 and the parameter geneSets is specified.
#' @param includeLegend logical if TRUE (by default) include a legend, 
#' otherwise not
#' @param includeLineOrigin if TRUE (by default) include vertical line at 
#' x = 0 and horizontal line at y = 0
#' @param typePlot type of the plot returned, either 'static' (static) or 
#' interactive' (potentially interactive)
#' @param figInteractiveSize vector containing the size of the interactive plot, 
#' as [width, height]
#' by default: c(600, 400). This is passed to the \code{width} and 
#' \code{height} parameters of:
#' \itemize{
#'  \item{for rbokeh plots: }{the \code{bokeh::figure} function}
#'  \item{for ggvis plots: }{the \code{ggvis::set_options} function}
#' }
#' @param ggvisAdjustLegend logical, if TRUE (by default) adjust the legends in \code{ggvis} to avoid
#' overlapping legends when multiple legends
#' @param interactiveTooltip logical, if TRUE, add hoover functionality showing
#' sample annotation (variables used in the plot) in the plot
#' @param interactiveTooltipExtraVars name of extra variable(s) 
#' (in varLabels of the \code{eset}) to add in rbokehEsetPlot to label the samples,
#' empty by default
#' @param packageInteractivity if \code{typePlot} is 'interactive', 
#' package used for interactive plot,
#' either 'rbokeh' (by default) or 'ggvis'
#' @param returnTopElements logical, if TRUE return also the top elements
#' @param returnEsetPlot logical, if TRUE return also the \link{esetPlot} object
#' @examples
#' library(ALL)
#' data(ALL)
#' 
#' ## run one spectral map analysis
#' 
#' # create custom color palette
#' colorPalette <- c("dodgerblue", colorRampPalette(c("white","dodgerblue2", "darkblue"))(5)[-1], 
#' 	"red", colorRampPalette(c("white", "red3", "darkred"))(5)[-1])
#' 
#' # run the analysis
#' # with 'returnAnalysis' set to TRUE to have all objects required for the esetPlotWrapper
#' outputEsetSPM <- esetSpectralMap(eset = ALL, 
#' 	title = "Acute lymphoblastic leukemia dataset \n Spectral map complete",
#' 	colorVar = "BT", color = colorPalette,
#' 	shapeVar = "sex", shape = 15:16,
#' 	sizeVar = "age", sizeRange = c(2, 6),
#' 	symmetryAxes = "separate",
#' 	topGenes = 10, topGenesJust = c(1, 0), topGenesCex = 2, topGenesColor = "darkgrey",
#' 	topSamples = 15, topSamplesVar = "cod", topSamplesColor = "black",
#' 	topSamplesJust = c(1, 0), topSamplesCex = 3, returnAnalysis = TRUE)
#' 
#' # plot the biplot
#' print(outputEsetSPM$plot)
#' 
#' 
#' ## re-call the plot function, to change some visualizations parameters
#' esetPlotWrapper(
#' 	dataPlotSamples = outputEsetSPM$analysis$dataPlotSamples,
#' 	dataPlotGenes = outputEsetSPM$analysis$dataPlotGenes,
#' 	esetUsed = outputEsetSPM$analysis$esetUsed,
#' 	title = paste("Acute lymphoblastic leukemia dataset \n Spectral map"),
#' 	colorVar = "BT", color = colorPalette,
#' 	shapeVar = "relapse", 
#' 	sizeVar = "age", sizeRange = c(2, 6),
#' 	topSamplesVar = "cod", topGenesVar = "SYMBOL"
#' )
#' @return if \code{typePlot} is:
#' \itemize{
#' 	 \item{\code{static}: }
#'    \itemize{
#'     \item{if \code{returnTopElements} is TRUE, and top elements can be displayed, a list with:}
#'      \itemize{
#'       \item{'topElements': }{the top elements labelled in the plot}
#' 	     \item{'plot': }{the \code{ggplot} object}
#'       }
#'     \item{otherwise, the \code{ggplot} object only}
#'    }
#'   \item{\code{interactive}: }{a \code{ggvis} or \code{rbokeh} object, depending on the \code{packageInteractivity} parameter}
#' }
#' @import Biobase
#' @importFrom hexbin hexbin
#' @importFrom grDevices colorRamp colorRampPalette
#' @importFrom utils getFromNamespace
#' @author Laure Cougnaud
#' @export
esetPlotWrapper <- function(
	dataPlotSamples, 
	dataPlotGenes = data.frame(), esetUsed, 
	xlab = "", ylab = "",
	colorVar = character(0), 
	color = if(length(colorVar) > 0)	"black"	else	character(0),
	shapeVar = character(0), 
	shape = if(length(shapeVar) > 0)	15	else	numeric(0),
	sizeVar = character(0), 
	size = if(length(sizeVar) > 0)	2.5	else	numeric(0), 
	sizeRange = numeric(0), #c(1, 6),
	alphaVar = character(0), 
	alpha = if(length(alphaVar) > 0)	1	else	numeric(0), 
	alphaRange = numeric(0),
	title = "", 
	#assume that data are already log-transformed
	symmetryAxes = c("combine", "separate", "none"),
	cloudGenes = TRUE, cloudGenesColor = "black",
	cloudGenesNBins = if(nrow(dataPlotGenes) > 0)	sqrt(nrow(dataPlotGenes))	else	numeric(),
	cloudGenesIncludeLegend = FALSE, cloudGenesTitleLegend = "nGenes",
	packageTextLabel = c("ggrepel", "ggplot2"),
	topGenes = 10, topGenesCex = 2.5, topGenesVar = character(0), topGenesJust = c(0.5, 0.5), topGenesColor = "black",
	topSamples = 10, topSamplesCex = 2.5, topSamplesVar = character(0), topSamplesJust = c(0.5, 0.5), topSamplesColor = "black",
	geneSets = list(), geneSetsVar = character(0), geneSetsMaxNChar = numeric(0), topGeneSets = 10, 
	topGeneSetsCex = 2.5, topGeneSetsJust = c(0.5, 0.5), topGeneSetsColor = "black",
	includeLegend = TRUE, includeLineOrigin = TRUE,
	typePlot = c("static", "interactive"),
	figInteractiveSize  = c(600, 400), ggvisAdjustLegend = TRUE,
	interactiveTooltip = TRUE, interactiveTooltipExtraVars = character(0),
	packageInteractivity = c("rbokeh", "ggvis"),
	returnTopElements = FALSE, returnEsetPlot = FALSE){

	symmetryAxes <- match.arg(symmetryAxes)
	packageTextLabel <- match.arg(packageTextLabel)
	packageInteractivity <- match.arg(packageInteractivity)
	typePlot <- match.arg(typePlot)
	
	# get methods depending on the class of the object
#	esetMethods <- getMethodsInputObjectEsetVis(esetUsed)
	
	# add needed variables in the input data for the plot
	## Bug #7160 fix - added drop = FALSE to handle cases with only one 'Var'
		
	## functions common between the different type of plots available

	esetPlotArgs <- list(
		eset = esetUsed,
		dataPlotSamples = dataPlotSamples,
		dataPlotGenes = dataPlotGenes,
		cloudGenes = cloudGenes, 
		cloudGenesColor = cloudGenesColor, 
		cloudGenesNBins = cloudGenesNBins, 
		cloudGenesIncludeLegend = cloudGenesIncludeLegend,
		cloudGenesTitleLegend = cloudGenesTitleLegend,
		colorVar = colorVar, 
		color = color, 
		shapeVar = shapeVar, 
		shape = shape, 
		sizeVar = sizeVar , 
		size = size, 
		alphaVar = alphaVar,
		alpha = alpha,
		sizeRange = sizeRange, 
		alphaRange = alphaRange,
		includeLineOrigin = includeLineOrigin,
		title = title,
		xlab = xlab, ylab = ylab, 
		symmetryAxes = symmetryAxes,
		topGenes = topGenes, 
		topGenesVar = topGenesVar, 
		topGenesCex = topGenesCex, 
		topGenesJust = topGenesJust, 
		topGenesColor = topGenesColor, 
		topSamples = topSamples, 
		topSamplesCex = topSamplesCex, 
		topSamplesVar = topSamplesVar, 
		topSamplesJust = topSamplesJust, 
		topSamplesColor = topSamplesColor,
		geneSets = geneSets, 
		geneSetsVar = geneSetsVar, 
		geneSetsMaxNChar = geneSetsMaxNChar, 
		topGeneSets = topGeneSets,
		topGeneSetsCex = topGeneSetsCex, 
		topGeneSetsJust = topGeneSetsJust, 
		topGeneSetsColor = topGeneSetsColor,
		includeLegend = includeLegend,
		packageTextLabel = packageTextLabel
	)
	
	esetPlotArgs <- c(esetPlotArgs,
			
		switch(typePlot,
	
			'static' = {
				list(returnTopElements = returnTopElements)
			},
				
			'interactive' = {
			
				c(
						
					list(
						includeTooltip = interactiveTooltip, 
						tooltipVars = interactiveTooltipExtraVars, 
						sizePlot = figInteractiveSize
					),
				
					if(packageInteractivity == "ggvis")
						list(adjustLegend = ggvisAdjustLegend)
					
				)
				
			}
		
		)
		
	)
	
	esetPlotFct <- switch(typePlot,
		"static" = "ggplotEsetPlot",
		"interactive" = switch(packageInteractivity,
			'rbokeh' = "rbokehEsetPlot",
			'ggvis' = "ggvisEsetPlot"
		)
	)
	inputPlot <- do.call(esetPlotFct, esetPlotArgs)
	
	res <- plotEset(inputPlot, returnEsetPlot = returnEsetPlot)

	return(res)
	
}
