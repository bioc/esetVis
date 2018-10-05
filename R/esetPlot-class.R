#' An S4 class to represent \code{esetPlot} object
#' expressionSet with visualization data from
#' dimension-reduction methods
#'
#' @slot dataPlotSamples data.frame with columns 'X', 'Y' with coordinates 
#' for the samples and with rownames which 
#' should correspond and be in the same order as the sampleNames of
#' \code{esetUsed}
#' @slot dataPlotGenes data.frame with two columns 'X' and 'Y' with 
#' coordinates for the genes
#' @slot eset expressionSet (or SummarizedExperiment) object with data
#' @slot colorVar name of variable (in varLabels of the \code{eset}) used 
#' for coloring, empty by default
#' @slot color character or factor with specified color(s) for the points, 
#' replicated if needed. This is used only if \code{colorVar} is empty.
#' By default: 'black' if \code{colorVar} is not specified and default 
#' \code{ggplot} palette otherwise
#' @slot shapeVar name of variable (in varLabels of the \code{eset}) 
#' used for the shape, empty by default
#' @slot shape character or factor with specified shape(s) (pch) for the points, 
#' replicated if needed. This is used only if \code{shapeVar} is empty.
#' By default: '15' (filled square) if \code{shapeVar} is not specified 
#' and default \code{ggplot} shape(s) otherwise
#' @slot sizeVar name of variable (in varLabels of the \code{eset}) 
#' used for the size, empty by default
#' @slot size size character or factor with specified size(s) (cex) for the points, 
#' replicated if needed.
#' This is used only if \code{sizeVar} is empty.
#' By default: '2.5' if \code{sizeVar} is not specified and default 
#' \code{ggplot} size(s) otherwise
#' @slot sizeRange, size (cex) range used in the plot, possible only 
#' if the \code{sizeVar} is 'numeric' or 'integer'
#' @slot alphaVar name of variable (in varLabels of the \code{eset}) 
#' used for the transparency, empty by default.
#' @slot alpha alpha character or factor with specified transparency(s) for the points,
#' replicated if needed. This is used only if \code{shapeVar} is empty. 
#' By default: '1' if \code{alphaVar} is not specified and default 
#' \code{ggplot} alpha otherwise.
#' @slot alphaRange transparency (alpha) range used in the plot, 
#' possible only if the \code{alphaVar} is 'numeric' or 'integer'
#' This parameter is not available for rbokeh plot.
#' @slot symmetryAxes set symmetry for axes, either:
#' \itemize{
#'  \item{'combine' (by default): }{both axes are symmetric and with the same limits}
#'  \item{'separate': }{each axis is symmetric and has its own limits}
#'  \item{'none': }{axes by default (plot limits)}
#' }
#' @slot cloudGenes logical, if TRUE (by default), include the 
#' cloud of genes in the spectral map
#' @slot cloudGenesColor if \code{cloudGenes} is TRUE, 
#' color for the cloud of genes, black by default
#' @slot cloudGenesNBins number of bins to used for the clouds of genes,
#' by default the square root of the number of genes
#' @slot cloudGenesIncludeLegend logical, if TRUE (FALSE by default) 
#' include the legend for the cloud of genes (in the top position if multiple legends)
#' @slot cloudGenesTitleLegend string with title for the legend for the cloud of genes
#' 'nGenes' by default
#' @slot packageTextLabel package used to label the outlying genes/samples/gene sets,
#' either \code{ggrepel} (by default, only used if package \code{ggrepel} is available),
#' or \code{ggplot2}
#' @slot topGenes numeric indicating which percentile (if <1) or number (if >=1) of genes
#' most distant to the origin of the plot to annotate, by default: 10 genes are selected
#' If no genes should be annotated, set this parameter to 0
#' Currently only available for static plot.
#' @slot topGenesCex cex for gene annotation (used when \code{topGenes} > 0)
#' @slot topGenesVar variable of the featureData used to label the genes, 
#' by default: empty, the featureNames are used for labelling (used when \code{topGenes} > 0)
#' @slot topGenesJust text justification for the genes 
#' (used when \code{topGenes} > 0 and if \code{packageTextLabel} is \code{ggplot2}),
#' by default: c(0.5, 0.5) so centered
#' @slot topGenesColor text color for the genes 
#' (used when \code{topGenes} > 0), black by default
#' @slot topSamples numeric indicating which percentile (if <1) or number (if >=1) 
#' of samples most distant to the origin of the plot to annotate, 
#' by default: 10 samples are selected
#' If no samples should be annotated, set this parameter to 0.
#' Currently available for static plot.
#' @slot topSamplesCex cex for sample annotation (used when \code{topSamples} > 0)
#' @slot topSamplesVar variable of the phenoData used to label the samples, 
#' by default: empty, the sampleNames are used for labelling 
#' (used when \code{topSample}s > 0)
#' @slot topSamplesJust text justification for the samples 
#' (used when \code{topSamples} > 0 and if \code{packageTextLabel} is \code{ggplot2}),
#' by default: c(0.5, 0.5) so centered
#' @slot topSamplesColor text color for the samples 
#' (used when \code{topSamples} > 0), black by default
#' @slot geneSets list of gene sets/pathways, each containing 
#' identifiers of genes contained in the set.
#' E.g. pathways from Gene Ontology databases output from the 
#' \code{\link{getGeneSetsForPlot}} function or any custom list of pathways.
#' The genes identifiers should correspond to the variable 
#' \code{geneSetsVar} contained in the phenoData, if not specified
#' the featureNames are used.
#' If several gene sets have the same name, they will be 
#' combine to extract the top gene sets.
#' @slot geneSetsVar variable of the featureData used to 
#' match the genes contained in geneSets,
#' most probably ENTREZID, if not specified the featureNames
#'  of the eSet are used
#' Only used when \code{topGeneSets} > 0 and the parameter 
#' geneSets is specified.
#' @slot geneSetsMaxNChar maximum number of characters for pathway names, 
#' by default keep entire names
#' Only used when \code{topGeneSets} > 0 and the parameter \code{geneSets} is specified.
#' @slot topGeneSets numeric indicating which percentile (if <=1) or number (if >1) of gene sets
#' most distant to the origin of the plot to annotate, by default: 10 gene sets are selected
#' If no gene sets should be annotated, set this parameter to 0.
#' Currently available for static plot.
#' Only used when \code{topGeneSets} > 0 and the parameter geneSets is specified.
#' @slot topGeneSetsCex cex for gene sets annotation
#' Only used when \code{topGeneSets} > 0 and the parameter geneSets is specified.
#' @slot topGeneSetsJust text justification for the gene sets 
#' by default: c(0.5, 0.5) so centered
#' Only used when \code{topGeneSets} > 0, the parameter \code{geneSets} 
#' is specified and if \code{packageTextLabel} is \code{ggplot2}.
#' @slot topGeneSetsColor color for the gene sets 
#' (used when \code{topGeneSets} > 0 and \code{geneSets} is specified), black by default
#' Only used when \code{topGeneSets} > 0 and the parameter geneSets is specified.
#' @slot includeLegend logical if TRUE (by default) include a legend, 
#' otherwise not
#' @slot includeLineOrigin if TRUE (by default) include vertical line at 
#' x = 0 and horizontal line at y = 0
#' @return S4 object of class \code{esetPlot}
#' @import Biobase
#' @name esetPlot-class
#' @rdname esetPlot-class
#' @importFrom methods new
#' @export
esetPlot <- setClass("esetPlot", 
	slots = c(
		eset = "ExpressionSet",
		dataPlotSamples = "data.frame",
		dataPlotGenes = "data.frame",
		cloudGenes = "logical", 
		cloudGenesColor = "character", 
		cloudGenesNBins = "numeric", 		
		cloudGenesIncludeLegend = "logical",
		cloudGenesTitleLegend = "character",
		colorVar = "character", 
		color = "character", 
		shapeVar = "character", 
		shape = "numeric", 
		sizeVar = "character", 
		size = "numeric",
		alphaVar = "character",
		alpha = "numeric",
		sizeRange = "numeric", 
		alphaRange = "numeric",
		includeLineOrigin = "logical",
		symmetryAxes = "character",
		topGenes = "numeric", 
		topGenesVar = "character", 
		topGenesCex = "numeric", 
		topGenesJust = "numeric", 
		topGenesColor = "character", 
		topSamples = "numeric", 
		topSamplesCex = "numeric", 
		topSamplesVar = "character", 
		topSamplesJust = "numeric", 
		topSamplesColor = "character",
		geneSets = "list", 
		geneSetsVar = "character", 
		geneSetsMaxNChar = "numeric", 
		topGeneSets = "numeric",
		topGeneSetsCex = "numeric", 
		topGeneSetsJust = "numeric", 
		topGeneSetsColor = "character",
		includeLegend = "logical",
		packageTextLabel = "character"
	),
	prototype = list(
		cloudGenes = TRUE, 
		cloudGenesColor = "black", 
		cloudGenesIncludeLegend = FALSE,
		cloudGenesTitleLegend = "character",
		color = "black", 
		shape = 15, 
		alpha = 1,
		size = 2.5,
		includeLineOrigin = TRUE,
		symmetryAxes = "combine",
		topGenes = 10, 
		topGenesCex = 2.5, 
		topGenesJust = c(0.5, 0.5), 
		topGenesColor = "black",
		topSamples = 10, 
		topSamplesCex = 2.5, 
		topSamplesJust = c(0.5, 0.5), 
		topSamplesColor = "black",
		topGeneSets = 10,
		topGeneSetsCex = 2.5, 
		topGeneSetsJust = c(0.5, 0.5), 
		topGeneSetsColor = "black",
		packageTextLabel = "ggrepel", 
		includeLegend = TRUE
	)
)

#' constructor for \link{esetPlot}
#' 
#' Constructor of the \link{esetPlot} class
#' @name esetPlot
#' @aliases initialize,esetPlot-method
#' @rdname esetPlot-class
#' @param .Object \link{esetPlot} object
#' @param ... additional class arguments
#' @importFrom methods callNextMethod
setMethod("initialize", "esetPlot", function(.Object, ...) {
			
	.Object <- callNextMethod(.Object, ...)
	if(length(.Object@cloudGenesNBins) == 0 && nrow(.Object@dataPlotGenes) > 0){
		.Object@cloudGenesNBins <- sqrt(nrow(.Object@dataPlotGenes))
	}
	return(.Object)
			
})



#' @title S4 Class Union with character/expression/call
#' @description This is used for the definition of the title/axes labels for the ggplot2 version
setClassUnion(name = "characterORexpressionOrCall", members = c("character", "expression", "call"))

#' a S4 class to represent \code{ggplot} plots
#' @slot returnTopElements logical, if TRUE (FALSE by default) 
#' return the outlying elements labelled in the plot (if any)
#' @slot title string or expression with plot title, '' by default
#' @slot xlab string or expression with label for the x axis
#' @slot ylab string or expression with label for the y axis
#' @return S4 object of class \code{ggplotEsetPlot}
#' @author Laure Cougnaud
#' @name ggplotEsetPlot-class
#' @rdname ggplotEsetPlot-class
#' @export
#' @importFrom methods new
ggplotEsetPlot <- setClass("ggplotEsetPlot", 
	slots = c(
		returnTopElements = "logical",
		title = "characterORexpressionOrCall",
		xlab = "characterORexpressionOrCall", 
		ylab = "characterORexpressionOrCall"
	),
	prototype = list(
		returnTopElements = FALSE,
		title = "",	xlab = "", ylab = ""
	),
	contains = "esetPlot"
)

#' a S4 class to represent interactive plots
#' @slot includeTooltip logical, if TRUE, add hoover functionality showing
#' sample annotation (variables used in the plot) in the plot
#' @slot tooltipVars name of extra phenotypic variable(s) 
#' to add in rbokehEsetPlot to label the samples
#' @slot sizePlot vector containing the size of the interactive plot, 
#' as [width, height], by default: c(600, 400).
#' @slot title string plot title, '' by default
#' @slot xlab string label for the x axis
#' @slot ylab string label for the y axis
#' @return S4 object of class \code{esetPlotInteractive}
#' @author Laure Cougnaud
#' @name esetPlotInteractive-class
#' @rdname esetPlotInteractive-class
#' @export
#' @importFrom methods new
esetPlotInteractive <- setClass("esetPlotInteractive", 
	slots = c(
		includeTooltip = "logical",
		tooltipVars = "character",
		sizePlot = "numeric",
		title = "character", xlab = "character", ylab = "character"
	),
	prototype = list(
		includeTooltip = TRUE,
		sizePlot = c(600, 400),
		title = "",	xlab = "", ylab = ""
	),
	contains = "esetPlot"
)

#' a S4 class to represent \code{rbokeh} plots
#' @slot size specified size(s) (cex) for the points, replicated if 
#' needed, used only if \code{sizeVar} is empty, a factor or character
#' by default: '5' if \code{sizeVar} is not specified and default 
#' \code{ggplot} size(s) otherwise
#' @return S4 object of class \code{rbokehEsetPlot}
#' @author Laure Cougnaud
#' @name rbokehEsetPlot-class
#' @rdname rbokehEsetPlot-class
#' @export
#' @importFrom methods new
rbokehEsetPlot <- setClass("rbokehEsetPlot", 
	slots = list(size = "numeric"),
	prototype = list(size = 5),
	contains = "esetPlotInteractive"
)

#' a S4 class for \code{ggvis} plot
#' @slot adjustLegend logical, if TRUE (by default) adjust the legends in \code{ggvis} to avoid
#' overlapping legends when multiple legends
#' @slot alphaVar name of numeric variable (in varLabels of the \code{eset}) 
#' used for the transparency, empty by default.
#' @slot alphaRange transparency (alpha) range used in the plot, 
#' c(0.1, 1) by default.
#' @return S4 object of class \code{ggvisEsetPlot}
#' @author Laure Cougnaud
#' @name ggvisEsetPlot-class
#' @rdname ggvisEsetPlot-class
#' @export
#' @importFrom methods new
ggvisEsetPlot <- setClass("ggvisEsetPlot", 
	slots = c("adjustLegend" = "logical"),
	prototype = list(
		adjustLegend = TRUE,
		alphaRange = c(0.1, 1)
	),
	contains = "esetPlotInteractive"
)