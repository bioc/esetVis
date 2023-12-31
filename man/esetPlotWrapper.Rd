% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/esetPlotWrapper.R
\name{esetPlotWrapper}
\alias{esetPlotWrapper}
\title{wrapper for biplot of features/samples contained in a \linkS4class{eSet} object}
\usage{
esetPlotWrapper(
  dataPlotSamples,
  dataPlotGenes = data.frame(),
  esetUsed,
  xlab = "",
  ylab = "",
  colorVar = character(0),
  color = if (length(colorVar) == 0) "black" else character(0),
  shapeVar = character(0),
  shape = if (length(shapeVar) == 0) 15 else numeric(0),
  sizeVar = character(0),
  size = if (length(sizeVar) == 0) {
     ifelse(typePlot[1] == "interactive" &&
    packageInteractivity[1] == "plotly", 20, 2.5)
 } else {
     numeric()
 },
  sizeRange = numeric(0),
  alphaVar = character(0),
  alpha = if (length(alphaVar) == 0) 1 else numeric(0),
  alphaRange = numeric(0),
  title = "",
  symmetryAxes = c("combine", "separate", "none"),
  cloudGenes = TRUE,
  cloudGenesColor = "black",
  cloudGenesNBins = if (nrow(dataPlotGenes) > 0) sqrt(nrow(dataPlotGenes)) else numeric(),
  cloudGenesIncludeLegend = FALSE,
  cloudGenesTitleLegend = "nGenes",
  packageTextLabel = c("ggrepel", "ggplot2"),
  topGenes = 10,
  topGenesCex = ifelse(typePlot[1] == "interactive" && packageInteractivity[1] ==
    "plotly", 10, 2.5),
  topGenesVar = character(0),
  topGenesJust = c(0.5, 0.5),
  topGenesColor = "black",
  topSamples = 10,
  topSamplesCex = 2.5,
  topSamplesVar = character(0),
  topSamplesJust = c(0.5, 0.5),
  topSamplesColor = "black",
  geneSets = list(),
  geneSetsVar = character(0),
  geneSetsMaxNChar = numeric(0),
  topGeneSets = 10,
  topGeneSetsCex = 2.5,
  topGeneSetsJust = c(0.5, 0.5),
  topGeneSetsColor = "black",
  includeLegend = TRUE,
  includeLineOrigin = TRUE,
  typePlot = c("static", "interactive"),
  figInteractiveSize = c(600, 400),
  ggvisAdjustLegend = TRUE,
  interactiveTooltip = TRUE,
  interactiveTooltipExtraVars = character(0),
  packageInteractivity = c("plotly", "ggvis"),
  returnTopElements = FALSE,
  returnEsetPlot = FALSE
)
}
\arguments{
\item{dataPlotSamples}{data.frame with columns 'X', 'Y' with coordinates 
for the samples and with rownames which 
should correspond and be in the same order as the sampleNames of
\code{esetUsed}}

\item{dataPlotGenes}{data.frame with two columns 'X' and 'Y' with 
coordinates for the genes}

\item{esetUsed}{expressionSet (or SummarizedExperiment) object with data}

\item{xlab}{label for the x axis}

\item{ylab}{label for the y axis}

\item{colorVar}{name of variable (in varLabels of the \code{eset}) used 
for coloring, empty by default}

\item{color}{character or factor with specified color(s) for the points, 
replicated if needed. This is used only if \code{colorVar} is empty.
By default: 'black' if \code{colorVar} is not specified and default 
\code{ggplot} palette otherwise}

\item{shapeVar}{name of variable (in varLabels of the \code{eset}) 
used for the shape, empty by default}

\item{shape}{character or factor with specified shape(s) (pch) for the points, 
replicated if needed. This is used only if \code{shapeVar} is empty.
By default: '15' (filled square) if \code{shapeVar} is not specified 
and default \code{ggplot} shape(s) otherwise}

\item{sizeVar}{name of variable (in varLabels of the \code{eset}) 
used for the size, empty by default}

\item{size}{character or factor with specified size(s) (cex) for the points, 
replicated if needed.
This is used only if \code{sizeVar} is empty.
By default: '2.5' if \code{sizeVar} is not specified (20 for a plotly plot)
and default \code{ggplot} size(s) otherwise}

\item{sizeRange, }{size (cex) range used in the plot, possible only 
if the \code{sizeVar} is 'numeric' or 'integer'}

\item{alphaVar}{name of variable (in varLabels of the \code{eset}) 
used for the transparency, empty by default.
This parameter is currently only available for static plot and ggvis (only numeric in this case).}

\item{alpha}{character or factor with specified transparency(s) for the points,
replicated if needed. This is used only if \code{shapeVar} is empty. 
By default: '1' if \code{alphaVar} is not specified and default 
\code{ggplot} alpha otherwise
This parameter is currently only available for static and ggvis.}

\item{alphaRange}{transparency (alpha) range used in the plot, 
possible only if the \code{alphaVar} is 'numeric' or 'integer'
This parameter is currently only available for static and ggvis plot.}

\item{title}{plot title, '' by default}

\item{symmetryAxes}{set symmetry for axes, either:
\itemize{
 \item{'combine' (by default): }{both axes are symmetric and with the same limits}
 \item{'separate': }{each axis is symmetric and has its own limits}
 \item{'none': }{axes by default (plot limits)}
}}

\item{cloudGenes}{logical, if TRUE (by default), include the 
cloud of genes in the plot}

\item{cloudGenesColor}{if \code{cloudGenes} is TRUE, 
color for the cloud of genes, black by default}

\item{cloudGenesNBins}{number of bins to used for the clouds of genes,
by default the square root of the number of genes}

\item{cloudGenesIncludeLegend}{logical, if TRUE (FALSE by default) 
include the legend for the cloud of genes (in the top position if multiple legends)}

\item{cloudGenesTitleLegend}{string with title for the legend for the cloud of genes
'nGenes' by default}

\item{packageTextLabel}{package used to label the outlying genes/samples/gene sets,
either \code{ggrepel} (by default, only used if package \code{ggrepel} is available),
or \code{ggplot2}}

\item{topGenes}{numeric indicating which percentile (if <1) or number (if >=1) of genes
most distant to the origin of the plot to annotate, by default: 10 genes are selected
If no genes should be annotated, set this parameter to 0
Currently only available for static plot.}

\item{topGenesCex}{cex for gene annotation (used when \code{topGenes} > 0)}

\item{topGenesVar}{variable of the featureData used to label the genes, 
by default: empty, the featureNames are used for labelling (used when \code{topGenes} > 0)}

\item{topGenesJust}{text justification for the genes 
(used when \code{topGenes} > 0 and if \code{packageTextLabel} is \code{ggplot2}),
by default: c(0.5, 0.5) so centered}

\item{topGenesColor}{text color for the genes 
(used when \code{topGenes} > 0), black by default}

\item{topSamples}{numeric indicating which percentile (if <1) or number (if >=1) 
of samples most distant to the origin of the plot to annotate, 
by default: 10 samples are selected
If no samples should be annotated, set this parameter to 0.
Currently available for static plot.}

\item{topSamplesCex}{cex for sample annotation (used when \code{topSamples} > 0)}

\item{topSamplesVar}{variable of the phenoData used to label the samples, 
by default: empty, the sampleNames are used for labelling 
(used when \code{topSample}s > 0)}

\item{topSamplesJust}{text justification for the samples 
(used when \code{topSamples} > 0 and if \code{packageTextLabel} is \code{ggplot2}),
by default: c(0.5, 0.5) so centered}

\item{topSamplesColor}{text color for the samples 
(used when \code{topSamples} > 0), black by default}

\item{geneSets}{list of gene sets/pathways, each containing 
identifiers of genes contained in the set.
E.g. pathways from Gene Ontology databases output from the 
\code{\link{getGeneSetsForPlot}} function or any custom list of pathways.
The genes identifiers should correspond to the variable 
\code{geneSetsVar} contained in the phenoData, if not specified
the featureNames are used.
If several gene sets have the same name, they will be 
combine to extract the top gene sets.}

\item{geneSetsVar}{variable of the featureData used to 
match the genes contained in geneSets,
most probably ENTREZID, if not specified the featureNames
 of the eSet are used
Only used when \code{topGeneSets} > 0 and the parameter 
geneSets is specified.}

\item{geneSetsMaxNChar}{maximum number of characters for pathway names, 
by default keep entire names
Only used when \code{topGeneSets} > 0 and the parameter \code{geneSets} is specified.
If \code{returnAnalysis} is set to TRUE and \code{geneSetsMaxNChar} specified, 
the top pathways will be returned in the output object, named with the identifiers used in the plot 
(so with maximum \code{geneSetsMaxNChar} number of characters)}

\item{topGeneSets}{numeric indicating which percentile (if <=1) or number (if >1) of gene sets
most distant to the origin of the plot to annotate, by default: 10 gene sets are selected
If no gene sets should be annotated, set this parameter to 0.
Currently available for static plot.
Only used when \code{topGeneSets} > 0 and the parameter geneSets is specified.}

\item{topGeneSetsCex}{cex for gene sets annotation
Only used when \code{topGeneSets} > 0 and the parameter geneSets is specified.}

\item{topGeneSetsJust}{text justification for the gene sets 
by default: c(0.5, 0.5) so centered
Only used when \code{topGeneSets} > 0, the parameter \code{geneSets} 
is specified and if \code{packageTextLabel} is \code{ggplot2}.}

\item{topGeneSetsColor}{color for the gene sets 
(used when \code{topGeneSets} > 0 and \code{geneSets} is specified), black by default
Only used when \code{topGeneSets} > 0 and the parameter geneSets is specified.}

\item{includeLegend}{logical if TRUE (by default) include a legend, 
otherwise not}

\item{includeLineOrigin}{if TRUE (by default) include vertical line at 
x = 0 and horizontal line at y = 0}

\item{typePlot}{type of the plot returned, either 'static' (static) or 
interactive' (potentially interactive)}

\item{figInteractiveSize}{vector containing the size of the interactive plot, 
as [width, height]
by default: c(600, 400). This is passed to the \code{width} and 
\code{height} parameters of:
\itemize{
 \item{for plotly plots: }{the \code{\link[plotly]{ggplotly}} function}
 \item{for ggvis plots: }{the \code{ggvis::set_options} function}
}}

\item{ggvisAdjustLegend}{logical, if TRUE (by default) adjust the legends in \code{ggvis} to avoid
overlapping legends when multiple legends}

\item{interactiveTooltip}{logical, if TRUE, add hoover functionality showing
sample annotation (variables used in the plot) in the plot}

\item{interactiveTooltipExtraVars}{name of extra variable(s) 
(in varLabels of the \code{eset}) to add in plotlyEsetPlot to label the samples,
empty by default}

\item{packageInteractivity}{if \code{typePlot} is 'interactive', 
package used for interactive plot,
either 'plotly' (by default) (by default) or 'ggvis'.}

\item{returnTopElements}{logical, if TRUE return also the top elements}

\item{returnEsetPlot}{logical, if TRUE return also the \link{esetPlot} object}
}
\value{
if \code{typePlot} is:
\itemize{
	 \item{\code{static}: }
   \itemize{
    \item{if \code{returnTopElements} is TRUE, and top elements can be displayed, a list with:}
     \itemize{
      \item{'topElements': }{the top elements labelled in the plot}
	     \item{'plot': }{the \code{ggplot} object}
      }
    \item{otherwise, the \code{ggplot} object only}
   }
  \item{\code{interactive}: }{a \code{ggvis} or \code{plotly} object,
		depending on the \code{packageInteractivity} parameter}
}
}
\description{
Wrapper function used for all plots of the visualizations contained in the package.
}
\examples{
library(ALL)
data(ALL)

## run one spectral map analysis

# create custom color palette
colorPalette <- c("dodgerblue", colorRampPalette(c("white","dodgerblue2", "darkblue"))(5)[-1], 
	"red", colorRampPalette(c("white", "red3", "darkred"))(5)[-1])

# run the analysis
# with 'returnAnalysis' set to TRUE to have all objects required for the esetPlotWrapper
outputEsetSPM <- esetSpectralMap(eset = ALL, 
	title = "Acute lymphoblastic leukemia dataset \n Spectral map complete",
	colorVar = "BT", color = colorPalette,
	shapeVar = "sex", shape = 15:16,
	sizeVar = "age", sizeRange = c(2, 6),
	symmetryAxes = "separate",
	topGenes = 10, topGenesJust = c(1, 0), topGenesCex = 2, topGenesColor = "darkgrey",
	topSamples = 15, topSamplesVar = "cod", topSamplesColor = "black",
	topSamplesJust = c(1, 0), topSamplesCex = 3, returnAnalysis = TRUE)

# plot the biplot
print(outputEsetSPM$plot)


## re-call the plot function, to change some visualizations parameters
esetPlotWrapper(
	dataPlotSamples = outputEsetSPM$analysis$dataPlotSamples,
	dataPlotGenes = outputEsetSPM$analysis$dataPlotGenes,
	esetUsed = outputEsetSPM$analysis$esetUsed,
	title = paste("Acute lymphoblastic leukemia dataset \n Spectral map"),
	colorVar = "BT", color = colorPalette,
	shapeVar = "relapse", 
	sizeVar = "age", sizeRange = c(2, 6),
	topSamplesVar = "cod", topGenesVar = "SYMBOL"
)
}
\author{
Laure Cougnaud
}
