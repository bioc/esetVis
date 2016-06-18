#' wrapper for all the plots of the package
#' @param dataPlotSamples data.frame with columns 'X', 'Y' with coordinates for the samples and with rownames which 
#' should correspond and be in the same order as the sampleNames of
#' \code{esetUsed}
#' @param dataPlotGenes data.frame with two columns 'X' and 'Y' with coordinates for the genes
#' @param esetUsed expressionSet (or SummarizedExperiment) object with data
#' @param xlab label for the x axis
#' @param ylab label for the y axis
#' @param colorVar name of variable (in varLabels of the \code{eset}) used for coloring, NULL by default
#' @param color specified color(s) for the points, replicated if needed, used only if \code{colorVar} is NULL, a factor or character
#' by default: 'black' if \code{colorVar} is not specified and default \code{ggplot} palette otherwise
#' @param shapeVar name of variable (in varLabels of the \code{eset}) used for the shape, NULL by default
#' @param shape specified shape(s) (pch) for the points, replicated if needed, used only if \code{shapeVar} is NULL, a factor or character
#' by default: '15' (filled square) if \code{shapeVar} is not specified and default \code{ggplot} shape(s) otherwise
#' @param sizeVar name of variable (in varLabels of the \code{eset}) used for the size, NULL by default
#' @param size specified size(s) (cex) for the points, replicated if needed, used only if \code{sizeVar} is NULL, a factor or character
#' by default: '2.5' if \code{sizeVar} is not specified and default \code{ggplot} size(s) otherwise
#' @param sizeRange, size (cex) range used in the plot, possible only if the \code{sizeVar} is 'numeric' or 'integer'
#' @param alphaVar name of variable (in varLabels of the \code{eset}) used for the transparency, NULL by default.
#' This parameter is currently only available for static plot.
#' @param alpha specified transparency(s) for the points, replicated if needed, used only if \code{shapeVar} is NULL, a factor or character
#' by default: '1' if \code{alphaVar} is not specified and default \code{ggplot} alpha otherwise
#' This parameter is currently only available for static plot.
#' @param alphaRange transparency (alpha) range used in the plot, possible only if the \code{alphaVar} is 'numeric' or 'integer'
#' This parameter is currently only available for static plot.
#' @param title plot title, '' by default
#' @param symmetryAxes set symmetry for axes, either:
#' \itemize{
#'  \item{'combine' (by default): }{both axes are symmetric and with the same limits}
#'  \item{'separate': }{each axis is symmetric and has its own limits}
#'  \item{'none': }{axes by default (plot limits)}
#' }
#' @param cloudGenes logical, if TRUE (by default), include the cloud of genes in the spectral map
#' @param cloudGenesColor if \code{cloudGenes} is TRUE, color for the cloud of genes, black by default
#' @param cloudGenesNBins number of bins to used for the clouds of genes,
#' by default the square root of the number of genes
#' @param cloudGenesIncludeLegend logical, if TRUE (FALSE by default) 
#' include the legend for the cloud of genes (in the top position if multiple legends)
#' @param cloudGenesTitleLegend string with title for the legend for the cloud of genes
#' 'nGenes' by default
#' @param packageTextLabel package used to label the outlying genes/samples/gene sets,
#' either \code{ggrepel} (by default, only used if package \code{ggrepel} is available),
#' or \code{ggplot2}
#' @param topGenes numeric indicating which percentile (if <=1) or number (if >1) of genes
#' most distant to the origin of the plot to annotate, by default: 10 genes are selected
#' If no genes should be annotated, set this parameter to 0
#' Currently only available for static plot.
#' @param topGenesCex cex for gene annotation (used when \code{topGenes} > 0)
#' @param topGenesVar variable of the featureData used to label the genes, 
#' by default: NULL, the featureNames are used for labelling (used when \code{topGenes} > 0)
#' @param topGenesJust text justification for the genes (used when \code{topGenes} > 0 and if \code{packageTextLabel} is \code{ggplot2}),
#' by default: c(0.5, 0.5) so centered
#' @param topGenesColor text color for the genes (used when \code{topGenes} > 0), black by default
#' @param topSamples numeric indicating which percentile (if <=1) or number (if >1) of samples
#' most distant to the origin of the plot to annotate, by default: 10 samples are selected
#' If no samples should be annotated, set this parameter to 0.
#' Currently available for static plot.
#' @param topSamplesCex cex for sample annotation (used when \code{topSamples} > 0)
#' @param topSamplesVar variable of the phenoData used to label the samples, 
#' by default: NULL, the sampleNames are used for labelling (used when \code{topSample}s > 0)
#' @param topSamplesJust text justification for the samples (used when \code{topSamples} > 0 and if \code{packageTextLabel} is \code{ggplot2}),
#' by default: c(0.5, 0.5) so centered
#' @param topSamplesColor text color for the samples (used when \code{topSamples} > 0), black by default
#' @param geneSets list of gene sets/pathways, each containing identifiers of genes contained in the set.
#' E.g. pathways from Gene Ontology databases output from the \code{\link{getGeneSetsForPlot}} function or any custom list of pathways.
#' The genes identifiers should correspond to the variable \code{geneSetsVar} contained in the phenoData, if not specified
#' the featureNames are used.
#' If several gene sets have the same name, they will be combine to extract the top gene sets.
#' @param geneSetsVar variable of the featureData used to match the genes contained in geneSets,
#' most probably ENTREZID, if not specified the featureNames of the eSet are used
#' Only used when \code{topGeneSets} > 0 and the parameter geneSets is specified.
#' @param geneSetsMaxNChar maximum number of characters for pathway names, by default keep entire names
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
#' Only used when \code{topGeneSets} > 0, the parameter \code{geneSets} is specified and if \code{packageTextLabel} is \code{ggplot2}.
#' @param topGeneSetsColor color for the gene sets (used when \code{topGeneSets} > 0 and \code{geneSets} is specified), black by default
#' Only used when \code{topGeneSets} > 0 and the parameter geneSets is specified.
#' @param includeLegend logical if TRUE (by default) include a legend, otherwise not
#' @param includeLineOrigin if TRUE (by default) include vertical line at x = 0 and horizontal line at y = 0
#' @param typePlot type of the plot returned, either 'static' (static) or interactive' (potentially interactive)
#' @param figInteractiveSize vector containing the size of the interactive plot, as [width, height]
#' by default: c(600, 400). This is passed to the \code{width} and \code{height} parameters of:
#' \itemize{
#'  \item{for rbokeh plots: }{the \code{bokeh::figure} function}
#'  \item{for ggvis plots: }{the \code{ggvis::set_options} function}
#' }
#' @param ggvisAdjustLegend logical, if TRUE (by default) adjust the legends in \code{ggvis} to avoid
#' overlapping legends when multiple legends
#' @param interactiveTooltip logical, if TRUE, add hoover functionality showing
#' sample annotation (variables used in the plot) in the plot
#' @param interactiveTooltipExtraVars name of extra variable(s) (in varLabels of the \code{eset}) to add in tooltip to label the samples,
#' NULL by default
#' @param packageInteractivity if \code{typePlot} is 'interactive', package used for interactive plot,
#' either 'rbokeh' (by default) or 'ggvis'
#' @param returnTopElements logical, if TRUE (FALSE by default) return the outlying elements labelled in the plot (if any)
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
#'       \item{the top elements labelled in the plot}
#' 	     \item{the \code{ggplot} object}
#'       }
#'     \item{otherwise, the \code{ggplot} object only}
#'    }
#'   \item{\code{interactive}: }{a \code{ggvis} or \code{rbokeh} object, depending on the \code{packageInteractivity} parameter}
#' }
#' @import Biobase
#' @importFrom grid unit
#' @importFrom hexbin hexbin
#' @importFrom grDevices colorRamp colorRampPalette
#' @importFrom utils getFromNamespace
#' @author Laure Cougnaud
#' @export
esetPlotWrapper <- function(
	dataPlotSamples, dataPlotGenes = NULL, esetUsed, 
	xlab = "", ylab = "",
	colorVar = NULL, #color = "black", colorVarValues = NULL, 
	color = if(is.null(colorVar))	"black"	else	NULL,
	shapeVar = NULL, 
	shape = if(is.null(shapeVar))	15	else	NULL,
	sizeVar = NULL, 
	size = if(is.null(sizeVar))	2.5	else	NULL, 
	sizeRange = NULL, #c(1, 6),
	alphaVar = NULL, 
	alpha = if(is.null(alphaVar))	1	else	NULL, 
	alphaRange = NULL,
	title = "", 
	#assume that data are already log-transformed
	symmetryAxes = c("combine", "separate", "none"),
	cloudGenes = TRUE, cloudGenesColor = "black",
	cloudGenesNBins = if(!is.null(dataPlotGenes))	sqrt(nrow(dataPlotGenes))	else	NULL,
	cloudGenesIncludeLegend = FALSE, cloudGenesTitleLegend = "nGenes",
	packageTextLabel = c("ggrepel", "ggplot2"),
	topGenes = 10, topGenesCex = 2.5, topGenesVar = NULL, topGenesJust = c(0.5, 0.5), topGenesColor = "black",
	topSamples = 10, topSamplesCex = 2.5, topSamplesVar = NULL, topSamplesJust = c(0.5, 0.5), topSamplesColor = "black",
	geneSets = list(), geneSetsVar = NULL, geneSetsMaxNChar = NULL, topGeneSets = 10, 
	topGeneSetsCex = 2.5, topGeneSetsJust = c(0.5, 0.5), topGeneSetsColor = "black",
	includeLegend = TRUE, includeLineOrigin = TRUE,
	typePlot = c("static", "interactive"),
	figInteractiveSize  = c(600, 400), ggvisAdjustLegend = TRUE,
	interactiveTooltip = TRUE, interactiveTooltipExtraVars = NULL,
	packageInteractivity = c("rbokeh", "ggvis"),
	returnTopElements = FALSE){

	symmetryAxes <- match.arg(symmetryAxes)
	packageTextLabel <- match.arg(packageTextLabel)
	packageInteractivity <- match.arg(packageInteractivity)
	typePlot <- match.arg(typePlot)
	
	# get methods depending on the class of the object
	esetMethods <- getMethodsInputObjectEsetVis(esetUsed)
	
	# add needed variables in the input data for the plot
	## Bug #7160 fix - added drop = FALSE to handle cases with only one 'Var'
	sampleAnnotation <- esetMethods$pData(esetUsed)[, c(colorVar, shapeVar, sizeVar, alphaVar), drop = FALSE]
	if(typePlot == "interactive" & interactiveTooltip & packageInteractivity == "ggvis")	
		sampleAnnotation <- cbind(sampleAnnotation, keyggvis = esetMethods$sampleNames(esetUsed))
	dataPlotWithAnnotation <- cbind(dataPlotSamples, sampleAnnotation)
	
	colnames(dataPlotWithAnnotation) <- c(colnames(dataPlotSamples), c(colorVar, shapeVar, sizeVar, alphaVar))
		
	## functions common between the different type of plots available
	
	# set fixed aesthetic (e.g. color, shape, size 'palette')
	setFixElement <- function(typeVar, valVar)	is.null(typeVar) & !is.null(valVar)
	
	# set manual aesthetic
	#only if variable, values are specified and if the variable is not numeric or integer (doesn't work with ggplot2)
	setManualScaleTest <- function(typeVar, valVar)	!is.null(typeVar) & !is.null(valVar) & 
		!class(dataPlotWithAnnotation[, typeVar]) %in% c("numeric", "integer")
	formatManualScale <- function(valVar, nameVar){
		values <- rep(valVar, length.out = nlevels(factor(sampleAnnotation[, nameVar])))
		names(values) <- NULL #cannot provide named argument for colors
		values
	}
	
	# get axes limits, depending on the 'symmetryAxes' parameter
	getAxesLimits <- function(){
		#if no genes plotted, only samples, otherwise use genes
		dataXY <- if((cloudGenes | topGenes > 0) & !is.null(dataPlotGenes))	
			rbind(dataPlotSamples[, c('X', 'Y')], dataPlotGenes[, c('X', 'Y')])	else	dataPlotSamples[, c('X', 'Y')]
		#get maximum by axis
		maxCoordByAxis <- apply(abs(dataXY), 2, max)
		getAxisLimit <- function(x) c(-x, x)
		#define axes limits
		axesLimits <- switch(symmetryAxes, 
			'separate' = {res <- sapply(maxCoordByAxis, getAxisLimit); 
				colnames(res) <- c("x", "y"); res},
			'combine' = sapply(c("x", "y"), function(x) getAxisLimit(max(maxCoordByAxis))),
			'none' = {res <- apply(dataXY, 2, range);colnames(res) <- c("x", "y");res})
	}

	plotGeneSets <- 
		if(length(geneSets) > 0 & topGeneSets > 0){
			if(is.null(geneSetsVar)){
				warning(paste("No gene sets are plotted because the variable",
					"describing the mapping of genes IDs between the 'geneSets' and the 'eset'",
					"objects is not specified in the 'geneSetsVar' parameter."
				))
				FALSE
			}else TRUE
		}else FALSE

	## frameplot
	switch(typePlot,
			
		'static' = {
			
			if(!requireNamespace("ggplot2", quietly = TRUE))
				stop(paste("The package 'ggplot2' need to be loaded to create static plots."))
			
			g <- ggplot2::ggplot()
			
			## gene plot first
			if(cloudGenes & !is.null(dataPlotGenes)){
				
				# implementation with alpha, but cannot then used alpha for the samples
#				g <- g + stat_binhex(aes_string(x = 'X', y = 'Y', alpha = '..count..'), 
#					data = dataPlotGenes, fill = cloudGenesColor, bins = cloudGenesNBins)#, colour = "white")
#				g <- g + guides(fill = FALSE, alpha = FALSE) 
				
				# implementation for fill
				baseFillColor <- do.call("rgb", c(as.list(
					colorRamp(c("white", cloudGenesColor))(0.2)), list(maxColorValue = 255))
				)
				
				# because parameter/function names changed in ggplot2
				if (utils::packageVersion("ggplot2") < "2.0.0"){
					fillAes <- '..count..'
					statBinHexFct <- ggplot2::stat_binhex
				}else{
					# ggplot2 2.0.0: stat_binhex() has been renamed to stat_bin_hex()
					statBinHexFct <- ggplot2::stat_bin_hex
					# ggplot2 2.1.0: output of stat_bin_hex() is now value instead of count.
					fillAes <- if (utils::packageVersion("ggplot2") < "2.1.0")	'..count..'	else	'..value..'
				}
				
				g <- g + statBinHexFct(ggplot2::aes_string(x = 'X', y = 'Y', fill = fillAes),
					data = dataPlotGenes, bins = cloudGenesNBins) +
					ggplot2::scale_fill_gradientn(colours = c(baseFillColor, cloudGenesColor))#c(0.5, 0.8)
				#if(includeLegend)	nLegendsSamples <- length(c(colorVar, shapeVar, sizeVar, alphaVar))
				# set legend if specified, and in the first position
				g <- g + ggplot2::guides(
					fill = if(!cloudGenesIncludeLegend)
						FALSE	else	
						ggplot2::guide_legend(order = 1, 
							# change legend title, by default 'count'
							title = cloudGenesTitleLegend))
					
			}
			
			## sample plot
			aesArgSamplePlot <- c(list(x = 'X', y = 'Y'),
				if(!is.null(colorVar))	list(color = colorVar),
				if(!is.null(shapeVar))	list(shape = shapeVar),
				if(!is.null(sizeVar))	list(size = sizeVar),
				if(!is.null(alphaVar))	list(alpha = alphaVar)
			)
						
			geomPointArgs <- c(
				list(data = dataPlotWithAnnotation,
					mapping = do.call(ggplot2::aes_string, aesArgSamplePlot)), 
				if(setFixElement(colorVar, color))	list(color = color),
				if(setFixElement(shapeVar, shape))	list(shape = shape),
				if(setFixElement(sizeVar, size))	list(size = size),
				if(setFixElement(alphaVar, alpha))	list(alpha = alpha)
			)
			
			g <- g + do.call(ggplot2::geom_point, geomPointArgs)
			
			#horizontal/vertical lines
			if(includeLineOrigin){
				g <- g + ggplot2::geom_vline(xintercept = 0, linetype = 'dashed')			  # X axis
				g <- g + ggplot2::geom_hline(yintercept = 0, linetype = 'dashed')			  # Y axis
			}
			
			#manual specifications: custom scales
			#only if variable, values are specified and if the variable is not numeric or integer (doesn't work with ggplot2)
			setManualScaleStatic <- function(typeVar, nameVar, valVar){
				values <- formatManualScale(valVar, nameVar)
				do.call(getFromNamespace(paste("scale", typeVar, "manual", sep = "_"), ns = "ggplot2"),
					list(values = values) )#, name = nameVar, 
			}
			if (setManualScaleTest(colorVar, color))	g <- g + setManualScaleStatic("color", colorVar, color)
			if (setManualScaleTest(shapeVar, shape))	g <- g + setManualScaleStatic("shape", shapeVar, shape)
			if (setManualScaleTest(sizeVar, size))	g <- g + setManualScaleStatic("size", sizeVar, size)
			if (setManualScaleTest(alphaVar, alpha))	g <- g + setManualScaleStatic("alpha", alphaVar, alpha)
						
			## axes 
			
			#labels
		
			#plot axes
			if(xlab != "")	g <- g + ggplot2::xlab(paste0("\n", xlab))
			if(ylab != "")	g <- g + ggplot2::ylab(paste(ylab, "\n"))

			# increase margin between ticks and axes labels
			argsTheme <- c(
				list(
					# labels y axes centered on ticks
					axis.text.y = ggplot2::element_text(vjust = 5),
					# title x axis further away from the axes labels
					axis.title.x = ggplot2::element_text(vjust = -1)
				),
				# code for previous ggplot2 version
				if (utils::packageVersion("ggplot2") < "2.0.0"){
					list(axis.ticks.margin = unit(0.5, "lines"))
				# code for more recent version
				}else{
					list(
						axis.text.x = ggplot2::element_text(margin = unit(0.5, "lines")),
						axis.text.y = ggplot2::element_text(margin = unit(0.5, "lines"))
					)
				}
			)
			g <- g + do.call(ggplot2::theme, argsTheme)
			
			#symmetry
			if(symmetryAxes != "none"){
				
				#define axes limits
				axesLimits <- getAxesLimits()
				#set axes limits
				g <- g + 
					ggplot2::scale_x_continuous(limits = axesLimits[, "x"]) +
					ggplot2::scale_y_continuous(limits = axesLimits[, "y"])
				
			}
			
			#custom size range, works only if size variable is numeric, or integer
			if(class(dataPlotWithAnnotation[, sizeVar]) %in% c("numeric", "integer")
				& !is.null(sizeRange))	g <- g + ggplot2::scale_size(range = sizeRange)
		
			#custom transparency range, works only if alpha variable is numeric, or integer
			if(class(dataPlotWithAnnotation[, alphaVar]) %in% c("numeric", "integer")
				& !is.null(alphaRange))	g <- g + ggplot2::scale_alpha(range = alphaRange)
			
			#add title
			if(title != "")	g <- g + ggplot2::ggtitle(title)
			
			#theme
			g <- g + ggplot2::theme_bw()
			
			## annotation top genes/samples, at the end to avoid overlapping with plot
			
			topElements <- NULL
			
			# top genes
			if(topGenes > 0 & !is.null(dataPlotGenes)){
				outputTopGenes <- plotTopElements(top = topGenes, type = "gene", 
					var = topGenesVar, cex = topGenesCex, just = topGenesJust,
					color = topGenesColor,
					dataPlotGenes = dataPlotGenes, esetUsed =  esetUsed,
					returnTopElements = returnTopElements,
					packageTextLabel = packageTextLabel)
				if(returnTopElements){
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
			if(topSamples > 0){
				outputTopSamples <- plotTopElements(top = topSamples, type = "sample", 
					var = topSamplesVar, cex = topSamplesCex, just = topSamplesJust,
					color = topSamplesColor,
					dataPlotSamples = dataPlotSamples, esetUsed =  esetUsed,
					returnTopElements = returnTopElements,
					packageTextLabel = packageTextLabel)
				if(returnTopElements){
					topElements <- c(topElements, list(topSamples = outputTopSamples$topElements))
					g <- g + outputTopSamples$geomText
				}else	g <- g + outputTopSamples
			}
			
			#top gene sets
			if(plotGeneSets){
				outputTopGeneSets <- plotTopElements(
					top = topGeneSets, type = "geneSets", 
					cex = topGeneSetsCex, just = topGeneSetsJust,
					color = topGeneSetsColor,
					geneSets = geneSets, esetUsed =  esetUsed, 
					dataPlotGenes = dataPlotGenes,
					geneSetsVar = geneSetsVar,
					geneSetsMaxNChar = geneSetsMaxNChar,
					returnTopElements = returnTopElements,
					packageTextLabel = packageTextLabel)
				if(returnTopElements){
					topElements <- c(topElements, list(topGeneSets = outputTopGeneSets$topElements))
					g <- g + outputTopGeneSets$geomText
				}else	g <- g + outputTopGeneSets
			}
		
			if(!includeLegend)	g <- g + ggplot2::theme(legend.position = "none")
		
			#if(!is.null(nColLegend))	guide_legend(ncol = nColLegend)

			res <- if(returnTopElements & typePlot == "static" & !is.null(topElements)){
				list(topElements = topElements, plot = g)
			}else g
		
		},
	
		'interactive' = {
			
			switch(packageInteractivity,
					
				'ggvis' = {
					
					if(!requireNamespace("ggvis", quietly = TRUE))
						stop(paste("The package 'ggvis' need to be loaded to create",
							"interactive plots with ggvis."))
		
					if(!is.null(alphaVar))
						warning("The transparency aesthetic (alpha) is not ",
							"yet implemented for a ggvis interactive plot.")
					
					## sample plot
					getProp <- function(type, typeVar)
						if(!is.null(typeVar))	list(ggvis::prop(type, as.name(typeVar)))	else	NULL
					
					ggvisArgsSamplePlot <- c(
						list(data = dataPlotWithAnnotation),
						getProp("x", "X"),
						getProp("y", "Y"),
						getProp("fill", colorVar),
						getProp("shape", shapeVar),
						getProp("size", sizeVar)
						#getProp("opacity", alphaVar)
						#if(interactiveTooltip)	list(props(key := ~"keyggvis"))
					)
								
					## gene plot first
					if(cloudGenes & !is.null(dataPlotGenes)){
						
						hexbinGeneData <- hexbin(dataPlotGenes$X, dataPlotGenes$Y, xbins = cloudGenesNBins)
						hexbinGeneDataDf <- data.frame(xGene = hexbinGeneData@xcm, yGene = hexbinGeneData@ycm, geneCol = hexbinGeneData@count)
						
						g <- ggvis::ggvis(x = ~xGene, y= ~yGene, data = hexbinGeneDataDf)
						g <- ggvis::layer_points(vis = g, fillOpacity = ~geneCol, fill = cloudGenesColor)
						g <- ggvis::hide_legend(vis = g, "fillOpacity") 
						g <- ggvis::scale_numeric(vis = g, property = "opacity", trans = "sqrt")
					
		#				g %>% layer_points(data = dataPlotWithAnnotation, 
		#					x=~X, y=~Y, fill =~getPropcolorVar, shape =~ shapeVar, size =~sizeVar)
						
						# add sample plot to the gene plot
						ggvisArgsSamplePlotWithGenePlot <- c(list(vis = g), ggvisArgsSamplePlot)
						g <- do.call(ggvis::layer_points, ggvisArgsSamplePlotWithGenePlot)
						
					}else	g <- ggvis::layer_points(vis = do.call(ggvis::ggvis, ggvisArgsSamplePlot))
								
					#TODO: horizontal/vertical lines
		#			g <- g + geom_vline(xintercept = 0, linetype = 'dashed')			  # X axis
		#			g <- g + geom_hline(yintercept = 0, linetype = 'dashed')			  # Y axis
					
					# add fixed elements
					if(setFixElement(colorVar, color))
						g <- ggvis::layer_points(g, fill := color)
					if(setFixElement(shapeVar, shape))
						g <- ggvis::layer_points(g, shape := shape)
					if(setFixElement(sizeVar, size))
						g <- ggvis::layer_points(g, size := size)
					
					setManualScaleGgvis <- function(g, typeVar, nameVar, valVar = NULL, typeScale = "logical"){
						values <- if(!is.null(valVar))	formatManualScale(valVar, nameVar)	else NULL
						do.call(
							getFromNamespace(paste0("scale_", typeScale), ns = "ggvis"), 
							list(vis = g, property = typeVar, range = values))
					}	
			
					#manual specifications: custom scales
					if (setManualScaleTest(colorVar, color))	
						g <- setManualScaleGgvis(g, typeVar = "fill", nameVar = colorVar, valVar = color)
					if (setManualScaleTest(shapeVar, shape))	
						g <- setManualScaleGgvis(g, typeVar = "shape", nameVar = shapeVar, valVar = shape)
					if (setManualScaleTest(sizeVar, size))	
						g <- setManualScaleGgvis(g, typeVar = "size", nameVar = sizeVar, valVar = size)
					
					#plot axes labels and title
					g <- ggvis::add_axis(vis = g, "x", title = xlab); g <- ggvis::add_axis(vis = g, "y", title = ylab)
					
					#symmetry
					if(symmetryAxes != "none"){
						
						#define axes limits
						axesLimits <- getAxesLimits()
						#set axes limits
						g <- ggvis::scale_numeric(vis = g, "x", trans = "linear", domain = axesLimits[, "x"])
						g <- ggvis::scale_numeric(vis = g, "y", trans = "linear", domain = axesLimits[, "y"])

					}
					
					#TODO: custom size range, works only if size variable is numeric, or integer
		#			if(class(dataPlotWithAnnotation[, sizeVar]) %in% c("numeric", "integer")
		#				& !is.null(sizeRange))	g <- g + scale_size(range = sizeRange)
					
					#add title
					if (title != ""){
						g <- ggvis::add_axis(
							vis = g, "x", orient = "top", ticks = 0, title = title,
							properties = ggvis::axis_props(
								axis = list(stroke = "white"),
								labels = list(fontSize = 2)
							)
						)
					}
					
					if(length(figInteractiveSize) != 2){
						figInteractiveSize <- c(600, 400)
						warning(paste("The size of the ggvis window should contain two",
							"elements: width and height. The default values c(600, 400) are used."))
						# default options
						# str(default_options())
					}
					
					g <- ggvis::set_options(vis = g, width = figInteractiveSize[1], height = figInteractiveSize[2])
		
		#			orderLegendLog <- c(!is.null(colorVar), !is.null(shapeVar), !is.null(sizeVar))
		#			names(orderLegendLog) <- c("fill", "shape", "size")
		#			if(sum(orderLegendLog) > 0){
		#				legendToSet <- names(orderLegendLog)[orderLegendLog]
		#				argsAddLegend <- c(as.list(legendToSet), list(orient = "right"))
		#				g <- g %>% do.call("add_legend", argsAddLegend)
		#			}
					
					#g <- g %>% set_options(duration = 0)
		
					# adjust legend manually because overlap if several present
					# only if the window size is numeric (not auto)
					# this will be fixed in future version of ggvis? 
		
					orderLegendLog <- c(!is.null(colorVar), !is.null(shapeVar), !is.null(sizeVar))
					names(orderLegendLog) <- c("fill", "shape", "size")
					legendToSet <- names(orderLegendLog)[orderLegendLog]
					
					if(is.numeric(figInteractiveSize) & includeLegend & ggvisAdjustLegend & length(legendToSet) > 0){
						
						# for legend side by side in y direction
						setLegendPos <- function(vis, typeVar, y)
							ggvis::add_legend(vis = vis, scales = typeVar, 
								properties = ggvis::legend_props(legend = list(y = y)))
					
						for(i in 1:length(legendToSet))
							g <- setLegendPos(vis = g, typeVar = legendToSet[i],
								y = (i-1) * figInteractiveSize[2]/(length(legendToSet) + 1))
						
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
					if((!includeLegend) & length(legendToSet) > 0)	
						g <- ggvis::hide_legend(vis = g, scales = legendToSet)
				
					if(interactiveTooltip){
						
						interactiveTooltipExtraVars <- 
							interactiveTooltipExtraVars[interactiveTooltipExtraVars %in% esetMethods$varLabels(esetUsed)]
					
						if(length(interactiveTooltipExtraVars) == 0)
							interactiveTooltipExtraVars <- NULL
						
						tooltipLabelsFct <- function(x)
							return(
								if(is.null(x))	NULL	else{
									#str(x);
									# hoover on gene
									if(any(c("xGene", "yGene") %in% names(x)))	NULL	else{
										
										formatCoor <- function(x)
											paste0("(", paste(formatC(x, width = 3), collapse = ", "), ")")
								
										xy <- as.numeric(x[, c("X", "Y")])
										xyST <- formatCoor(xy)
										
										# use bracket rather than subset ta avoid having notes in R CMD check
										dataPoints <- dataPlotSamples[
											round(dataPlotSamples$X, 3) == round(xy[1], 3) & 
												round(dataPlotSamples$Y, 3) == round(xy[2], 3), 
										]
#											subset(dataPlotSamples, 
#											round(X, 3) == round(xy[1], 3) & round(Y, 3) == round(xy[2], 3))
							
										if(nrow(dataPoints) > 0 & !is.null(interactiveTooltipExtraVars)){
											addedAnnot <- esetMethods$pData(esetUsed)[dataPoints[1, "sampleName"], interactiveTooltipExtraVars, drop = FALSE]
											#str(addedAnnot)
											colnames(addedAnnot) <- interactiveTooltipExtraVars
											x <- data.frame(x, addedAnnot)
										}
										#str(x)
								
										extraCols <- colnames(x)[!colnames(x) %in% c("X", "Y")]
									
										if(length(extraCols) > 0){
											sampleAnnot <- paste0(names(x[, extraCols]), ": ", x[, extraCols], collapse = "<br />")
											paste0(xyST, "<br />", sampleAnnot)
										}else	xyST
									}
								}
							)
						
						g <- ggvis::add_tooltip(vis = g, tooltipLabelsFct, "hover") # can use 'click' too
						
						## TODO: annotation top genes/samples, at the end to avoid overlapping with plot
						
					}
					
						
				},#,
				
				# rbokeh implementation
				
				'rbokeh' = {
					
					if(!requireNamespace("rbokeh", quietly = TRUE))
						stop(paste("The package 'rbokeh' need to be loaded to create",
							"interactive plots with rbokeh."))
		
					if(!is.null(alphaVar))
						if(is.factor(esetMethods$pData(esetUsed)[, alphaVar])){
							warning("The transparency aesthetic (alpha) is not ",
								"yet implemented for factors in rbokeh interactive plot,",
								"so no transparency is used.")
							alphaVar <- NULL; alpha <- 1
						}
					
					#define axes limits
					axesLimits <- getAxesLimits()
		
					## create empty figure
					argsFigure <- c(
						if(is.numeric(figInteractiveSize)) list(width = figInteractiveSize[1], height = figInteractiveSize[2])	else	NULL,
						if(title != "") list(title = title)	else	NULL,
						if(xlab != "")	list(xlab = xlab)	else	NULL,
						if(ylab != "")	list(ylab = ylab)	else	NULL,
						if(!is.null(axesLimits))	list(xlim = axesLimits[, "x"], ylim = axesLimits[, "y"]),
						list(xgrid = FALSE, ygrid = FALSE)
					)
					#if(is.null(argsFigure))	argsFigure <- list(NULL)
					g <- do.call(rbokeh::figure, argsFigure)
					
					# add axes limits
					
					## gene plot first
					if(cloudGenes & !is.null(dataPlotGenes)){
						
						# can only specify color as palette for color ramp: don't use white at lower palette
						paletteFctUsed <- colorRampPalette(c(colorRampPalette(c("white", cloudGenesColor))(10)[2], cloudGenesColor))
						g <- rbokeh::ly_hexbin(fig = g, data = dataPlotGenes, x = "X", y = "Y",
							xbins = cloudGenesNBins, 
							palette = paletteFctUsed, trans = sqrt)#, alpha = 0.8
						
					}
					
					## samples plot
					
					# need to remove NA values for glyph?
					idxRowsKept <- rowSums(is.na(dataPlotWithAnnotation[, c(colorVar, shapeVar, sizeVar), drop = FALSE])) == 0
					dataPlotWithAnnotationWthtNA <- dataPlotWithAnnotation[idxRowsKept, ]
					
					# possible to specify a data.frame for the hoover, 
					# so add additional variables if requested

					# remove variable if already present in the data
					# (so should have the same name several times)
					if(!is.null(interactiveTooltipExtraVars)){
						interactiveTooltipExtraVars <- interactiveTooltipExtraVars[
							!interactiveTooltipExtraVars %in%  colnames(dataPlotWithAnnotationWthtNA)]
						if(length(interactiveTooltipExtraVars) == 0) interactiveTooltipExtraVars <- NULL
					}
				
					hoverDf <- if(interactiveTooltip)	
						if(!is.null(interactiveTooltipExtraVars)){
							hoverDf <- data.frame(dataPlotWithAnnotationWthtNA, 
								esetMethods$pData(esetUsed)[as.character(dataPlotWithAnnotationWthtNA$sampleName), interactiveTooltipExtraVars]
							)
							colnames(hoverDf) <- c(colnames(dataPlotWithAnnotationWthtNA), interactiveTooltipExtraVars)
							hoverDf
						}else	dataPlotWithAnnotationWthtNA	else	NULL
					# remove redundant column
					hoverDf <- as.data.frame(t(unique(t(hoverDf))))
					
					# samples plot
					
					color <- if(!is.null(colorVar))	colorVar	else	color
					glyph <- if(!is.null(shapeVar))	shapeVar	else	shape
					size <- if(!is.null(sizeVar))	sizeVar	else	size
					alpha <- if(!is.null(alphaVar))	alphaVar	else	alpha
					g <- rbokeh::ly_points(
						fig = g,
						data = dataPlotWithAnnotationWthtNA,
						x = "X", y = "Y",
						hover = hoverDf,
						#lname = "sample",
						legend = includeLegend,
						color = color,
						glyph = glyph,
						size = size,
						alpha = alpha
					)
					
					## add horizontal/vertical lines
					if(includeLineOrigin)
						g <- rbokeh::ly_abline(fig = g, v = 0, type = 2)
						g <- rbokeh::ly_abline(fig = g, h = 0, type = 2)
						
					# TODO: add custom palettes when will be available in rbokeh

				})
			
			#)

		res <- g
					
	})

	# Return
	return(res)
	
}
