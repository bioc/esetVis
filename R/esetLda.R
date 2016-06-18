#' @title plot a biplot of a linear discriminant analysis of an \linkS4class{eSet} object
#' @description \code{esetLda} reduces the dimension of the data contained in the \linkS4class{eSet} via a linear discriminant analysis
#' on the specified grouping variable with the \code{lda} function and plot the subsequent biplot,
#'  possibly with sample annotation and gene annotation contained in the eSet.
#' @param eset expressionSet (or SummarizedExperiment) object with data
#' @param ldaVar name of variable (in varLabels of the \code{eset}) used for grouping for lda, NULL by default
#' @param psids featureNames of genes to include in the plot, all by default
#' @param dim dimensions of the analysis to represent, first two dimensions by default
#' @param returnAnalysis logical, if TRUE (FALSE by default), return also the output of the analysis,
#' and the outlying samples in the topElements element if any, otherwise only the plot object
#' @inheritParams esetPlotWrapper
#' @references Fisher, R. A. (1936). The Use of Multiple Measurements in
#' Taxonomic Problems. Annals of Eugenics, 7 (2), 179--188
#' @seealso the function used internally: \code{\link[MASS]{lda}}
#' @return if \code{returnAnalysis} is TRUE, return a list:
#' \itemize{
#'  \item{analysis: }{output of the spectral map analysis, can be given as input to the \code{\link{esetPlotWrapper}} function}
#'    \itemize{
#' 		\item{dataPlotSamples: }{coordinates of the samples}
#' 		\item{dataPlotGenes: }{coordinates of the genes}
#' 		\item{esetUsed: }{expressionSet used in the plot}
#'      \item{axisLabels: }{axes labels indicating percentage of variance explained by the selected axes}
#' 	  }
#'   \item{topElements: }{list with top outlying elements if any, possibly genes, samples and gene sets}
#'   \item{plot: }{the plot output}
#' }
#' otherwise return only the plot
#' @examples
#' # load data
#' library(ALL)
#' data(ALL)
#' 
#' # specify several variables in ldaVar (this might take a few minutes to run...)
#' 
#' # sample subsetting: currently cannot deal with missing values
#' samplesToRemove <- which(apply(pData(ALL)[, c("sex", "BT")], 1, anyNA)) 
#' 
#' # extract random features, because analysis is quite time consuming
#' retainedFeatures <- sample(featureNames(ALL), size = floor(nrow(ALL)/5))
#' 
#' # create the plot
#' esetLda(eset = ALL[retainedFeatures, -samplesToRemove], 
#'   ldaVar = "BT", colorVar = "BT", shapeVar = "sex", sizeVar = "age",
#'   title = "Linear discriminant analysis on the ALL dataset")
#' @import Biobase
#' @importFrom stats cor
#' @importFrom MASS lda
#' @author Laure Cougnaud
#' @export
esetLda <- function(eset, ldaVar,
	psids = 1:nrow(eset),
	dim = c(1, 2),
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
	symmetryAxes = c("combine", "separate", "none"),
	packageTextLabel = c("ggrepel", "ggplot2"),
	cloudGenes = TRUE, cloudGenesColor = "black", cloudGenesNBins = sqrt(length(psids)), 
	cloudGenesIncludeLegend = FALSE, cloudGenesTitleLegend = "nGenes",
	topGenes = 10, topGenesCex = 2.5, topGenesVar = NULL, topGenesJust = c(0.5, 0.5), topGenesColor = "black",
	topSamples = 10, topSamplesCex = 2.5, topSamplesVar = NULL, topSamplesJust = c(0.5, 0.5), topSamplesColor = "black",
	geneSets = list(), geneSetsVar = NULL, geneSetsMaxNChar = NULL, topGeneSets = 10, 
	topGeneSetsCex = 2.5, topGeneSetsJust = c(0.5, 0.5), topGeneSetsColor = "black",
	includeLegend = TRUE, includeLineOrigin = TRUE,
	typePlot = c("static", "interactive"), packageInteractivity = c("rbokeh", "ggvis"),
	figInteractiveSize  = c(600, 400), ggvisAdjustLegend = TRUE, interactiveTooltip = TRUE, interactiveTooltipExtraVars = NULL,
	returnAnalysis = FALSE){
	
	packageTextLabel <- match.arg(packageTextLabel)
	symmetryAxes <- match.arg(symmetryAxes)
	packageInteractivity <- match.arg(packageInteractivity)
	
	if(length(psids) <= 1)
		stop(paste("At least two genes should be used for the visualization.",
			"Please change the 'psids' argument accordingly."))

	# get methods depending on the class of the object
	esetMethods <- getMethodsInputObjectEsetVis(eset)

	## Extract exprsMat with specified psids
	esetUsed <- eset[psids, ]
	
	errorNAGroupingVariable <- paste("The current implementation cannot deal with",
		"missing values in the grouping variable.",
		"Please remove the corresponding samples from the data."
	)
	
	ldaVariable <- if(length(ldaVar) > 1){
		mesCombinedVar <- "The variables used for lda are combined."
		ldaVariableDf <- esetMethods$pData(esetUsed)[, ldaVar] # pData(esetUsed)
		ldaVariable <- apply(ldaVariableDf, 1, paste, collapse = "-")
		idxNA <- rowSums(is.na(ldaVariableDf)) > 0
		if(sum(idxNA) > 0){
			stop(errorNAGroupingVariable)
#			mesCombinedVar <- paste(mesCombinedVar, 
#				sum(idxNA), "samples with NA in at least one of the variable are set to NA in the combined variable.")
#			ldaVariable[idxNA] <- NA
		}
		message(mesCombinedVar);factor(ldaVariable)
	}else esetMethods$pData(esetUsed)[, ldaVar]
	
	if(!is.factor(ldaVariable))
		stop("The grouping variable for the lda ",
			"should be a factor.")
	
	if(nlevels(ldaVariable) <= 2)
		stop("The current function only deal with ",
			"a grouping variable with at least 3 levels.")

	if(anyNA(ldaVariable))	stop(errorNAGroupingVariable)
		
	## Run lda
	inputLda <- t(esetMethods$exprs(esetUsed)) # exprs(esetUsed)
	system.time(ldaOutput <- lda(inputLda, 
		grouping = ldaVariable))#MASS package

	## reformat results of lda function
	propTraces <- round((ldaOutput$svd^2)/sum(ldaOutput$svd^2)*100)

	if(max(dim) > length(propTraces)){
		warning("The dimensions to represent should be less than",
			" the number of levels in the grouping variable minus 1.",
			" The dimensions c(1, 2) will be used.")
		dim <- c(1, 2)
	}
		
	labelsAxes <- paste0("LD", dim, " (", propTraces[dim],"%)")
	
	#scores of the variables (here genes)
	scoreVar <- ldaOutput$scaling[, dim]
	
	#coordinates of the individuals
	means <- colMeans(ldaOutput$means)
	coordInd <- scale(inputLda, center = means, scale = FALSE) %*% scoreVar
	
	#compute correlation between coordinates of samples and input data
	corrVarInit <- cor(inputLda, coordInd, use = "pairwise")
	# scale it to the maximum limits
	corrVar <- t(t(corrVarInit) * apply(coordInd, 2, max))
#	mostInfluentGenes <- names(sort(apply(corrVar , 1, function(x) sum(abs(x)^2)), decreasing = TRUE)[1:nTopGenes])
#	names(mostInfluentGenes) <- fData(eset)[mostInfluentGenes, "SYMBOL"]
#	corrVarWithoutMostInfluentGenes <- corrVar[!rownames(corrVar) %in% mostInfluentGenes, ]

	## Extract data for final plot
	dataPlotSamples <- data.frame(coordInd, rownames(coordInd))
	colnames(dataPlotSamples) <- c("X", "Y", "sampleName")
	
	dataPlotGenes <- data.frame(corrVar, rownames(corrVar))
	colnames(dataPlotGenes) <- c("X", "Y", "geneName")
	
	plot <- esetVis::esetPlotWrapper(dataPlotSamples = dataPlotSamples,
		dataPlotGenes = dataPlotGenes, 
		esetUsed = esetUsed, 
		colorVar = colorVar, color = color,
		shapeVar = shapeVar, shape = shape,
		sizeVar = sizeVar, size = size, sizeRange = sizeRange,
		alphaVar = alphaVar, alpha = alpha, alphaRange = alphaRange,
		title = title, symmetryAxes = symmetryAxes,
		cloudGenes = cloudGenes, cloudGenesColor = cloudGenesColor, 
		cloudGenesNBins = cloudGenesNBins, 
		cloudGenesIncludeLegend = cloudGenesIncludeLegend, cloudGenesTitleLegend = cloudGenesTitleLegend,
		topGenes = topGenes, topGenesCex = topGenesCex, 
		topGenesVar = topGenesVar, topGenesJust = topGenesJust, topGenesColor = topGenesColor,
		topSamples = topSamples, topSamplesCex = topSamplesCex, 
		topSamplesVar = topSamplesVar, topSamplesJust = topSamplesJust, topSamplesColor = topSamplesColor,
		geneSets = geneSets, geneSetsVar = geneSetsVar, geneSetsMaxNChar = geneSetsMaxNChar, 
		topGeneSets = topGeneSets, topGeneSetsCex = topGeneSetsCex, topGeneSetsJust = topGeneSetsJust,
		topGeneSetsColor = topGeneSetsColor,
		includeLegend = includeLegend, includeLineOrigin = includeLineOrigin,
		typePlot = typePlot, packageInteractivity = packageInteractivity,
		figInteractiveSize = figInteractiveSize, 
		ggvisAdjustLegend = ggvisAdjustLegend, interactiveTooltip = interactiveTooltip, interactiveTooltipExtraVars = interactiveTooltipExtraVars,
		returnTopElements = returnAnalysis,
		packageTextLabel = packageTextLabel)

	someTopElementsReturned <- any(!class(plot) == "ggplot")

	res <- if(!returnAnalysis)	plot else
		c(
			list(analysis = list(
				dataPlotSamples = dataPlotSamples,
				dataPlotGenes = dataPlotGenes, 
				esetUsed = esetUsed
				#axisLabels = axisLabels
			)),
		if(someTopElementsReturned)	list(topElements = plot$topElements)	else NULL,
		list(plot = plot$plot)
		)

	return(res)

}
