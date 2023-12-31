#' plot top elements for a static plot
#'
#' This create \code{geom_text} object with top genes/sample/pathways
#' @param packageTextLabel package used to label the outlying genes/samples/gene sets,
#' either 'ggrepel' (by default, only used if package \code{ggrepel} is available),
#' or 'ggplot2'
#' @param cex cex of text in the plot
#' @param just justification of elements in the plot, only use if \code{packageTextLabel} is 'ggplot2'
#' @param color color for the elements in the plot
#' @param returnTopElements logical if TRUE (FALSE by default) return the outlying elements
#' @inheritParams getTopElements
#' @return 
#' \itemize{
#'  \item{if the \code{elements} are present in the data: if \code{returnTopElements} is: }
#'	 \itemize{
#'    \item{TRUE: }{return a list with two arguments:} 
#'      \itemize{
#'       \item{topElements: }{string with top elements labelled in the plot}
#'       \item{geomText: }{output of \code{geom_text}}
#'      }
#'     \item{FALSE: }{only return the output of \code{geom_text}} 
#' 	 }
#'  \item{if not, return \code{NULL}}
#' }
#' @keywords internal
#' @author Laure Cougnaud
plotTopElements <- function(
	packageTextLabel = c("ggrepel", "ggplot2"),
	cex = 1, just = c(0.5, 0.5), color = "black",
	returnTopElements = FALSE, ...){

	packageTextLabel <- match.arg(packageTextLabel)
	dataPlotText <- getTopElements(...)

	# issue if all labels are NA
	geomText <- if(all(is.na(dataPlotText$labels)))	NULL	else{

		argsGeomTextFct <- list(
			mapping = ggplot2::aes_string(x = 'X', y = 'Y', label = 'labels'),
			data = dataPlotText, color = color, size = cex
		)
		if(packageTextLabel	== "ggrepel" & requireNamespace("ggrepel", quietly = TRUE)){
			geomTextFct <- ggrepel::geom_text_repel
		}else{
			geomTextFct <- ggplot2::geom_text
			argsGeomTextFct <- c(argsGeomTextFct, list(hjust = just[1], vjust = just[2]))
		}
		do.call(geomTextFct, argsGeomTextFct)

	}

	res <- if(returnTopElements){

		topElements <- dataPlotText$labels
		names(topElements) <- rownames(dataPlotText)
		list(topElements = topElements, geomText = geomText)

	}else	geomText

	return(res)

}

#' create \code{geom_text} object with top genes/sample/pathways
#' @param top numeric, number of top elements
#' @param type type of elements to plot, either 'gene', 'sample', or 'geneSets'
#' @param var variable used to annotate the elements, only used for 'gene' and 'sample'
#' @param dataPlotGenes data.frame with two columns 'X' and 'Y' with coordinates for the genes
#' @param dataPlotSamples data.frame with two columns 'X' and 'Y' with coordinates for the samples
#' @param esetUsed expressionSet (or SummarizedExperiment) object with data
#' @param geneSets list of gene sets, e.g. gene pathways, output from the 'getGeneSets' function in MLP
#' the genes IDs must correspond to the sampleNames in the eset.
#' If several gene sets have the same name, they will be combine to extract the top gene sets.
#' @param geneSetsVar variable of the featureData used to match the genes contained in geneSets,
#' most probably ENTREZID, if not specified the featureNames of the eSet are used
#' @param geneSetsMaxNChar maximum number of characters for pathway names, by default keep entire names
#' @param returnTopElements logical if TRUE (FALSE by default) return the outlying elements
#' @return Data.frame with coordinates and labels of the top elements
#' @keywords internal
#' @author Laure Cougnaud
getTopElements <- function(top,
	type = c("gene", "sample", "geneSets"),
	var = character(),
	dataPlotGenes = data.frame(), dataPlotSamples = data.frame(), esetUsed, 
	geneSets = list(), geneSetsVar = character(), geneSetsMaxNChar = numeric()){

	type <- match.arg(type)
	
	# get methods depending on the class of the object
	esetMethods <- getMethodsInputObjectEsetVis(esetUsed)
	
	switch(type,
			
		'geneSets' = if(length(geneSets) == 0|nrow(dataPlotGenes) == 0)
			stop("`geneSets' and some 'dataPlotGenes' should be specified."),
		
		'gene' = if(nrow(dataPlotGenes) == 0)
			stop("'dataPlotGenes' should be specified."),
		
		'sample' = if(nrow(dataPlotSamples) == 0)
			stop("'dataPlotSamples' should be specified.")

	)
	
	coor <- switch(type, 
			
		#if gene sets, take the mean coordinates for each gene set
		#TODO: improve speed
		'geneSets' = getCoordGeneSets(
			dataPlotGenes = dataPlotGenes, geneSets = geneSets, 
			esetUsed = esetUsed, geneSetsVar = geneSetsVar
		),
		
		#if gene or samples, take directly the coordinates
		switch(type, 'sample' = dataPlotSamples, 'gene' = dataPlotGenes)[, c("X", "Y")]

	)
	
	# for gene sets, if no mapped gene is found
	dataPlotText <- if(nrow(coor) > 0){
	
		nElements <- nrow(coor)
		distToOrigin <- sqrt(rowSums(coor ^ 2))
		idxElementsSorted <- order(distToOrigin, decreasing = TRUE)
		
		#if top less than 1, percentage, otherwise absolute number
		idxTop <- 1:(if(top < 1)	(top * nElements)	else	min(top, nElements))
		idxElementsKept <- idxElementsSorted[idxTop]
		coorElementsKept <- coor[idxElementsKept, ]
	
		labels <- if(type != "geneSets"){
			varInAnnot <- ifelse(length(var) == 0, "", var) %in% 
				switch(type, 'gene' = esetMethods$fvarLabels, 'sample' = esetMethods$varLabels)(esetUsed)
			if(!varInAnnot)	rownames(coorElementsKept)	else
				switch(type, 'gene' = esetMethods$fData, 'sample' = esetMethods$pData)(esetUsed)[
					rownames(coorElementsKept), var]
		}else{
			if(length(geneSetsMaxNChar) > 0){
				res <- rownames(coorElementsKept)
				names(res) <- substr(rownames(coorElementsKept), 1, geneSetsMaxNChar)
				res
			}else	rownames(coorElementsKept)
		}
		
		data.frame(coor[idxElementsKept, ], labels, stringsAsFactors = FALSE)
		
	}else{
		warning(paste("No labels for the", type, "are found."))
		NULL
	}

	return(dataPlotText)
	
}

#' extract coordinates gene sets
#' @param dataPlotGenes data.frame with two columns 'X' and 'Y' with coordinates for the genes
#' @param geneSets geneSets list of gene sets, e.g. gene pathways, output from the 'getGeneSets' function in MLP
#' the genes IDs must correspond to the sampleNames in the eset
#' @param geneSetsVar variable of the featureData used to match the genes contained in geneSets,
#' most probably ENTREZID, if NULL the featureNames of the eSet are used
#' @param ... Any parameters passed to the \code{\link{getTopElements}}
#' @param esetUsed expressionSet (or SummarizedExperiment) object with data
#' @return data.frame with two columns 'X' and 'Y' with coordinates for the gene sets
#' @author Laure Cougnaud
#' @import Biobase
#' @keywords internal
#' @author Laure Cougnaud
getCoordGeneSets <- function(dataPlotGenes, geneSets, esetUsed, geneSetsVar = list()){
	
	# get methods depending on the class of the object
	esetMethods <- getMethodsInputObjectEsetVis(esetUsed)
	
	#previously only, but very slow
#system.time(test <- as.data.frame(t(sapply(geneSets, function(x) 
#		colMeans(dataPlotGenes[matchGeneID(x), c("X", "Y")], na.rm = TRUE)))))
	
	if(any(duplicated(names(geneSets))))
		warning("Some gene sets have the same name, they will be combined ",
				"for the extraction of the top gene sets.")
	
	# convert to vector
	system.time(geneSetsVect <- unlist(geneSets, recursive = FALSE, use.names = TRUE)) # 1s
	names(geneSetsVect) <- rep(names(geneSets), times = sapply(geneSets, length))
	
	# match with gene ID
	useVarIDToMatch <- if(length(geneSetsVar) > 0)	
		geneSetsVar %in% esetMethods$fvarLabels(esetUsed)	else	FALSE
	system.time(geneSetsVectInEset <- esetMethods$featureNames(esetUsed)[
		match(geneSetsVect, 
			if(useVarIDToMatch)	esetMethods$fData(esetUsed)[, geneSetsVar]	else
				esetMethods$featureNames(esetUsed)
		)
	])
	names(geneSetsVectInEset) <- names(geneSetsVect)
	
	# filtered the ones not in eset
	geneSetsVectInEsetFiltered <- geneSetsVectInEset[!is.na(geneSetsVectInEset)]
	
	# extract coordinates
	system.time(dataPlotGenesGS <- dataPlotGenes[geneSetsVectInEsetFiltered, c("X", "Y")]) # 150s
	
	# means by coordinates
	getCoordGeneSets <- function(colName)
		tapply(dataPlotGenesGS[, colName], names(geneSetsVectInEsetFiltered), mean)
	
	system.time(coordGeneSets <- cbind(X = getCoordGeneSets("X"), Y = getCoordGeneSets("Y")))
	
	return(coordGeneSets)
}

