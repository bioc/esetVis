#' @title get gene sets for plot of \linkS4class{eSet} object.
#' @description get and format gene sets to be used as \code{geneSets} for the functions:
#' \code{\link{esetSpectralMap}}, \code{\link{esetLda}}, or \code{\link{esetPlotWrapper}}
#' Use the \code{\link[MLP]{getGeneSets}} function to get the gene sets,
#' combine all databases, and format the gene sets name if required.
#' @param entrezIdentifiers string with Entrez Gene identifiers of the genes of interest
#' @param species species to use, given to the \code{\link[MLP]{getGeneSets}} function
#' @param geneSetSource gene set source, either 'GOBP', 'GOMF', 'GOCC' or 'KEGG'.
#' Multiple choices are available
#' @param useDescription logical, if TRUE (by default) use the description to 
#' label the gene sets, otherwise use the original gene set identifiers
#' Function 'substr' is used.
#' @param trace logical, if TRUE (by default) 
#' a few extra information are printed during the process
#' @seealso the function used internally: \code{\link[MLP]{getGeneSets}}
#' @return list with gene sets, each element is a gene set and 
#' contains the ENTREZ IDs of the genes contained in this set.
#' If \code{useDescription} is:
#' \itemize{
#' 	\item{FALSE: }{pathways are labelled with identifiers (Gene Ontology IDs for GOBP, GOMF and GOCC, KEGG IDs for KEGG)}
#'  \item{TRUE: }{pathways are labelled with gene sets descriptions}
#' }
#' @examples 
#' # example dataset
#' library(ALL)
#' data(ALL)
#' 
#' # get gene annotation from probe IDs
#' library("hgu95av2.db")
#' probeIDs <- featureNames(ALL)
#' geneInfo <- select(hgu95av2.db, probeIDs,"ENTREZID", "PROBEID")
#' 
#' # get pathway annotation for the genes contained in the ALL dataset (can take a few minutes)
#' geneSets <- getGeneSetsForPlot(entrezIdentifiers = geneInfo$ENTREZID, species = "Human", 
#' 	geneSetSource = 'GOBP',
#' 	useDescription = FALSE, trace = TRUE)
#' head(geneSets) # returns a pathway list of genes
#' 
#' # gene sets labelled with gene sets description
#' geneSets <- getGeneSetsForPlot(entrezIdentifiers = geneInfo$ENTREZID, species = "Human", 
#' 	geneSetSource = 'GOBP', useDescription = TRUE, trace = TRUE)
#' head(geneSets) # returns a pathway list of genes
#' 
#' # see also vignette for an example of the use of this function as input for the esetSpectralMap, esetLda or esetPlotWrapper functions
#' @importFrom MLP getGeneSets
#' @author Laure Cougnaud
#' @export
getGeneSetsForPlot <- function(entrezIdentifiers, species = "Human", 
	geneSetSource = c('GOBP', 'GOMF', 'GOCC', 'KEGG'),
	useDescription = TRUE, trace = TRUE){
	
	geneSetSource <- match.arg(geneSetSource, 
		choices = c('GOBP', 'GOMF', 'GOCC', 'KEGG'), several.ok = TRUE)

	# function from the MLP package, 36 s for 4 databases
	timeGetGeneSets <- system.time(geneSetsList <- lapply(geneSetSource, function(db) 
		getGeneSets(species = species, geneSetSource = db,
			entrezIdentifiers = entrezIdentifiers)))["elapsed"]
	names(geneSetsList) <- geneSetSource
	
	if(trace) message(paste("The extraction of the gene sets for the", length(geneSetSource),
		"databases and", length(entrezIdentifiers), "features", "took:", round(timeGetGeneSets,2), "s."))
	
	if(trace){
		nPathwaysPerDb <- sapply(geneSetsList, length)
		message(paste0("The number of available pathways per database is: ",
			paste(names(nPathwaysPerDb), ": ", nPathwaysPerDb, collapse = ", ", sep = "")
		), ".")
	}
	
	geneSets <- unlist(geneSetsList, recursive = FALSE, use.names = FALSE)
	names(geneSets) <- unlist(sapply(geneSetsList, names))
	
	if(useDescription){
	
		# sometimes more description than names?
		geneSetsDescriptionsList <-  sapply(geneSetsList,
			function(x) attributes(x)$descriptions[attributes(x)$names])
		#sapply(geneSetsDescriptions, length)
		names(geneSets) <- unlist(geneSetsDescriptionsList) 
		# not at this level, otherwise, some pathways can have same names
		#substr(unlist(geneSetsDescriptionsList), 1, nMaxChar)
		
		
	}
	
	return(geneSets)
	
}
