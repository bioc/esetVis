#' wrapper to extract useful functions, depending if the object is 
#' an ExpressionSet or a SummarizedExperiment.
#' 
#' This returns an error is \code{x} is not of the correct class.
#' The package \code{SummarizedExperiment} should be available if \code{x} is of class \code{SummarizedExperiment}.
#' @param x object
#' @return if the object is an ExpressionSet or a SummarizedExperiment, 
#' returns a list with the functions specific of the class of \code{x},
#' and equivalent of the ExpressionSet functions: 'sampleNames', 'featureNames',
#' 'fData', 'pData', 'exprs'
#' \itemize{
#' 	 \item{sampleNames: }{sample names}
#'   \item{featureNames: }{feature names}
#'   \item{fData: }{feature annotation}
#'   \item{pData: }{sample annotation}
#'   \item{exprs: }{data matrix}
#'   \item{varLabels: }{sample annotation variables}
#'   \item{fvarLabels: }{feature annotation variables}
#' }
#' @author Laure Cougnaud
#' @import Biobase
#' @keywords internal
getMethodsInputObjectEsetVis <- function(x){
	
	if(inherits(x, what = c('ExpressionSet'))){
		
		return(
			list(
				sampleNames = sampleNames,
				featureNames = featureNames,
				fData = fData, 
				pData = pData,
				exprs = exprs,
				varLabels = varLabels,
				fvarLabels = fvarLabels
			)
		)
		
	}else if(inherits(x, what = c('SummarizedExperiment', 'RangedSummarizedExperiment'))){
		
		if(!requireNamespace("SummarizedExperiment", quietly = TRUE))
			stop(paste("The package 'SummarizedExperiment' need to be loaded to use",
				"the functionalities of esetVis in an 'SummarizedExperiment' object."))
		
		return(
			list(
				sampleNames = colnames,
				featureNames = rownames,
				fData = SummarizedExperiment::rowData, 
				pData = SummarizedExperiment::colData,
				exprs = SummarizedExperiment::assay,
				varLabels = function(x) colnames(SummarizedExperiment::colData(x)),
				fvarLabels = function(x) colnames(SummarizedExperiment::rowData(x))
			)
		)
		
	}else stop("Input object is not of the correct class: 'ExpressionSet' or 'SummarizedExperiment'.")
	
}

