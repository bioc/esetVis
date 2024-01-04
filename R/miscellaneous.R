#' check if the aesthetic is fixed (e.g. color, shape, size 'palette')
#' @param typeVar name of variable for aesthetic
#' @param valVar fixed value of variable of aesthetic
#' @return logical, if TRUE the element is fixed
#' @author Laure Cougnaud
#' @keywords internal
setFixElement <- function(typeVar, valVar)	(length(typeVar) == 0) & (length(valVar) > 0)

#' check if manual aesthetic should be set
#' 
#' This is the case only if \code{typeVar} and \code{valVar} are specified,
#' and if the variable is not numeric or integer (doesn't work with ggplot2)
#' @param x data.frame with \code{typeVar}
#' @param typeVar name of variable for aesthetic
#' @param valVar fixed value of variable of aesthetic
#' @return logical, if TRUE the manual scale should be set
#' @author Laure Cougnaud
#' @keywords internal
setManualScale <- function(x, typeVar, valVar)	
	(length(typeVar)) > 0 & (length(valVar) > 0) & !class(x[, typeVar]) %in% c("numeric", "integer")

#' extend manual scale values if required
#' @param x data.frame with \code{nameVar}
#' @param valVar fixed value of variable of aesthetic
#' @param nameVar name of variable for aesthetic
#' @return vector of manual scales
#' @author Laure Cougnaud
#' @keywords internal
formatManualScale <- function(x, valVar, nameVar){
	values <- rep(valVar, length.out = nlevels(factor(x[, nameVar])))
	names(values) <- NULL #cannot provide named argument for colors
	values
}

#' capitalize the first letter of a word
#' @param x string
#' @return string with first letter capitalized
#' @keywords internal
simpleCap <- function(x) {
	
	paste0c <- function(...) paste(..., sep = "", collapse = " ")
	paste0c(toupper(substring(x, 1, 1)), substring(x, 2))
	
}

#' format output of \link{plotEset} function
#' @param res result of specific \link{plotEset} function
#' @param object \code{esetPlot} object or extended class
#' @param type string type of plot
#' @param returnEsetPlot logical, should the object be returned in the output function?
#' @return result
#' @author Laure Cougnaud
#' @keywords internal
formatOutput <- function(res, object, type, returnEsetPlot){
	if(returnEsetPlot)
		res <- c(
			if(inherits(res, type))	list(plot = res)	else	res,
			list(esetPlot = object)
		)
	return(res)
}	