# Mortality process -------------------------------------------------------

#' @title Mortality Process
#' @description This functions performs the 'mortality' process over an object,
#' decreasing the number of individuals. It is a generic, S3 methods can be specified
#' for a particular specification of the population. 
#' @param object The population object, containing the information about individuals.
#' @param rate The mortality rate or rates.
#' @details The rate can be a single value or a value for each individual calculated
#' externally. No recycling is allowed.
mortality = function(object, rates, ...) {
  UseMethod("mortality")  
}


#' @export
#' @method mortality default
mortality.default = function(object, rates) {
  n = if(is.matrix(object)) nrow(object) else length(object)
  pop = seq_len(n)
  rates = .checkRates(rates, n)
  die = which(rbinom(n = n, size = 1, prob=rates)==1)
  return(die)
}

