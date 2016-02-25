# Reproduction process ----------------------------------------------------


#' @title Reproduction Process
#' @description This functions performs the 'reproduction' process over an object,
#' increasing the number of individuals. It is a generic, S3 methods can be specified
#' for a particular specification of the population. 
#' @param object The population object, containing the information about individuals.
#' @param rate The reproduction rate or rates.
#' @details The rate can be a single value or a value for each individual calculated
#' externally. No recycling is allowed.
reproduction = function(object, rates, ...) {
  UseMethod("reproduction")
}

#' @export
#' @method reproduction default
reproduction.default = function(object, rates, random=TRUE, newborns=1) {
  n = if(is.matrix(object)) nrow(object) else length(object)
  pop = seq_len(n)
  rates = .checkRates(rates, n)
  if(isTRUE(random)) {
    pop = sample(pop, replace=TRUE)
    newborns = 1
  }
  new = rbinom(n = n, size = newborns, prob=rates)
  new = rep(pop, times=new)
  return(new)
}
