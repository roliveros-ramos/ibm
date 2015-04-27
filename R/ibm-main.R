# ibm package: Individual Based Models in R -------------------------------

#' Individual Based Models in R
#' 
#' Individual Based Models in R 
#' 
#' \tabular{ll}{ Package: \tab ibm\cr Type: \tab Package\cr Version: \tab
#' 0.1\cr Date: \tab 2015-04-27\cr License: \tab GPL-2\cr } ibm()
#' 
#' @name ibm-package
#' @aliases ibm-package ibm
#' @docType package
#' @author Ricardo Oliveros-Ramos Maintainer: Ricardo Oliveros-Ramos
#' <ricardo.oliveros@@gmail.com>
#' @references ibm
#' @keywords ibm
#' @examples
#' 
#' LV.ibm()
NULL

LV.ibm = function(pop, natality, mortality, sd, K, dim, dd=TRUE, boundaries="toro") {
  N = length(pop)/dim
  births = 0
  deaths = 0
  if(N!=0) 
  {
    
    if(length(natality)==1) natality = rep(natality,len=N)
    if(length(mortality)==1) mortality = rep(mortality,len=N)
    if(length(natality)!=N) stop("natality rates and population size don't match") 
    if(length(mortality)!=N) stop("mortality rates and population size don't match") 
    boundaries = match.fun(boundaries)
    brownian = match.fun("brownian")
    
    natality.N   = sapply(X=natality,FUN=muestreo)
    natN 		= which(natality.N==1)
    births		= sum(natality.N)
    
    mortality.N 	= sapply(X=mortality,FUN=muestreo)
    morN 		= which(mortality.N==1)
    deaths		= sum(mortality.N)
    
    if(!all(length(pop)>=K,dd)) pop = rbind(pop,pop[natN,,drop=FALSE])
    if(length(morN)!=0) pop = pop[-morN,,drop=FALSE]
    
    pop = brownian(pos=pop,sd=sd)
    pop = boundaries(x=pop)	
    
  }
    return(list(pop=pop,births=births,deaths=deaths))
  }




