# Local integration -------------------------------------------------------

local.int = function(x, y, R1, R2=R1) {
    #    if(ncol(x)!=ncol(y)) stop("Dimensions of x and y must agree.")
    #    if(ncol(x)>3) stop("Dimensions must be less than 3.")
    Nx = nrow(x) 
    Px = nrow(y)
    
    if(all(Nx!=0,Px!=0)) { 
      
      int    = local.int.cpp(posx=x, posy=y, R1=R1, R2=R2)
      intX   = int$x
      intY   = int$y
      
    } else {
      
      intX = rep(0, length=Nx)
      intY = rep(0, length=Px)	
      
    }
    
    return(list(x=intX, y=intY))
  }


# Movement ----------------------------------------------------------------

toro = function(x) {
  x = ifelse(x<0, x+1, x)
  x = ifelse(x>1, x-1, x)
  return(x)
  }


# Brownian motion ---------------------------------------------------------


brownian1D = function(object, sd) {
  out 	= object + rnorm(length(object), mean=0, sd=sd)
  return(out)
}

brownian2D = function(object, sd) {
  N	= nrow(object)
  angulo	= runif(N)
  dist	= rnorm(N, mean=0, sd=sd)
  out 	= object + dist*cbind(cos(pi*angulo),sin(pi*angulo))
  return(out)
}

brownian3D = function(object, sd) {
  N	= nrow(object)
  angulo1	= runif(N)
  angulo2	= runif(N)
  dist	= rnorm(N, mean=0, sd=sd)
  out 	= object + dist*cbind(sin(pi*angulo2)*cos(pi*angulo1), 
                             sin(pi*angulo2)*sin(pi*angulo1), 
                             cos(pi*angulo2))
  return(out)
}


# sampling ----------------------------------------------------------------

.checkRates = function(rates, n) {
  if(length(rates)==1) rates = rep(rates, length=n)
  if(length(rates)!=n) stop("Rates and population size do not match.") 
  rates = pmin(pmax(rates, 0),1)
  if(any(is.na(rates))) warning("NA rates are being taken as zero.")
  rates[is.na(rates)] = 0
  return(rates)
}



