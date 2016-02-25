
# Brownian motion ---------------------------------------------------------

#' @title Brownian diffusion of a set of particles
#' @description This funtions performs a brownian difussion over a set of particles.
#' The dimension is automatically calculated from the number of columns of the object.
#' @param object The positions of the particles, dimension is taken from the number
#' of columns or assumed to be 1 is no columns.
#' @param sd Standard deviation for the gaussian jump, for dynamics models should be
#' set proportional to \code{sqrt(dt)}.
#' @details This functions apply a brownian diffusion to a set of point coordinates.
#' @export
brownian = function(object, sd) {
  UseMethod("brownian")
}

#' @export
brownian.default = brownian1D

#' @export
brownian.matrix = function(object, sd) {
  n = ncol(x)
  if(n<1 | n>3) stop("Dimension for brownian diffusion must be 1, 2 or 3.")
  out = switch (n,
                "1" = brownian1D(object=object, sd=sd),
                "2" = brownian2D(object=object, sd=sd),
                "3" = brownian3D(object=object, sd=sd)
  )
  return(out)
}

#' @export
brownian.array = function(object, sd) {
  
  out  = apply(X=object, MARGIN=seq_along(dim(object))[-c(1,2)], FUN=brownian, sd=sd)
  dim(out) = dim(object)
  return(out)
  
}
