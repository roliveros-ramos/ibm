# Local integration -------------------------------------------------------

local.int = function(posx,posy,R1,R2=R1) {
    #    if(ncol(posx)!=ncol(posy)) stop("Dimensions of posx and posy must agree.")
    #    if(ncol(posx)>3) stop("Dimensions must be less than 3.")
    Nx = nrow(posx) 
    Px = nrow(posy)
    if(all(Nx!=0,Px!=0)) { 
      
      int     = local.int.cpp(posx=posx, posy=posy, R1=R1, R2=R2)
      int.x   = int$x
      int.y   = int$y
      
    } else {
      
      int.x = rep(0,len=Nx)
      int.y = rep(0,len=Px)	
      
    }
    
    return(list(x=int.x,y=int.y))
  }


# Movement ----------------------------------------------------------------

toro = function(x) {
    x[x<0]=x[x<0]+1
    x[x>1]=x[x>1]-1
    return(x)
  }



# Brownian motion ---------------------------------------------------------

brownian.1D = function(pos,sd) {
    out 	= pos + rnorm(length(pos), mean=0, sd=sd)
    return(out)
  }

brownian.2D = function(pos,sd) {
    N	= length(pos)/2
    angulo	= runif(N)
    dist	= rnorm(N, mean=0, sd=sd)
    out 	= pos + dist*cbind(cos(pi*angulo),sin(pi*angulo))
    return(out)
  }

brownian.3D = function(pos,sd) {
    N	= length(pos)/3
    angulo1	= runif(N)
    angulo2	= runif(N)
    dist	= rnorm(N, mean=0, sd=sd)
    out 	= pos + dist*cbind(sin(pi*angulo2)*cos(pi*angulo1), sin(pi*angulo2)*sin(pi*angulo1), cos(pi*angulo2))
    return(out)
  }


# sampling ----------------------------------------------------------------

muestreo =  function(prob) {
    x = sample(c(1,0), 1, prob=c(prob,1-prob))
    return(x)
  }

rates = function(x) {
    x[x<0]=0
    x[x>1]=1
    return(x)
  }


# Initialization ----------------------------------------------------------

pop.initial = function(N,n,min=0,max=1) {
    out = matrix(runif(n*N,min=min,max=max),ncol=n,nrow=N)
    return(out)
  }



# Visualization -----------------------------------------------------------

plot.pop = function(N,P,pop.N,pop.P,dim,delay=0,prey.col="blue",pred.col="red") {
    layout(matrix(c(1,2),nrow=2), heights=c(3,1))
    par(mar=c(3,3,1,1),oma=c(1,1,1,1))
    if(dim==1)    {
      plot(density(pop.N,from=0,to=1),col=prey.col)
      lines(density(pop.P,from=0,to=1),col=pred.col)
      plot(N,type="l",ylim=c(0,max(N,P,na.rm=TRUE)),col=prey.col)
      lines(P,col=pred.col)
    }
    if(dim==2)    {
      plot(pop.N[,1],pop.N[,2],xlim=c(0,1),ylim=c(0,1),col=prey.col,pch=19)
      points(pop.P[,1],pop.P[,2],col=pred.col,pch=19)
      plot(N,type="l",ylim=c(0,max(N,P,na.rm=TRUE)),col=prey.col)
      lines(P,col=pred.col)
    }
    Sys.sleep(delay)
}






