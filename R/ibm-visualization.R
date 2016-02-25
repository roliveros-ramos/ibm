
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






