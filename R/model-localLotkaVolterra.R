
# Main function -----------------------------------------------------------


.PredatorPreyModel = function(par, T) {
  if(!requireNamespace("deSolve", quietly = TRUE)) 
    stop("You need to install the 'deSolve' package.")# check on hadley
  # par is a list with 'alpha', 'beta' 'gamma', 'sd' and 'mu_ini'.
  LV = function(t, y, parms, ...) {
    r = parms$r
    l = parms$l
    alpha = parms$alpha
    gamma = parms$gamma
    K = parms$K
    dN = r*y[1]*(1-(y[1]/K)) - alpha*y[1]*y[2]
    dP = -l*y[2] + gamma*alpha*y[1]*y[2]
    return(list(c(dN, dP)))
  }
  times = seq(0, T)
  y0 = c(par$initial$N, par$initial$P)
  sol = deSolve::ode(y=y0, times=times, func=LV, parms=par, method="ode45")
  out = as.list(as.data.frame(sol[,-1]))
  names(out) = c("prey", "predator")
  out$prey[is.na(out$prey)] = 0
  out$predator[is.na(out$predator)] = 0
  return(out)
}



.localLotkaVolterra = function(par, dim=1, periodic=TRUE, plot=FALSE, delay=0.1) {
  
  par = list(initial=list(N=100, P=10), 
             D=list(N, P) )
  
  K 		= constants[1]
  T 		= constants[2]

  r0 	= parameters[1,1] # maximal natality rate - prey
  m0 	= parameters[2,1] # natural mortality rate - predator
  
  alpha 	= 10^(-parameters[1,2])
  beta 	= 10^(-parameters[2,2])
  
  DN 	= parameters[1,3] # diffusion coefficient prey
  DP 	= parameters[2,3] # diffusion coefficiente predator 
  
  R1 	= parameters[1,4] # prey
  R2 	= parameters[2,4] # predator
  
  N.initial = as.integer(parameters[1,5])
  P.initial = as.integer(parameters[2,5])
  
  sdN = sqrt(2*DN)
  sdP = sqrt(2*DP)
  
  brownian = match.fun(paste("brownian.", dim,"D",sep=""))
  
  boundaries = if(periodic) "toro" else "reflejo"
  
  # Initializing population vectors
  
  N = numeric(T)
  P = numeric(T)
  
  N[1] = N.initial
  P[1] = P.initial
  
  # Initializing individuals positions
  pop.N.ini = pop.initial(N=N.initial,n=n)
  pop.P.ini = pop.initial(N=P.initial,n=n)
  
  pop.N = pop.N.ini
  pop.P = pop.P.ini
  
  if(isTRUE(plot)) plot.pop(N=N,P=P,pop.N=pop.N,pop.P=pop.P,dim=n,delay=delay)
  
  dem = matrix(nrow=T, ncol=4)
  
  for(t in 2:T) {
    pred = local.int(posx=pop.N,posy=pop.P,R1=R1,R2=R2)
    NP = pred$x 		# number of predators near to each prey
    PN = pred$y		# number of preys near to each predator
    
    r = rates(r0) 		# natality prey
    s = rates(alpha*NP)	# mortality prey
    l = rates(beta*PN) 	# natality predator
    m = rates(m0) 		# mortality predator
    
    sim.N = LV.ibm(pop=pop.N, natality=r, mortality=s, sd=sdN, K=K, dim=n, dd=TRUE)		# prey dynamics
    sim.P = LV.ibm(pop=pop.P, natality=l, mortality=m, sd=sdP, K=K, dim=n, dd=FALSE)		# predator dynamics
    
    pop.N = sim.N$pop
    pop.P = sim.P$pop 
    
    N[t] = length(pop.N)/n  	
    P[t] = length(pop.P)/n  	
  
    diffusion = function(x, sd, simetric=TRUE) {
      out = brownian(x, sd=sd)
      if(isTRUE(simetric)) out = toro(out)
      return(out)
    } 
    
    pop = brownian(pos=pop, sd=sd)
    pop = toro(x=pop)	
    
    if(isTRUE(plot)) plot.pop(N=N,P=P,pop.N=pop.N,pop.P=pop.P,dim=n,delay=delay)
    
  }
  
  ######## OUTPUTS ############
  
  output = list(#
    prey 		= N, 		# 1
    predator 	= P, 		# 2
    # end of calibration variables
    pop.N		= pop.N,	# 3
    pop.P		= pop.P	# 4
  )
  
  return(output)
  
}


# Generate simulated data -------------------------------------------------

.generatePredatorPreyModel = function(path, r=0.5, l=0.2, alpha=0.1, gamma=0.1, K=100, T=100, 
                                      N0=10, P0=1, ...) {
  
  # 'real' parameters
  par_real = list(r=r, l=l, K=K, alpha=alpha, gamma=gamma, initial=list(N=N0, P=P0))
  
  pop = .PredatorPreyModel(par=par_real, T=T)
  
  # observed abundances
  n = rapply(pop, f=jitter, how = "list") 
  
  main.folder   = file.path(path, "PredatorPreyDemo")
  data.folder   = file.path(main.folder, "data")
  
  if(!file.exists(data.folder)) dir.create(data.folder, recursive=TRUE)
  
  for(i in c("prey", "predator")) {
    ifile = paste0(i, ".csv")
    dat = matrix(n[[i]], ncol=1)
    colnames(dat) = i
    write.csv(dat, file.path(data.folder, ifile))
  }
  
  # parInfo.csv
  
  parInfo = list()
  parInfo$guess = list(r=0.1, l=0.1, K=1.1*max(n$prey), alpha=0.05, gamma=0.1, initial=list(N=n$prey[1], P=n$predator[1]))
  parInfo$lower = list(r=0, l=0, K=0.25*max(n$prey), alpha=0, gamma=0, initial=list(N=0.5*n$prey[1], P=0.5*n$predator[1]))
  parInfo$upper = list(r=2, l=2, K=5*max(n$prey), alpha=1, gamma=1, initial=list(N=1.5*n$prey[1], P=1.5*n$predator[1]))
  parInfo$phase = list(r=1, l=1, K=1, alpha=1, gamma=1, initial=list(N=NA, P=NA))
  
  # calibrationInfo.csv
  
  calibrationInfo = list()
  calibrationInfo$variable  = c("prey", "predator")
  calibrationInfo$type      = "lnorm2"
  calibrationInfo$calibrate = TRUE
  calibrationInfo$weights   = 1
  calibrationInfo$useData   = TRUE
  
  calibrationInfo = as.data.frame(calibrationInfo)
  
  write.csv(calibrationInfo, file.path(main.folder, "calibrationInfo.csv"), row.names=FALSE)
  
  constants = list(T=T)
  
  output = c(list(path=main.folder, par=par_real), constants, parInfo)
  
  return(output)
  
}

# Auxiliar functions ------------------------------------------------------



