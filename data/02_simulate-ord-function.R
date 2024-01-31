# simulate data under open robust design

simulate_ord = function(n.years, n.marked, fixed.parms){
 
  n.sec = fixed.parms$n.sec
  
  ### hyperparameters ###
  mean.phi = runif(1, fixed.parms$min.phi, fixed.parms$max.phi)
  mean.gammaII = runif(1, fixed.parms$min.gammaII, fixed.parms$max.gammaII)
  mean.gammaOI = runif(1, fixed.parms$min.gammaOI, fixed.parms$max.gammaOI)
  mean.tau = runif(1, fixed.parms$min.tau, fixed.parms$max.tau)
  mean.p = runif(1, fixed.parms$min.p, fixed.parms$max.p)
  mean.psi = runif(n.sec, fixed.parms$min.psi, fixed.parms$max.psi)

  ### convert to logit scale ###
  mu.phi <- qlogis(mean.phi)
  mu.gammaII <- qlogis(mean.gammaII)
  mu.gammaOI <- qlogis(mean.gammaOI)
  mu.tau <- qlogis(mean.tau)
  mu.p = qlogis(mean.p)
  mu.psi = qlogis(mean.psi)

  
  ### draw random effects ###
  
  # year-specific parameters
  eps.phi = rnorm(n.years, 0, fixed.parms$sig.phi)
  phi = plogis(mu.phi + eps.phi)
  
  eps.gammaII = rnorm(n.years, 0, fixed.parms$sig.gammaII)
  gammaII = plogis(mu.gammaII + eps.gammaII)
  gammaII[1] = gammaII[n.years] = mean.gammaII
  
  eps.gammaOI = rnorm(n.years, 0, fixed.parms$sig.gammaOI)
  gammaOI = plogis(mu.gammaOI + eps.gammaOI)
  gammaOI[1] = gammaOI[2] = gammaOI[n.years] = mean.gammaOI
  
  eps.tau = rnorm(n.years, 0, fixed.parms$sig.tau)
  tau = plogis(mu.tau + eps.tau)
  
  
  # year- and period-specific parameters
  eps.psi = matrix(rnorm(n.years*n.sec, 0, fixed.parms$sig.psi),
                   nrow = n.years, ncol = n.sec)
  psi = matrix(plogis(mu.psi + eps.psi),
               nrow = n.years, ncol = n.sec, byrow = T)

  eps.p = matrix(rnorm(n.years*n.sec, 0, fixed.parms$sig.p),
                 nrow = n.years, ncol = n.sec)
  p = matrix(plogis(mu.p + eps.p),
             nrow = n.years, ncol = n.sec, byrow = T)
  
  delta = ldelta = matrix(nrow = n.years, ncol = n.sec)
  for(j in 1:n.years){
    for(t in 1:n.sec){
      ldelta[j,t] = rgamma(1, fixed.parms$dprior[t], 1)
    }
  }
  for(j in 1:n.years){
    delta[j,] = ldelta[j,]/sum(ldelta[j,1:n.sec])
  }
  # delta[,n.sec] <- 0

  
  # probability of non-reencounter
  q = 1-p
  
  # persistence probability after first entry
  psi.t = apply(psi, 2, FUN = function(x) x*tau)
  
  
  # probability of being available during each secondary occasion
  # alpha = prob present
  # omega = prob present AND not detected
  
  o = k = alpha = omega = z = pi = matrix(nrow = n.years, ncol = n.sec)
  for(i in 1:n.years){
    o[i,1] <- delta[i,1]

    k[i,1] <- 0
    k[i,2] <- delta[i,1]*psi[i,1] + delta[i,2]

    for(j in 2:n.sec){
      k[i,j] <- k[i,j-1]*psi[i,j-1] + delta[i,j]
      o[i,j] <- o[i,j-1]*psi[i,j-1]*(1-p[i,j-1]) + delta[i,j]
    }

    alpha[i,1] <- delta[i,1]
    omega[i,1] <- delta[i,1]
    for(j in 2:n.sec){
      alpha[i,j] <- (tau[i]*(k[i,j]-delta[i,j])) + delta[i,j]
      omega[i,j] <- (tau[i]*(o[i,j]-delta[i,j])) + delta[i,j]
    }

    for(j in 1:n.sec){
      # with temp emi
      z[i,j] <- gammaII[i]*alpha[i,j] + gammaOI[i]*alpha[i,j]
      # make sure pi < 1
      pi[i,j] <- min(1, z[i,j])

    }
  }
  
  pstar = apply(omega*p, 1, sum)
  qstar = 1-pstar
  
  
  ### capture-recapture model ###
  
  # these matrices/vectors simplify (kind of) the m-array setup
  # open robust design formulation from Kendall and Bjorkland 2001
  # G = prob of being unavailable (1) or available but not seen (2) given state last year
  # a = prob of being available if unavailable (1) or available (2) last year
  # f = prob of being unavailable (1) or available but not seen (2) if available last year
  # use j+1 indexing for f and G because e.g. column 1 refers to year 2
  
  arr_prod = function(x){
    nmat <- dim(x)[2]
    y <- x[,1,]
    
    for(i in 2:nmat){
      y <- y %*% x[,i,]
    }
    return(y)
  }
  
  G = array(NA, dim = c(2, n.years, 2))
  f = matrix(NA, nrow = n.years, ncol = 2)
  a = matrix(NA, nrow = 2, ncol = n.years)
  U = matrix(NA, nrow = n.years-1, ncol = n.years)
  
  for(j in 1:n.years){
    
    G[1,j,1] <- 1-gammaOI[j]
    G[1,j,2] <- gammaOI[j]*qstar[j]
    G[2,j,1] <- 1-gammaII[j]
    G[2,j,2] <- gammaII[j]*qstar[j]
    
    f[j,1] <- 1-gammaII[j]
    f[j,2] <- gammaII[j]*qstar[j]
    
    a[1,j] <- gammaOI[j]
    a[2,j] <- gammaII[j]
  }

  # primary period m-array
  # need custom nimbleFunction arr_prod for matrix multiplication within arrays
  
  prL1 = matrix(NA, nrow = n.years-1, ncol = n.years)
  
  # main diagonal
  for(j in 1:(n.years-1)){
    prL1[j,j] <- phi[j] * gammaII[j+1] * pstar[j+1]
  }
  
  # above main diagonal
  for(j in 1:(n.years-3)){
    prL1[j,j+1] <- prod(phi[j:(j+1)]) * f[j+1,1:2] %*% a[1:2,j+2] * pstar[j+2]
    prL1[j,j+2] <- prod(phi[j:(j+2)]) * f[j+1,1:2] %*% G[1:2,j+2,1:2] %*% a[1:2,j+3] * pstar[j+3]
  }
  
  # upper triangle
  for(j in 1:(n.years-4)){
    for(t in (j+3):(n.years-1)){
      prL1[j,t] <- prod(phi[j:t]) * f[j+1,1:2] %*% arr_prod(G[1:2,(j+2):t,1:2]) %*% a[1:2,t+1] * pstar[t+1]
    }
  }
  prL1[(n.years-2), (n.years-1)] <- prod(phi[(n.years-2):(n.years-1)]) * f[(n.years-1),1:2] %*% a[1:2,n.years] * pstar[n.years]
  
  # below main diagonal
  for(j in 2:(n.years-1)){
    for(k in 1:(j-1)){
      prL1[j,k] <- 0
    }
  }
  
  # last column - probability of non-recapture
  for(j in 1:(n.years-1)){
    prL1[j,n.years] <- 1-sum(prL1[j,1:(n.years-1)])
  }
  
  
  
  # seconday period m-arrays
  prL2a = matrix(NA, nrow = n.years, ncol = n.sec)
  prL2c.r = prL2c.t = array(NA, dim = c(n.years, n.sec-1, n.sec))
  
  for(j in 1:n.years){
    
    # first encounter within primary period
    for(t in 1:n.sec){
      prL2a[j,t] <- omega[j,t]*p[j,t]
    }
    
    for(t in 1:(n.sec-1)){
      
      prL2c.t[j,t,t] <- psi.t[j,t]*p[j,t]
      prL2c.r[j,t,t] <- psi[j,t]*p[j,t]
      
      for(k in (t+1):(n.sec-1)){
        prL2c.t[j,t,k] <- psi.t[j,t]*prod(psi[j,(t+1):k]) * prod(q[j,t:(k-1)]) * p[j,k]
        prL2c.r[j,t,k] <- prod(psi[j,t:k]) * prod(q[j,t:(k-1)]) * p[j,k]
      }
      
      for(k in 1:(t-1)){
        prL2c.t[j,t,k] <- 0
        prL2c.r[j,t,k] <- 0
      }
    }
    
    # last column - never re-encountered
    for(t in 1:(n.sec-1)){
      prL2c.t[j,t,n.sec] <- 1-sum(prL2c.t[j,t,1:(n.sec-1)])
      prL2c.r[j,t,n.sec] <- 1-sum(prL2c.r[j,t,1:(n.sec-1)])
    }
  }
  
  rel = rep(n.marked, n.years-1)
  marr.p = matrix(NA, nrow = n.years-1, ncol = n.years)
  new.s = matrix(NA, nrow = n.years, ncol = n.sec)
  
  # across years
  for(j in 1:(n.years-1)){
    
    # L1 - encounters across primary periods
    rel[j] = ifelse(j>1, sum(marr.p[,j-1], rel[j], na.rm = T), rel[j])
    marr.p[j,1:n.years] = rmultinom(1, rel[j], prL1[j,])
    
  }
  
  rel[n.years] = sum(marr.p[,n.years])
  
  for(j in 1:n.years){
    #L2a - new encounters each secondary period
    new.s[j,1:n.sec] = rmultinom(1, rel[j], prL2a[j,])
  }
  
  # within season re-encounters
  all.marr.r = all.marr.t = array(NA, dim = c(n.years, n.sec-1, n.sec))
  rel.s = matrix(NA, nrow = n.years, ncol = n.sec-1)
  
  for(j in 1:n.years){
    for(t in 1:(n.sec-1)){
      # first resight of new individuals (possible transients)
      all.marr.t[j,t,1:n.sec] = rmultinom(1, new.s[j,t], prL2c.t[j,t,])
      
      # resighting inds previously seen that year
      rel.s[j,t] = ifelse(t > 1, sum(all.marr.r[j,,t-1], all.marr.t[j,,t-1],
                                     na.rm = T), 0)
      all.marr.r[j,t,1:n.sec] = rmultinom(1, rel.s[j,t], prL2c.r[j,t,])
    }
  }
  
  
  
  ### outputs ###
  jagsdat = list(marr.p = marr.p,
                 new.s = new.s,
                 marr.s.r = all.marr.r,
                 marr.s.t = all.marr.t,
                 rel = rel,
                 rel.sec.r = rel.s,
                 rel.sec.t = new.s,
                 n.years = n.years,
                 n.sec = n.sec)
  
  oparms <- tibble(parm = c("mean.phi", "mean.gammaII", "mean.gammaOI", "mean.tau", 
                            "sig.phi", "sig.gammaII", "sig.gammaOI", "sig.tau", 
                            "sig.psi", "sig.p", "mean.p"),
                   truth = c(mean.phi, mean.gammaII, mean.gammaOI, mean.tau, 
                             sig.phi, sig.gammaII, sig.gammaOI, sig.tau, 
                             sig.psi, sig.p, mean.p))
  
  yparms <- as.tibble(
    expand.grid(year = c(1:n.years),
                parm = c("phi", "gammaII", "gammaOI", "tau", "pstar")))
  yparms$truth = c(phi, gammaII, gammaOI, tau, pstar)
  
  occparms <- as.tibble(
    expand.grid(occ = c(1:n.sec),
                parm = c("mean.psi")))
  occparms$truth = c(mean.psi)
  
  yoparms <- as.tibble(
    expand.grid(year = c(1:n.years),
                occ = c(1:n.sec),
                parm = c("psi", "p", "delta"))
  )
  yoparms$truth = c(c(psi), c(p), c(delta))
  
  true.params <- oparms %>% 
    full_join(yparms) %>% 
    full_join(occparms) %>% 
    full_join(yoparms)
  
  list(jagsdat = jagsdat,
       true.params = true.params)
  
}



