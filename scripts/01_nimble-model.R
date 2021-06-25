# nimble code for ORD model

# function to perform matrix multiplication across elements of an array
arr_prod <- nimbleFunction(    
  run = function(x = double(3)){
    
    nmat <- dim(x)[2]
    y <- x[,1,]
    
    for(i in 2:nmat){
      y <- y %*% x[,i,]
    }
    
    returnType(double(2))
    return(y)
  })


ordCodeBeta <- nimbleCode({ 
  
  # -- PARAMETERS -- #
  # indexing
  # j = primary period (year)
  # t = secondary period (3-day sampling occ)
  
  # primary periods 
  # gammaOI[j] = probability of being available in year j if unavailable j-1
  # gammaII[j] = probability of being available in year j if available in j-1
  # phi[j] = survival probability year j to j+1
  # pstar[j] = detected at least once in year j if available
  
  # secondary periods
  # psi[j,t] = probability of remaining from t to t+1 in year j
  # delta[j,t] = probability of entering just before occ t of year j
  # p[j,t] = detection probability year j occ t
  # tau[j] = residency probability year j (remaining for at least 2 sampling occs)
  
  
  # -- LIKELIHOOD -- #
  # three parts
  
  for(j in 1:(n.years-1)){
    
    # L1 - encounters across primary periods
    marr.p[j,1:n.years] ~ dmulti(prL1[j,1:n.years], rel[j])
  }
  
  for(j in 1:n.years){  
    #L2 - new encounters each secondary period
    new.s[j,1:n.sec] ~ dmulti(prL2[j,1:n.sec], tot[j])
    
    #L3 - subsequent encounters across secondary periods 
    for(t in 1:(n.sec-1)){
      
      # possible transients 
      marr.s.t[j,t,1:n.sec] ~ dmulti(prL3.t[j,t,1:n.sec], rel.sec.t[j,t])
      
      # known residents
      marr.s.r[j,t,1:n.sec] ~ dmulti(prL3.r[j,t,1:n.sec], rel.sec.r[j,t])
    }
    
  }
  
  # -- GOODNESS OF FIT -- #
  # calculate expected cell frequencies and Freeman-Tukey statistic for real and simulated data

  # for(j in 1:(n.years-1)){
  # 
  #   # L1
  #   marr.new.p[j,1:n.years] ~ dmulti(prL1[j,1:n.years], rel[j])
  # 
  #   for(k in 1:n.years){
  #     expmarr.p[j,k] <- rel[j]*prL1[j,k]
  #     E.org.p[j,k] <- pow((pow(marr.p[j,k], 0.5) - pow(expmarr.p[j,k], 0.5)), 2)
  #     E.new.p[j,k] <- pow((pow(marr.new.p[j,k], 0.5) - pow(expmarr.p[j,k], 0.5)), 2)
  #   }
  # }
  # 
  # for(j in 1:n.years){
  # 
  #   # L2
  #   new.s.new[j,1:n.sec] ~ dmulti(prL2[j,1:n.sec], tot[j])
  #   for(t in 1:n.sec){
  #     exp.s[j,t] <- tot[j] * prL2[j,t]
  #     E.org.new[j,t] <- pow((pow(new.s[j,t], 0.5) - pow(exp.s[j,t], 0.5)), 2)
  #     E.new.new[j,t] <- pow((pow(new.s.new[j,t], 0.5) - pow(exp.s[j,t], 0.5)), 2)
  #   }
  # 
  #   # #L3
  #   for(t in 1:(n.sec-1)){
  #     marr.new.s.r[j,t,1:n.sec] ~ dmulti(prL3.r[j,t,1:n.sec], rel.sec.r[j,t])
  #     marr.new.s.t[j,t,1:n.sec] ~ dmulti(prL3.t[j,t,1:n.sec], rel.sec.t[j,t])
  # 
  #     for(l in 1:n.sec){
  #       expmarr.s.r[j,t,l] <- rel.sec.r[j,t]*prL3.r[j,t,l]
  #       expmarr.s.t[j,t,l] <- rel.sec.t[j,t]*prL3.t[j,t,l]
  # 
  #       E.org.s.r[j,t,l] <- pow((pow(marr.s.r[j,t,l], 0.5) - pow(expmarr.s.r[j,t,l], 0.5)), 2)
  #       E.new.s.r[j,t,l] <- pow((pow(marr.new.s.r[j,t,l], 0.5) - pow(expmarr.s.r[j,t,l], 0.5)), 2)
  # 
  #       E.org.s.t[j,t,l] <- pow((pow(marr.s.t[j,t,l], 0.5) - pow(expmarr.s.t[j,t,l], 0.5)), 2)
  #       E.new.s.t[j,t,l] <- pow((pow(marr.new.s.t[j,t,l], 0.5) - pow(expmarr.s.t[j,t,l], 0.5)), 2)
  #     }
  #   }
  # 
  #   fitL3r.y[j] <- sum(E.org.s.r[j, 1:(n.sec-1), 1:n.sec])
  #   fit.newL3r.y[j] <- sum(E.new.s.r[j, 1:(n.sec-1), 1:n.sec])
  # 
  #   fitL3t.y[j] <- sum(E.org.s.t[j, 1:(n.sec-1), 1:n.sec])
  #   fit.newL3t.y[j] <- sum(E.new.s.t[j, 1:(n.sec-1), 1:n.sec])
  # 
  # }
  # 
  # fitL1 <- sum(E.org.p[1:(n.years-1), 1:n.years])
  # fit.newL1 <- sum(E.new.p[1:(n.years-1), 1:n.years])
  # 
  # fitL2 <- sum(E.org.new[1:n.years, 1:n.sec])
  # fit.newL2 <- sum(E.new.new[1:n.years, 1:n.sec])
  # 
  # fitL3r <- sum(fitL3r.y[1:n.years])
  # fit.newL3r <- sum(fit.newL3r.y[1:n.years])
  # 
  # fitL3t <- sum(fitL3t.y[1:n.years])
  # fit.newL3t <- sum(fit.newL3t.y[1:n.years])


  # -- PRIMARY PERIOD CELL PROBABILITIES -- #
  
  # these matrices/vectors simplify the m-array setup
  # G = prob of being unavailable (1) or available but not seen (2) in year j given state in year j-1
  # a = prob of being available in year j if unavailable (1) or available (2) in year j-1
  # f = prob of being unavailable (1) or available (2) in year j if unavailable in year j-1
  
  for(j in 1:n.years){

    G[1,j,1] <- 1-gammaOI[j]
    G[1,j,2] <- gammaOI[j]*qstar[j]
    G[2,j,1] <- 1-gammaII[j]
    G[2,j,2] <- gammaII[j]*qstar[j]
    
    f[j,1] <- 1-gammaII[j]
    f[j,2] <- gammaII[j]*qstar[j]
  }
  
  for(j in 1:n.years){
    a[1,j] <- gammaOI[j]
    a[2,j] <- gammaII[j]
  }
  
  ### TO RUN IN JAGS INSTEAD OF NIMBLE ###
  ## The following takes the place of the `arr_prod()` function that was used in Nimble.  
  ## replace `arr_prod(G[1:2,(j+2):t,1:2])` with `G.prod[s,,,(j+2),t])` in setup below
  
  # for(j in 1:(n.years-4)){
  #   for(t in (j+3):(n.years-1)){
  #     G.pre.prod[s,1:2,1:2,j+2,t,j+2] <- G[s,,j+2,]
  #     for (k in (j+3):t) {
  #       G.pre.prod[s,1:2,1:2,j+2,t,k] <- G.pre.prod[s,,,j+2,t,k-1] %*% G[s,,k,]
  #     }
  #     G.prod[s,1:2,1:2,j+2,t] <- G.pre.prod[s,,,j+2,t,t]
  #   }
  # }
  
  # prL1 = primary period m-array cell proabilities
  # need custom nimbleFunction arr_prod for matrix multiplication within arrays
  
  # main diagonal
  for(j in 1:(n.years-1)){
    prL1[j,j] <- phi[j] * gammaII[j+1] * pstar[j+1]
  }
  
  # above main diagonal
  for(j in 1:(n.years-3)){
    prL1[j,j+1] <- prod(phi[j:(j+1)]) * pstar[j+2] * (f[j+1,1:2] %*% a[1:2,j+2])[1,1] 
    prL1[j,j+2] <- prod(phi[j:(j+2)]) * pstar[j+3] * (f[j+1,1:2] %*% G[1:2,j+2,1:2] %*% a[1:2,j+3])[1,1]
  }
  
  # upper triangle
  for(j in 1:(n.years-4)){
    for(t in (j+3):(n.years-1)){
      prL1[j,t] <- prod(phi[j:t]) * pstar[t+1] * (f[j+1,1:2] %*% arr_prod(G[1:2,(j+2):t,1:2]) %*% a[1:2,t+1])[1,1]
    }
  }
  prL1[(n.years-2), (n.years-1)] <- prod(phi[(n.years-2):(n.years-1)]) * pstar[n.years] * (f[(n.years-1),1:2] %*% a[1:2,n.years])[1,1]
  
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
  
  
  
  # -- SECONDARY PERIOD CELL PROBABILITIES -- #
  
  # prL3.r = secondary period m-array for known residents - 2nd+subsequent resightings
  # prL3.t = secondary period m-array for possible transients - only first and 2nd (if applicable) resightings
  # prL2 = vector - prob first detected at secondary period t
  
  for(j in 1:n.years){
    
    # first encounter within primary period
    for(t in 1:n.sec){
      prL2[j,t] <- omega[j,t]*p[j,t]
    }
    
    for(t in 1:(n.sec-1)){
      
      prL3.t[j,t,t] <- psi.t[j,t]*p[j,t]
      prL3.r[j,t,t] <- psi[j,t]*p[j,t]
      
      for(k in (t+1):(n.sec-1)){
        prL3.t[j,t,k] <- psi.t[j,t]*prod(psi[j,(t+1):k]) * prod(q[j,t:(k-1)]) * p[j,k]
        prL3.r[j,t,k] <- prod(psi[j,t:k]) * prod(q[j,t:(k-1)]) * p[j,k]
      }
      
      for(k in 1:(t-1)){
        prL3.t[j,t,k] <- 0
        prL3.r[j,t,k] <- 0
      }
    }
    
    # last column - never re-encountered
    for(t in 1:(n.sec-1)){
      prL3.t[j,t,n.sec] <- 1-sum(prL3.t[j,t,1:(n.sec-1)])
      prL3.r[j,t,n.sec] <- 1-sum(prL3.r[j,t,1:(n.sec-1)])
    }
  }

  
  # -- DERIVED PARAMETERS -- #
  
  
  # probability of being available during each secondary occasion
  # alpha = prob present in occ t if available (using the site) in year j
  # omega = prob present and not detected OR just arrived if available in year j
  # z = prob present in occ t of year j
  
  for(j in 1:n.years){
    o[j,1] <- 0
    o[j,2] <- delta[j,1]*psi[j,1]*(1-p[j,1]) + delta[j,2]

    k[j,1] <- 0
    k[j,2] <- delta[j,1]*psi[j,1] + delta[j,2]

    for(t in 3:n.sec){
      k[j,t] <- k[j,t-1]*psi[j,t-1] + delta[j,t]
      o[j,t] <- o[j,t-1]*(psi[j,t-1]*(1-p[j,t-1])) + delta[j,t]
    }

    alpha[j,1] <- delta[j,1]
    omega[j,1] <- delta[j,1]
    for(t in 2:n.sec){
      alpha[j,t] <- (tau[j]*(k[j,t]-delta[j,t])) + delta[j,t]
      omega[j,t] <- (tau[j]*(o[j,t]-delta[j,t])) + delta[j,t]
    }

    for(t in 1:n.sec){
      # with temp emi
      z[j,t] <- gammaII[j]*alpha[j,t] + gammaOI[j]*alpha[j,t]
      # make sure pi < 1
      pi[j,t] <- min(1, z[j,t])

    }
    
    pstar[j] <- sum(omega[j,1:n.sec]*p[j,1:n.sec])
    qstar[j] <- 1-pstar[j]
  }
  
  # -- PRIORS AND CONSTRAINTS -- #
  
  ### apparent annual survival probability ###
  # random effect of year
  
  mean.phi ~ dbeta(5,2)
  mu.phi <- log(mean.phi/(1-mean.phi))
  sig.phi ~ dunif(0,5)
  tau.phi <- pow(sig.phi, -2)
  var.phi <- pow(sig.phi, 2) * pow(mean.phi, 2) * pow((1-mean.phi), 2)
  
  for(j in 1:(n.years-1)){
    eps.phi[j] ~ dnorm(0, tau.phi)
    logit(phi[j]) <- mu.phi + eps.phi[j]
  }
  
  
  
  ### temporary emigration probability ###
  # random effect of year
  mean.gammaII ~ dbeta(3,3)
  mu.gammaII <- log(mean.gammaII/(1-mean.gammaII))
  sig.gammaII ~ dunif(0,10)
  tau.gammaII <- pow(sig.gammaII, -2)
  var.gammaII <- pow(sig.gammaII, 2) * pow(mean.gammaII, 2) * pow((1-mean.gammaII), 2)
  
  mean.gammaOI ~ dbeta(3,3)
  mu.gammaOI <- log(mean.gammaOI/(1-mean.gammaOI))
  sig.gammaOI ~ dunif(0,10)
  tau.gammaOI <- pow(sig.gammaOI, -2)
  var.gammaOI <- pow(sig.gammaOI, 2) * pow(mean.gammaOI, 2) * pow((1-mean.gammaOI), 2)
  
  # gammas not identifiable for every year
  for(j in 2:(n.years-1)){
    eps.gammaII[j] ~ dnorm(0, tau.gammaII)
  }
  eps.gammaII[1] <- 0
  eps.gammaII[n.years] <- eps.gammaII[n.years-1]
  
  for(j in 3:(n.years-1)){
    eps.gammaOI[j] ~ dnorm(0, tau.gammaOI)
  }
  eps.gammaOI[1] <- 0
  eps.gammaOI[2] <- 0
  eps.gammaOI[n.years] <- eps.gammaOI[n.years-1]
  
  for(j in 1:n.years){
    logit(gammaII[j]) <- mu.gammaII + eps.gammaII[j]
    logit(gammaOI[j]) <- mu.gammaOI + eps.gammaOI[j]
  }
  
  
  ### residency probability ###
  # random effect of year
  
  mean.tau ~ dbeta(1,1)
  mu.tau <- log(mean.tau/(1-mean.tau))
  sig.tau ~ dunif(0,10)
  tau.tau <- pow(sig.tau, -2)
  var.tau <- pow(sig.tau, 2) * pow(mean.tau, 2) * pow((1-mean.tau), 2)
  
  for(j in 1:n.years){
    eps.tau[j] ~ dnorm(0, tau.tau)
    logit(tau[j]) <- mu.tau + eps.tau[j]
  }
  
  
  ### detection probability ###
  # random effect of year + occ
  
  mean.p ~ dbeta(3, 3)
  mu.p <- log(mean.p/(1-mean.p))
  sig.p ~ dunif(0, 10)
  tau.p <- pow(sig.p, -2)
  
  for (j in 1:n.years) {
    for (t in 2:(n.sec - 1)) {
      eps.p[j,t] ~ dnorm(0, tau.p)
    }
    
    # constrain p for first and last occasion
    eps.p[j,1] <- eps.p[j,2]
    eps.p[j,n.sec] <- eps.p[j,(n.sec-1)]
    
    for (t in 1:n.sec) {
      logit(p[j,t]) <- mu.p + eps.p[j,t]
      q[j,t] <- 1 - p[j,t]
    }
  }
  
  
  ### arrival probability ###
  # fixed effect of year + occ 
  
  # for(t in 1:(n.sec-1)){
  #   mean.delta[t] ~ dbeta(1,1)
  #   mu.delta[t] <- log(mean.delta[t]/(1-mean.delta[t]))
  # }
  # 
  # sig.delta ~ dunif(0,10)
  # tau.delta <- pow(sig.delta, -2)
  # 
  # for(j in 1:n.years){
  #   for(t in 1:(n.sec-1)){
  #     eps.delta[j,t] ~ dnorm(0, tau.delta)
  #     ldelta[j,t] <- mu.delta[t] + eps.delta[j,t]
  #     # ldelta[j,t] ~ dgamma(1,1)
  #   }
  # }
  # 
  # for(j in 1:n.years){
  #   for(t in 1:(n.sec-1)){
  #     # delta[j,t] <- ldelta[j,t]/sum(ldelta[j,1:(n.sec-1)])
  #     delta[j,t] <- exp(ldelta[j,t])/(1 + sum(exp(ldelta[j,1:(n.sec-1)])))
  #   }
  #   delta[j,n.sec] <- 1-sum(delta[j,1:(n.sec-1)])
  # }
  
  for(j in 1:n.years){
    delta[j,1:n.sec] ~ ddirch(dpriors[1:n.sec])
  }
  
  
  
  ### stopover persistence probability ###
  # fixed effect of year + occ
  
  # for(t in 1:(n.sec-1)){
  #   mean.psi[t] ~ dbeta(1,1)
  #   mu.psi[t] <- log(mean.psi[t]/(1-mean.psi[t]))
  # }
  # 
  # sig.psi ~ dunif(0,5)
  # tau.psi <- pow(sig.psi, -2)

  for (j in 1:n.years) {
    for (t in 1:(n.sec-1)) {
      # eps.psi[j,t] ~ dnorm(0, tau.psi)
      # logit(psi[j,t]) <- mu.psi[t] + eps.psi[j,t]
      psi[j,t] ~ dbeta(1,1)
      psi.t[j, t] <- psi[j, t] * tau[j]
    }
  }
  
  
})
