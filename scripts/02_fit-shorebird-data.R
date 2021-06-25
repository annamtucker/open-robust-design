# fit open robust design data to shorebird data
# set working directory to main file location (~/shorebird-ord)

# load packages and helper functions
source("scripts/00_setup.R")

# read in data
# in the data files:
# marr.p = primary period m-array
# new.s = number of newly-encountered birds in each period (columns) of each year (rows)
# marr.s.r = secondary period m-arrays for known residents
# marr.s.t = secondary period m-arrays for possible transients
# rel = number "released" each year (resightings and physical captures)
# rel.sec.r = number "released" in each secondary period, known residents
# rel.sec.t = number "released" in each secondary period, possible transients
# tot = total number of flagged birds encountered in each year 

# rekn = red knot
# rutu = ruddy turnstone
# sand = sanderling

rekn_dat <- readRDS("data/rekn_dat_ord.rds")
rutu_dat <- readRDS("data/rutu_dat_ord.rds")
sand_dat <- readRDS("data/sand_dat_ord.rds")

# build and compile model 

n.years = 14
n.sec = 9

constants = list(n.years = n.years,
                 n.sec = n.sec,
                 dpriors = c(1, 1, 2, 3, 3, 2, 1, 1, 1))

inits = function(){list(mean.phi = runif(1, 0.4, 0.9),
                        mean.gammaII = runif(1, 0.1, 0.9),
                        mean.gammaOI = runif(1, 0.1, 0.9),
                        mean.tau = runif(1, 0.1, 0.9),
                        mean.p = runif(1, 0.1, 0.9),
                        sig.phi = runif(1, 1, 5),
                        sig.gammaII = runif(1, 1, 5),
                        sig.gammaOI = runif(1, 1, 5),
                        sig.tau = runif(1, 1, 5),
                        sig.p = runif(1, 1, 5),
                        eps.tau = rnorm(n.years, 0, 0.05),
                        eps.gammaOI = rnorm(n.years, 0, 0.05), 
                        eps.gammaII = rnorm(n.years*2, 0, 0.05),
                        eps.phi = rnorm((n.years-1), 0, 0.05),
                        eps.p = matrix(rnorm(n.years*n.sec, 0, 0.05),
                                       nrow = n.years, ncol = n.sec),
                        psi = matrix(runif(n.years*(n.sec-1), 0.1, 0.9),
                                     nrow = n.years, ncol = n.sec-1),
                        delta = matrix(rep(rdirch(1, alpha = c(1, 1, 2, 3, 3, 2, 1, 1, 1)), n.years),
                                       nrow = n.years, ncol = n.sec, byrow = T)
)}


# model code
source("scripts/01_nimble-model.R")

# initialize with rekn data
ordMod <- nimbleModel(ordCodeBeta, 
                      constants = constants,
                      data = rekn_dat,
                      inits = inits())

# check that all nodes are initialized
ordMod$initializeInfo()

# check that logProbs are not NA
ordMod$logProb_marr.p
ordMod$logProb_new.s
ordMod$logProb_marr.s.r
ordMod$logProb_marr.s.t

# compile model
Cord <- compileNimble(ordMod, 
                      showCompilerOutput = T,
                      resetFunctions = T)

# configure MCMC
ordMCMCconfig <- configureMCMC(ordMod)

# nodes to monitor 
parms = c("phi", "p", "delta", "psi", "tau", "gammaII", "gammaOI", "pstar", "omega", "pi",
          "mean.phi", "sig.phi", "mean.tau", "sig.tau", "mean.gammaII", "sig.gammaII",
          "mean.gammaOI", "sig.gammaOI", "var.phi", "var.gammaII", "var.gammaOI", "var.tau")

ordMCMCconfig$resetMonitors()
ordMCMCconfig$addMonitors(parms)

# to add all stochastic nodes (needed for WAIC)
# snodes = ordMod$getNodeNames(stochOnly = T, includeData = F)
# ordMCMCconfig$addMonitors(c(snodes, parms))

# build MCMC
ordMCMC <- buildMCMC(ordMCMCconfig, enableWAIC = F)

# compile MCMC
CordMCMC <- compileNimble(ordMCMC,
                          project = ordMod,
                          showCompilerOutput = T,
                          resetFunctions = T)

# sample MCMC (red knot)
rekn <- runMCMC(CordMCMC, 
                inits = inits,
                niter = 500000,
                nburnin = 250000,
                thin = 25,
                nchains = 3,
                samplesAsCodaMCMC = TRUE)

# check convergence with Rhat and traceplots
rhat = Rhat(rekn)
hist(rhat)

par(mfrow = c(5,5))
traceplot(rekn)

# save model results
saveRDS(rekn, "output/rekn-ord.rds")



# refit compiled model with ruddy turnstone data

Cord$setData(rutu_dat)

rutu <- runMCMC(CordMCMC, 
                inits = inits,
                niter = 250000,
                nburnin = 150000,
                thin = 10,
                nchains = 3,
                samplesAsCodaMCMC = TRUE)

rhat = Rhat(rutu)
hist(rhat)

par(mfrow = c(5,5))
traceplot(rutu)

saveRDS(rutu, "output/rutu-ord.rds")


# refit compiled model with sanderling data

Cord$setData(sand_dat)

sand <- runMCMC(CordMCMC, 
                inits = inits,
                niter = 250000,
                nburnin = 150000,
                thin = 10,
                nchains = 3,
                samplesAsCodaMCMC = TRUE)

rhat = Rhat(sand)
hist(rhat)

par(mfrow = c(5,5))
traceplot(sand)

saveRDS(sand, "output/sand-ord.rds")
