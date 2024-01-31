# simulate open robust design data

# simulate data for ORD 
## count data = aerial survey
## MR data = open robust design

source("scripts/02_simulate-ord-function.R")

# simulation parameters ----
n.reps = 10
n.years = 5
n.sec = 3
n.marked = 20

sim = as.tibble(expand.grid(rep = c(1:n.reps),
                            n.years = n.years,
                            n.marked = n.marked))
sim %>% 
  group_by(rep) %>% 
  mutate(scenario = c(1:length(n.years))) %>% 
  ungroup() -> sim

# other input parameters constant across scenarios ----

# annual survival
min.phi = 0.6
max.phi = 0.9
sig.phi = plogis(0.01)

# temporary emigration
min.gammaII = 0.6
max.gammaII = 0.9
sig.gammaII = plogis(0.01)

min.gammaOI = 0.1
max.gammaOI = 0.5
sig.gammaOI = plogis(0.01)

# detection
min.p = 0.3
max.p = 0.7
sig.p = plogis(0.01)

# stopover persistence 
min.psi = 0.2
max.psi = 0.9
sig.psi = plogis(0.01)

# within season entry probs
dprior = c(1, 1, 2, 3, 3, 2, 1, 1, 1)

# transience
min.tau = 0.2
max.tau = 0.5
sig.tau = plogis(0.01)

fixed.parms = list(min.phi = min.phi, max.phi = max.phi, sig.phi = sig.phi,
                   min.gammaII = min.gammaII, max.gammaII = max.gammaII,
                   sig.gammaII = sig.gammaII, min.gammaOI = min.gammaOI, 
                   max.gammaOI = max.gammaOI, sig.gammaOI = sig.gammaOI,
                   min.p = min.p, max.p = max.p, sig.p = sig.p,
                   min.psi = min.psi, max.psi = max.psi, sig.psi = sig.psi,
                   min.tau = min.tau, max.tau = max.tau, sig.tau = sig.tau,
                   n.sec = n.sec, dprior = dprior)


# run simulation ----

# simulate data
sim %>%
  mutate(fixed.parms = list(fixed.parms),
         dat = pmap(list(n.years, n.marked, fixed.parms),
                    simulate_ord)) %>%
  select(-fixed.parms) -> sim
