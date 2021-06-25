# summarize open robust design model results

# load packages and functions
source("scripts/00_setup.R")

# read in nimble output
rekn = readRDS("output/rekn-ord.rds")
rutu = readRDS("output/rutu-ord.rds")
sand = readRDS("output/sand-ord.rds")

# convert to dataframes and combine into one to make figures
rekn.res <- mcmc_to_df(rekn, parms)
rutu.res <- mcmc_to_df(rutu, parms)
sand.res <- mcmc_to_df(sand, parms)

rekn.res %>% 
  mutate(species = "rekn") -> rekn.res

rutu.res %>% 
  mutate(species = "rutu") -> rutu.res

res <- sand.res %>% 
  mutate(species = "sand") %>% 
  rbind(rekn.res) %>% 
  rbind(rutu.res) %>% 
  mutate(year = c(2005:2018)[year])


# compare priors and posteriors 
priors = tibble(x = seq(0, 1, 0.01),
                phi = dbeta(seq(0, 1, 0.01), 5, 2),
                gammaII = dbeta(seq(0, 1, 0.01), 3, 3),
                gammaOI = dbeta(seq(0, 1, 0.01), 3, 3),
                tau = dbeta(seq(0, 1, 0.01), 1, 1),
                p = dbeta(seq(0, 1, 0.01), 3, 3))

res %>% 
  filter(parm == "mean.phi") %>% 
  ggplot() +
  geom_density(aes(x = val), lwd = 1, col = "steelblue4") +
  facet_wrap(~species) +
  xlim(0, 1) +
  geom_line(data = priors, aes(x = x, y = phi), lwd = 1)

res %>% 
  filter(parm == "mean.gammaII") %>% 
  ggplot() +  geom_density(aes(x = val), lwd = 1, col = "steelblue4") +
  facet_wrap(~species) +
  xlim(0, 1) +
  geom_line(data = priors, aes(x = x, y = gammaII), lwd = 1)

res %>% 
  filter(parm == "mean.gammaOI") %>% 
  ggplot() +
  geom_density(aes(x = val), lwd = 1, col = "steelblue4") +
  facet_wrap(~species) +
  xlim(0, 1) +
  geom_line(data = priors, aes(x = x, y = gammaOI), lwd = 1)

res %>% 
  filter(parm == "mean.tau") %>% 
  ggplot() +
  geom_density(aes(x = val), lwd = 1, col = "steelblue4") +
  facet_wrap(~species) +
  xlim(0, 1) +
  geom_line(data = priors, aes(x = x, y = tau), lwd = 1)


res_sum <- res %>% 
  group_by(species, year, occ, parm) %>% 
  mean_qi(est = val) %>% 
  ungroup()

# global averages
res_global <- res_sum %>% 
  filter(str_detect(parm, "mean")) %>% 
  mutate(parm = substr(parm, 6, nchar(parm))) %>% 
  mutate(species = c("rekn" = "Red knot",
                     "rutu" = "Ruddy turnstone",
                     "sand" = "Sanderling")[species])


# across-year parameters

res %>% 
  filter(parm == "phi") %>% 
  ggplot(aes(x = as.character(year), y = val)) +
  geom_violin() +
  facet_wrap(~species) +
  ggtitle("survival")

res %>% 
  filter(parm == "gammaII") %>% 
  ggplot(aes(x = as.character(year), y = val)) +
  geom_violin() +
  facet_wrap(~species) +
  ggtitle("temp emi in-in")

res %>% 
  filter(parm == "gammaOI") %>% 
  ggplot(aes(x = as.character(year), y = val)) +
  geom_violin() +
  facet_wrap(~species) +
  ggtitle("temp emi out-in")

res %>% 
  filter(parm == "tau") %>% 
  ggplot(aes(x = as.character(year), y = val)) +
  geom_violin() +
  facet_wrap(~species) +
  ggtitle("residency")


# within-year parameters

ggplot(res, aes(x = occ, group = occ, y = delta)) +
  geom_violin() +
  facet_grid(species~year) +
  ggtitle("arrival")

ggplot(res, aes(x = occ, group = occ, y = psi)) +
  geom_violin() +
  facet_grid(species~year) +
  ggtitle("persistence")


ggplot(res, aes(x = occ, group = occ, y = p)) +
  geom_violin() +
  facet_grid(species~year) +
  ggtitle("detection")


# check that posterior dists are nearly symmetrical (for Pop Eco AE comment)
reknL1 <- res %>% 
  filter(species == "rekn") %>% 
  select(-year, -occ) %>% 
  filter(parm %in% c("mean.phi", "var.phi", "mean.gammaII", "var.gammaII",
                     "mean.gammaOI", "var.gammaOI", "mean.tau", "var.tau")) %>% 
  spread(parm, val) %>% 
  mutate(sig.phi = sqrt(var.phi/(mean.phi^2 * (1-mean.phi)^2)),
         sig.gammaII = sqrt(var.gammaII/(mean.gammaII^2 * (1-mean.gammaII)^2)),
         sig.gammaOI = sqrt(var.gammaOI/(mean.gammaOI^2 * (1-mean.gammaOI)^2)),
         sig.tau = sqrt(var.tau/(mean.tau^2 * (1-mean.tau)^2))) %>% 
  gather(parm, val, 3:14) %>% 
  filter(!parm %in% c("var.phi", "var.gammaII", "var.gammaOI", "var.tau")) %>% 
  ggplot(aes(x = val)) +
  geom_density(fill = "gray80") +
  facet_wrap(~parm, scales = "free") +
  ggtitle("Red knot - L1 hyperparameters")
reknL1

res %>% 
  filter(species == "rekn") %>% 
  filter(parm %in% c("delta")) %>% 
  ggplot(aes(x = val)) +
  geom_density(fill = "gray80") +
  facet_grid(occ~year, scales = "free") +
  ggtitle("Red knot - delta")

res %>% 
  filter(species == "rekn") %>% 
  filter(parm %in% c("psi")) %>% 
  ggplot(aes(x = val)) +
  geom_density(fill = "gray80") +
  facet_grid(occ~year, scales = "free") +
  ggtitle("Red knot - psi")


res %>% 
  filter(species == "rekn") %>% 
  filter(parm %in% c("p")) %>% 
  ggplot(aes(x = val)) +
  geom_density(fill = "gray80") +
  facet_grid(occ~year, scales = "free") +
  ggtitle("Red knot - p")


rutuL1 <- res %>% 
  filter(species == "rutu") %>% 
  select(-year, -occ) %>% 
  filter(parm %in% c("mean.phi", "var.phi", "mean.gammaII", "var.gammaII",
                     "mean.gammaOI", "var.gammaOI", "mean.tau", "var.tau")) %>% 
  spread(parm, val) %>% 
  mutate(sig.phi = sqrt(var.phi/(mean.phi^2 * (1-mean.phi)^2)),
         sig.gammaII = sqrt(var.gammaII/(mean.gammaII^2 * (1-mean.gammaII)^2)),
         sig.gammaOI = sqrt(var.gammaOI/(mean.gammaOI^2 * (1-mean.gammaOI)^2)),
         sig.tau = sqrt(var.tau/(mean.tau^2 * (1-mean.tau)^2))) %>% 
  gather(parm, val, 3:14) %>% 
  filter(!parm %in% c("var.phi", "var.gammaII", "var.gammaOI", "var.tau")) %>% 
  ggplot(aes(x = val)) +
  geom_density(fill = "gray80") +
  facet_wrap(~parm, scales = "free") +
  ggtitle("Ruddy turnstone - L1 hyperparameters")
rutuL1

res %>% 
  filter(species == "rutu") %>% 
  filter(parm %in% c("delta")) %>% 
  ggplot(aes(x = val)) +
  geom_density(fill = "gray80") +
  facet_grid(occ~year, scales = "free") +
  ggtitle("Ruddy turnstone - delta")

res %>% 
  filter(species == "rutu") %>% 
  filter(parm %in% c("psi")) %>% 
  ggplot(aes(x = val)) +
  geom_density(fill = "gray80") +
  facet_grid(occ~year, scales = "free") +
  ggtitle("Ruddy turnstone - psi")


res %>% 
  filter(species == "rutu") %>% 
  filter(parm %in% c("p")) %>% 
  ggplot(aes(x = val)) +
  geom_density(fill = "gray80") +
  facet_grid(occ~year, scales = "free") +
  ggtitle("Ruddy turnstone - p")


sandL1 <- res %>% 
  filter(species == "sand") %>% 
  select(-year, -occ) %>% 
  filter(parm %in% c("mean.phi", "var.phi", "mean.gammaII", "var.gammaII",
                     "mean.gammaOI", "var.gammaOI", "mean.tau", "var.tau")) %>% 
  spread(parm, val) %>% 
  mutate(sig.phi = sqrt(var.phi/(mean.phi^2 * (1-mean.phi)^2)),
         sig.gammaII = sqrt(var.gammaII/(mean.gammaII^2 * (1-mean.gammaII)^2)),
         sig.gammaOI = sqrt(var.gammaOI/(mean.gammaOI^2 * (1-mean.gammaOI)^2)),
         sig.tau = sqrt(var.tau/(mean.tau^2 * (1-mean.tau)^2))) %>% 
  gather(parm, val, 3:14) %>% 
  filter(!parm %in% c("var.phi", "var.gammaII", "var.gammaOI", "var.tau")) %>% 
  ggplot(aes(x = val)) +
  geom_density(fill = "gray80") +
  facet_wrap(~parm, scales = "free") +
  ggtitle("Sanderling - L1 hyperparameters")
sandL1

res %>% 
  filter(species == "sand") %>% 
  filter(parm %in% c("delta")) %>% 
  ggplot(aes(x = val)) +
  geom_density(fill = "gray80") +
  facet_grid(occ~year, scales = "free") +
  ggtitle("Sanderling - delta")

res %>% 
  filter(species == "sand") %>% 
  filter(parm %in% c("psi")) %>% 
  ggplot(aes(x = val)) +
  geom_density(fill = "gray80") +
  facet_grid(occ~year, scales = "free") +
  ggtitle("Sanderling - psi")


res %>% 
  filter(species == "sand") %>% 
  filter(parm %in% c("p")) %>% 
  ggplot(aes(x = val)) +
  geom_density(fill = "gray80") +
  facet_grid(occ~year, scales = "free") +
  ggtitle("Ruddy turnstone - p")



# global averages ----
res_global <- res_sum %>% 
  filter(str_detect(parm, "mean")) %>% 
  mutate(parm = substr(parm, 6, nchar(parm))) %>% 
  mutate(species = c("rekn" = "Red knot",
                     "rutu" = "Ruddy turnstone",
                     "sand" = "Sanderling")[species])


res_global_mean <- res_sum %>% 
  filter(str_detect(parm, "mean")) %>% 
  mutate(parm = substr(parm, 6, nchar(parm))) %>% 
  mutate(species = c("rekn" = "Red knot",
                     "rutu" = "Ruddy turnstone",
                     "sand" = "Sanderling")[species]) %>% 
  distinct(species, parm, est, .lower, .upper) %>% 
  mutate(type = "mean")
res_global_mean

res_global_var <- res_sum %>% 
  filter(str_detect(parm, "var")) %>% 
  mutate(parm = substr(parm, 5, nchar(parm))) %>% 
  mutate(species = c("rekn" = "Red knot",
                     "rutu" = "Ruddy turnstone",
                     "sand" = "Sanderling")[species]) %>% 
  distinct(species, parm, est, .lower, .upper) %>% 
  mutate(type = "var")
res_global_var


res_global_var %>% arrange(-est)

# table 1 - hyperparameters ----
res_global_mean %>% 
  gather(var, val, 3:5) %>% 
  mutate(val = round(val, 2)) %>% 
  spread(var, val) %>% 
  mutate(lower = gsub(" ", "", paste(.lower,",")),
         mean.cri = paste(lower, .upper)) %>% 
  rename(mean = est) %>% 
  select(species, parm, mean, mean.cri) ->gmean

res_global_var %>% 
  gather(var, val, 3:5) %>% 
  mutate(val = round(val, 3)) %>% 
  spread(var, val) %>% 
  mutate(lower = gsub(" ", "", paste(.lower,",")),
         var.cri = paste(lower, .upper)) %>% 
  rename(var = est) %>% 
  select(species, parm, var, var.cri) %>% 
  full_join(gmean) %>% 
  select(species, parm, mean, mean.cri, var, var.cri) %>% 
  write.csv("tables/table3_hyperparms.csv", row.names = F)


# table s1 - annual estimates ----
res_sum %>% 
  filter(parm %in% c("phi", "gammaII", "gammaOI", "tau")) %>% 
  gather(var, val, 5:7) %>% 
  mutate(val = round(val, 2)) %>% 
  spread(var, val) %>% 
  mutate(lower = gsub(" ", "", paste("(", .lower,",")),
         upper = gsub(" ", "", paste(.upper, ")")),
         var.cri = paste(lower, upper),
         val = paste(est, var.cri)) %>% 
  select(species, year, parm, val) %>% 
  arrange(species, year) %>% 
  spread(parm, val) %>% 
  write.csv("tables/tableS1_annual-ests.csv", row.names = F)

# figure 2 - temp emi, residency, survival ----
# temporary emigration

res_sum %>% 
  filter(parm == "gammaOI") %>%
  filter(! year %in% c(2005, 2006, 2018)) -> gammaOI

tempemi <- res_sum %>% 
  filter(parm == "gammaII") %>%
  filter(year > 2005 & year < 2018) %>% 
  full_join(gammaOI) %>% 
  mutate(species = c("rekn" = "Red knot",
                     "rutu" = "Ruddy turnstone",
                     "sand" = "Sanderling")[species]) %>%
  ggplot(aes(shape = parm)) +
  geom_rect(data = res_global[res_global$parm == "gammaOI",],
            aes(ymin = .lower, ymax = .upper), xmin = 2005, xmax = 2018,
            alpha = 0.3, fill = "gray60") +
  geom_hline(data = res_global[res_global$parm == "gammaOI",],
             aes(yintercept = est), lty = 2, lwd = 1) +
  geom_rect(data = res_global[res_global$parm == "gammaII",],
            aes(ymin = .lower, ymax = .upper), xmin = 2005, xmax = 2018,
            alpha = 0.3, fill = "gray60") +
  geom_hline(data = res_global[res_global$parm == "gammaII",],
             aes(yintercept = est), lty = 2, lwd = 1) +
  geom_linerange(aes(x = year, ymin = .lower, ymax = .upper), lwd = 1) +
  geom_point(aes(x = year, y = est), size = 3) +
  facet_wrap(~species, scales = "free") +
  ylim(0, 1) +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(margin = margin(1,1,2,1),
                                  size = 16),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14)) +
  scale_x_continuous(breaks = seq(2006, 2018, 3), limits = c(2006, 2018)) +
  scale_shape_manual(labels = c(bquote(gamma^{II}), bquote(gamma^{OI})),
                       values = c(16, 15),
                       name = "") +
  xlab("") +
  ylab("Temporary emigration\nprobability")
tempemi  


# residency 
tau <- res_sum %>% 
  filter(parm == "tau") %>% 
  filter(year > 2005) %>% 
  mutate(species = c("rekn" = "Red knot",
                     "rutu" = "Ruddy turnstone",
                     "sand" = "Sanderling")[species]) %>% 
  ggplot(aes(x = year)) +
  geom_rect(data = res_global[res_global$parm == "tau",],
              aes(ymin = .lower, ymax = .upper), xmin = 2005, xmax = 2018,
              alpha = 0.3, fill = "gray60") +
  geom_hline(data = res_global[res_global$parm == "tau",],
            aes(yintercept = est), lty = 2, lwd = 1) +
  geom_linerange(aes(ymin = .lower, ymax = .upper), lwd = 1) +
  geom_point(aes(y = est), size = 3) +
  facet_wrap(~species, scales = "free") +
  ylim(0.4, 1) +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(margin = margin(1,1,2,1),
                                  size = 16),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) +
  scale_x_continuous(breaks = seq(2005, 2017, 3)) +
  xlab("") +
  ylab("Stopover residency\nprobability (\U1D70F)")
tau


phi <- res_sum %>% 
  filter(parm == "phi") %>% 
  mutate(species = c("rekn" = "Red knot",
                     "rutu" = "Ruddy turnstone",
                     "sand" = "Sanderling")[species]) %>% 
  ggplot(aes(x = year)) +
  geom_rect(data = res_global[res_global$parm == "phi",],
            aes(ymin = .lower, ymax = .upper), xmin = 2005, xmax = 2018,
            alpha = 0.3, fill = "gray60") +
  geom_hline(data = res_global[res_global$parm == "phi",],
             aes(yintercept = est), lty = 2, lwd = 1) +
  geom_linerange(aes(ymin = .lower, ymax = .upper), lwd = 1) +
  geom_point(aes(y = est), size = 3) +
  facet_wrap(~species, scales = "free") +
  ylim(0.65, 1) +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(margin = margin(1,1,2,1),
                                  size = 16),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) +
  scale_x_continuous(breaks = seq(2005, 2017, 3)) +
  xlab("Year") +
  ylab("Apparent annual\nsurvival probability (\U03D5)")
phi

fig2 <- plot_grid(tempemi, tau, phi, align = "v", axis = "lrtb",
                  ncol = 1, labels = c("(a)", "(b)", "(c)"))
fig2

# save_plot(fig2, file = "figs/fig2_temp-emi_residency_survival.tiff",
#           base_width = 10, base_height = 10)


 # arrival ---- 

res_sum %>% 
  filter(parm == "delta") %>% 
  mutate(species = c("rekn" = "Red knot",
                     "rutu" = "Ruddy turnstone",
                     "sand" = "Sanderling")[species]) %>% 
  ggplot(aes(x = occ, col = as.character(year), fill = as.character(year))) +
  # geom_ribbon(aes(ymin = .lower, ymax = .upper),
  #             alpha = 0.5) +
  # geom_path(aes(y = est), lty = 1) +
  geom_linerange(aes(ymin = .lower, ymax = .upper)) +
  geom_point(aes(y = est)) +
  facet_grid(~species, scales = "free") +
  ylim(0, 1) +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(margin = margin(1,1,2,1))) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  scale_x_continuous(breaks = seq(1, 9, 2)) +
  xlab("Sampling occasion") +
  ylab("Arrival probability")


# persistence ---- 

res_sum %>% 
  filter(parm == "psi" & occ < 9) %>% 
  mutate(species = c("rekn" = "Red knot",
                     "rutu" = "Ruddy turnstone",
                     "sand" = "Sanderling")[species]) %>% 
  ggplot(aes(x = occ)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              alpha = 0.5, fill = "gray60") +
  geom_path(aes(y = est), lty = 1) +
  facet_grid(species~year, scales = "free") +
  ylim(0, 1) +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(margin = margin(1,1,2,1))) +
  scale_x_continuous(breaks = seq(1, 9, 2)) +
  xlab("Sampling occasion") +
  ylab("Persistence probability")


# proportion present ---- 

# alpha  = proportion of pop that stopped there that year
# pi = proportion of total flyway pop

calc_pi = function(x){
  
  delta <- x$delta 
  
  psi <- x$psi
  
  tau <- unique(x$tau)[2]
  
  gammaII <- unique(x$gammaII)[2]
  
  gammaOI <- unique(x$gammaOI)[2]
  
  # k = probability of arrived in previous period and still there
  # g = probabiilty departed in that period
  
  k = alpha = numeric()
  
  k[1] <- delta[1]
  for(t in 2:length(delta)){
    k[t] <- k[t-1]*psi[t-1] + delta[t]
  }
  
  alpha[1] <- delta[1]
  for(t in 2:length(delta)){
    alpha[t] <- (tau*(k[t]-delta[t])) + delta[t]
  }
  
  pi <- gammaII*alpha + gammaOI*alpha
  
  return(pi)
  
}

pi <- res %>% 
  filter(parm %in% c("delta", "psi", "tau", "gammaII", "gammaOI")) %>% 
  spread(parm, val) %>% 
  group_by(species, year, draw) %>% 
  nest() %>% 
  mutate(pi = map(data, calc_pi)) %>% 
  unnest() %>% 
  ungroup() %>% 
  group_by(species, year, occ) %>% 
  mean_qi(est = pi)



figS2 <- pi %>% 
  ungroup() %>% 
  filter(year > 2005) %>% 
  mutate(species = c("rekn" = "Red knot",
                     "rutu" = "Ruddy turnstone",
                     "sand" = "Sanderling")[species]) %>% 
  ggplot(aes(x = occ, col = as.character(year), 
             fill = as.character(year))) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              alpha = 0.4) +
  geom_path(aes(y = est), lwd = 1) +
  ylim(0, 1) +
  facet_wrap(~species, scales = "free") +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(margin = margin(1,1,2,1),
                                  size = 14),
        legend.position = "top") +
  scale_color_viridis_d(name = "Year", direction = -1) +
  scale_fill_viridis_d(name = "Year", direction = -1) +
  scale_x_continuous(breaks = seq(2, 9, 2)) +
  xlab("Sampling occasion") +
  ylab("Proportion of the population present")

figS2


# save_plot(figS2, file = "figs/figS2_prop-present.jpg",
#           base_width = 10, base_height = 6)



# fig 4 ----
fig4A <- pi %>% 
  ungroup() %>% 
  filter(year > 2005) %>% 
  mutate(species = c("rekn" = "Red knot",
                     "rutu" = "Ruddy turnstone",
                     "sand" = "Sanderling")[species],
         occ =  6 + (occ*3)) %>% 
  ggplot(aes(x = occ, col = as.character(year), 
             fill = as.character(year))) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              alpha = 0.1, col = NA) +
  geom_line(aes(y = est), lwd = 1, alpha = 0.8) +
  geom_vline(xintercept = 26, lty = 2, lwd = 1)+
  geom_vline(xintercept = 28, lty = 2, lwd = 1)+
  #geom_path(aes(y = est), lwd = 1) +
  ylim(0, 1) +
  facet_wrap(~species, scales = "free") +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 14, margin = margin(1,1,2,1)),
        legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  scale_color_viridis_d(name = "Year", direction = -1) +
  scale_fill_viridis_d(name = "Year", direction = -1) +
  #scale_x_continuous(breaks = seq(2, 9, 2)) +
  #guides(fill = guide_legend(ncol = 7),
  #       col = guide_legend(ncol = 7)) +
  xlab("Day in May") +
  ylab("Probability of presence at\nstopover site (\U1D70B)")
fig4A

fig4B <- pi %>% 
  ungroup() %>% 
  filter(occ == 7 & year > 2005) %>% 
  ggplot(aes(x = year, y = est, col = as.character(year))) +
  geom_linerange(aes(ymin = .lower, ymax = .upper), lwd = 1.5) +
  geom_point(size = 3) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  facet_wrap(~species, scales = "free") +
  ylim(0, 1) +
  scale_x_continuous(breaks = seq(2006, 2018, 3)) +
  scale_color_viridis_d(name = "Year", direction = -1) +
  xlab("Year") +
  ylab("Probability of presence\nMay 26-28") +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(margin = margin(1,1,2,1))) 

legend <- get_legend(
  fig4A + 
    theme(legend.position = "top") +
    guides(fill = guide_legend(ncol = 7),
           col = guide_legend(ncol = 7))
  )

fig4 <- plot_grid(fig4A, fig4B, legend, ncol = 1,
                  rel_heights = c(1, 1, 0.2),
                  labels = c("(a)", "(b)", ""),
                  label_size = 16,
                  scale = 0.95)
fig4

# save_plot(fig4, base_width = 10, base_height = 8, file = "figs/fig4_pi.tiff")


