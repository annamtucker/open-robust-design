# posterior analysis to estimate covariate relationships

# read in covariate data
# continuous covariates have been scaled and centered
# hsc = relative horseshoe crab egg availability
# storm = extreme storm during stopover?
# rekn_mass_rate = red knot average rate of mass gain
# rekn_temp = number of days above 15C before 95% red knot arrival
# rutu_temp = number of days above 15C before 95% ruddy turnstone arrival

covars <- read_csv("data/covars_scaled.csv")

# functions to fit linear models to scaled parameters
fit_covs_phi_hsc = function(x){
  y = x$phi[1:13]
  lm(y ~  covars$hsc + covars$storm)
}

fit_covs_phi_rate = function(x){
  y = x$phi[1:13]
  lm(y ~ covars$rekn_mass_rate)
}


fit_covs_gammaII = function(x){
  y = x$gammaII[2:14]
  lm(y ~  covars$hsc)
}

fit_covs_tau = function(x){
  y = x$tau[1:13]
  lm(y ~ covars$hsc + covars$rekn_temp)
}



# red knot ----

rekn = readRDS("output/rekn-ord.rds")

rekn_res = mcmc_to_df(rekn, parms = c("phi", "gammaII", "tau"))

# fit linear models to each iteration
rekn_path <- rekn_res %>% 
  select(-occ) %>% 
  group_by(parm, draw) %>% 
  mutate(val_sc = (val-mean(val))/sd(val)) %>% 
  ungroup() %>% 
  select(-val) %>% 
  spread(parm, val_sc) %>% 
  group_by(draw) %>% 
  nest(year, phi, gammaII, tau) %>% 
  mutate(phi_hsc = map(data, fit_covs_phi_hsc),
         phi_rate = map(data, fit_covs_phi_rate),
         gammaII = map(data, fit_covs_gammaII),
         tau = map(data, fit_covs_tau)) 

# summarize results
rk_path_res <- rekn_path %>% 
  select(-data) %>% 
  gather(parm, mod, 2:5) %>% 
  mutate(res = map(mod, tidy)) %>% 
  select(-mod) %>% 
  unnest() %>% 
  ungroup()

# beta posterior distributions
rk_path_res %>% 
  filter(parm == "phi_hsc") %>% 
  ggplot(aes(x = estimate)) +
  geom_density() +
  geom_vline(xintercept = 0) +
  facet_wrap(~term, scales = "free") +
  ggtitle("rekn survival")

rk_path_res %>% 
  filter(parm == "phi_rate") %>% 
  ggplot(aes(x = estimate)) +
  geom_density() +
  geom_vline(xintercept = 0) +
  facet_wrap(~term, scales = "free") +
  ggtitle("rekn survival")

rk_path_res %>% 
  filter(parm == "gammaII") %>% 
  ggplot(aes(x = estimate)) +
  geom_density() +
  geom_vline(xintercept = 0) +
  facet_wrap(~term, scales = "free") +
  ggtitle("rekn temp emi")

rk_path_res %>% 
  filter(parm == "tau") %>% 
  ggplot(aes(x = estimate)) +
  geom_density() +
  geom_vline(xintercept = 0) +
  facet_wrap(~term, scales = "free") +
  ggtitle("rekn residency")


# ruddy turnstone ----

rutu = readRDS("output/rutu-ord.rds")

rutu_res = mcmc_to_df(rutu, parms = c("phi", "gammaII", "tau"))

# fit linear models to each iteration
rutu_path <- rutu_res %>% 
  select(-occ) %>% 
  group_by(parm, draw) %>% 
  mutate(val_sc = (val-mean(val))/sd(val)) %>% 
  ungroup() %>% 
  select(-val) %>% 
  spread(parm, val_sc) %>% 
  group_by(draw) %>% 
  nest(year, phi, gammaII, tau) %>% 
  mutate(phi = map(data, fit_covs_phi_hsc),
         gammaII = map(data, fit_covs_gammaII),
         tau = map(data, fit_covs_tau)) 

# summarize results
rt_path_res <- rutu_path %>% 
  select(-data) %>% 
  gather(parm, mod, 2:4) %>% 
  mutate(res = map(mod, tidy)) %>% 
  select(-mod) %>% 
  unnest() %>% 
  ungroup()

# beta posterior distributions
rt_path_res %>% 
  filter(parm == "phi") %>% 
  ggplot(aes(x = estimate)) +
  geom_density() +
  geom_vline(xintercept = 0) +
  facet_wrap(~term, scales = "free") +
  ggtitle("survival")

rt_path_res %>% 
  filter(parm == "gammaII") %>% 
  ggplot(aes(x = estimate)) +
  geom_density() +
  geom_vline(xintercept = 0) +
  facet_wrap(~term, scales = "free") +
  ggtitle("temp emi")

rt_path_res %>% 
  filter(parm == "tau") %>% 
  ggplot(aes(x = estimate)) +
  geom_density() +
  geom_vline(xintercept = 0) +
  facet_wrap(~term, scales = "free") +
  ggtitle("residency")



# sanderling ----

sand = readRDS("output/sand-ord.rds")

sand_res = mcmc_to_df(sand, parms = c("phi", "gammaII", "tau"))

# fit linear models to each iteration
sand_path <- sand_res %>% 
  select(-occ) %>% 
  group_by(parm, draw) %>% 
  mutate(val_sc = (val-mean(val))/sd(val)) %>% 
  ungroup() %>% 
  select(-val) %>% 
  spread(parm, val_sc) %>% 
  group_by(draw) %>% 
  nest(year, phi, gammaII, tau) %>% 
  mutate(phi = map(data, fit_covs_phi_hsc),
         gammaII = map(data, fit_covs_gammaII),
         tau = map(data, fit_covs_tau)) 

# summarize results
sa_path_res <- sand_path %>% 
  select(-data) %>% 
  gather(parm, mod, 2:4) %>% 
  mutate(res = map(mod, tidy)) %>% 
  select(-mod) %>% 
  unnest() %>% 
  ungroup()

# beta posterior distributions
sa_path_res %>% 
  filter(parm == "phi") %>% 
  ggplot(aes(x = estimate)) +
  geom_density() +
  geom_vline(xintercept = 0) +
  facet_wrap(~term, scales = "free") +
  ggtitle("survival")

sa_path_res %>% 
  filter(parm == "gammaII") %>% 
  ggplot(aes(x = estimate)) +
  geom_density() +
  geom_vline(xintercept = 0) +
  facet_wrap(~term, scales = "free") +
  ggtitle("temp emi")

sa_path_res %>% 
  filter(parm == "tau") %>% 
  ggplot(aes(x = estimate)) +
  geom_density() +
  geom_vline(xintercept = 0) +
  facet_wrap(~term, scales = "free") +
  ggtitle("residency")

# all ----

rk_path_res2 <- rk_path_res %>% 
  mutate(species = "rekn")

rt_path_res2 <- rt_path_res %>% 
  mutate(species = "rutu")

all_path_res <- sa_path_res %>% 
  mutate(species = "sand") %>% 
  bind_rows(rt_path_res2) %>% 
  bind_rows(rk_path_res2)

# saveRDS(all_path_res, "output/all_sp_covars.rds")
# all_path_res = readRDS("output/all_sp_covars.rds")


all_sum <- all_path_res %>% 
  filter(term != "(Intercept)") %>% 
  group_by(species, parm, term) %>% 
  summarize(est = mean(estimate),
            lower = quantile(estimate, 0.025),
            upper = quantile(estimate, 0.975),
            pos = mean(estimate > 0),
            neg = mean(estimate < 0),
            f = ifelse(est > 0, pos, neg)) %>% 
  select(-pos, -neg)

# table 2 - beta estimates 
all_sum %>% 
  arrange(-f) %>% 
  print(n = 50)

all_sum %>% 
  mutate_if(is.double, ~round(., 2)) %>% 
  mutate(term = substr(term, 5, nchar(term))) %>%
  arrange(species, -f) %>% 
  write.csv("tables/table2_covars.csv", row.names = F)


# predicted relationships ---- 

predict_parm = function(estimate){
  x = seq(-2, 2, 0.01)
  y = estimate*x
  return(data.frame(x = x, y = y))
}

toplot <- all_path_res %>% 
  mutate(parm = ifelse(str_detect(parm, "phi"), "phi", parm)) %>% 
  select(species, draw, parm, term, estimate) %>% 
  mutate(pred = map(estimate, predict_parm)) 

covs_weather_pres <- toplot %>% 
  select(-pred) %>% 
  filter(term %in% c("covars$n15", "covars$storm")) %>% 
  group_by(species, parm, term) %>% 
  mean_qi(est = estimate) %>% 
  ungroup() %>% 
  mutate(parmlab = c("phi" = "Effect of severe weather\non survival probability",
                     "tau" = "Effect of water temperature\non residency probability")[parm]) %>% 
  ggplot(aes(y = est, x = term, col = species)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_linerange(aes(ymin = .lower, ymax = .upper), lwd = 2,
                 position = position_dodge(width = 0.8)) +
  geom_point(size = 4, position = position_dodge(width = 0.8)) +
  coord_flip() +
  facet_wrap(~parmlab, scales = "free") +
  ylim(-1.5, 1.5) +
  ylab("Standardized regression coefficient") +
  xlab("") +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(margin = margin(1,1,2,1),
                                  size = 16),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18)) +
  scale_color_manual(values = c("#b95b22", "#436e79", "#cd9f5a"),
                     labels = c("Red knot", "Ruddy turnstone", "Sanderling"),
                     name = "")
covs_weather_pres

# save_plot(covs_weather_pres, file = "figs/covs_weather_pres.jpg",
#           base_width = 10, base_height = 6)

covs_hsc_pres <- toplot %>% 
  select(-pred) %>% 
  filter(term == "covars$hsc") %>% 
  group_by(species, parm, term) %>% 
  mean_qi(est = estimate) %>% 
  ungroup() %>% 
  mutate(parmlab = c("gammaII" = "Returning",
                     "phi" = "Survival",
                     "tau" = "Residency")[parm],
         covlab = c("covars$n15" = "Water\ntemperature",
                    "covars$hsc_scaled" = "Hsc egg\navailability",
                    "covars$rate" = "Rate of\nmass gain",
                    "covars$storm" = "Severe\nweather")[term]) %>% 
  ggplot(aes(y = est, x = covlab, col = species)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_linerange(aes(ymin = .lower, ymax = .upper), lwd = 2,
                 position = position_dodge(width = 0.8)) +
  geom_point(size = 4, position = position_dodge(width = 0.8)) +
  coord_flip() +
  facet_wrap(~parmlab, scales = "free") +
  scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-1.2, 1.2)) +
  ylab("Standardized regression coefficient") +
  xlab("") +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(margin = margin(1,1,2,1),
                                  size = 16),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 18),
        axis.title = element_text(size = 20),
        #legend.text = element_text(size = 18),
        legend.position = "none",
        panel.spacing = unit(2, "lines")) +
  scale_color_manual(values = c("#b95b22", "#436e79", "#cd9f5a"),
                     labels = c("Red knot", "Ruddy turnstone", "Sanderling"),
                     name = "") 
covs_hsc_pres

# save_plot(covs_hsc_pres, file = "figs/covs_hsc_pres.jpg", 
#           base_width = 12, base_height = 6)



# convert scaled y back to probabilities
rekn_res %>% 
  group_by(parm, draw) %>% 
  summarize(mean = mean(val),
            sd = sd(val)) %>% 
  ungroup() %>% 
  mutate(species = "rekn") -> rekn_convert

rutu_res %>% 
  group_by(parm, draw) %>% 
  summarize(mean = mean(val),
            sd = sd(val)) %>% 
  ungroup() %>% 
  mutate(species = "rutu") -> rutu_convert

sand_res %>% 
  group_by(parm, draw) %>% 
  summarize(mean = mean(val),
            sd = sd(val)) %>% 
  ungroup() %>% 
  mutate(species = "sand") -> sand_convert  

rekn_convert %>% 
  bind_rows(rutu_convert) %>% 
  bind_rows(sand_convert)  -> convert_val

toplot3 <- toplot %>% 
  filter(term != "(Intercept)") %>% 
  full_join(convert_val) %>% 
  select(species, parm, draw, term, pred, mean, sd) %>% 
  unnest() %>% 
  mutate(prob = (y*sd)+mean) %>% 
  group_by(species, parm, term, x) %>% 
  mean_qi(est = prob) %>% 
  ungroup()


# new fig3 ----


rekn_res %>%
  filter(parm %in% c("phi", "gammaII")) %>%
  group_by(year, parm) %>%
  mean_qi(est = val) %>%
  ungroup() %>%
  mutate(species = "rekn") -> rekn_est

rutu_res %>%
  filter(parm %in% c("phi", "gammaII")) %>%
  group_by(year, parm) %>%
  mean_qi(est = val) %>%
  ungroup() %>%
  mutate(species = "rutu") -> rutu_est

sand_res %>%
  filter(parm %in% c("phi", "gammaII")) %>%
  group_by(year, parm) %>%
  mean_qi(est = val) %>%
  ungroup() %>%
  mutate(species = "sand") -> sand_est

birds <- rekn_est %>%
  bind_rows(rekn_est) %>%
  bind_rows(rutu_est) %>%
  bind_rows(sand_est)

hsc_est <- covars %>% 
  dplyr::select(hsc) %>% 
  mutate(year = c(1:13)) %>% 
  full_join(birds) %>% 
  rename(x = hsc,
         parmest = est,
         parmlow = .lower,
         parmhigh = .upper) %>% 
  mutate(species = c("rekn" = "Red knot",
                     "rutu" = "Ruddy turnstone",
                     "sand" = "Sanderling")[species],
         parm = c("phi" = "Apparent annual\nsurvival probability",
                  "gammaII" = "Probability of returning")[parm]) %>% 
  distinct()

# figure 3 ---
fig3 <- toplot3 %>% 
  filter(parm %in% c("phi", "gammaII") & term == "covars$hsc") %>% 
  mutate(species = c("rekn" = "Red knot",
                     "rutu" = "Ruddy turnstone",
                     "sand" = "Sanderling")[species],
         parm = c("phi" = "Apparent annual\nsurvival probability",
                  "gammaII" = "Probability of returning")[parm]) %>% 
  # full_join(hsc_est) %>%
  ggplot(aes(x = x)) + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = "gray80", alpha = 0.5) + 
  geom_line(aes(y = est), lwd = 1, col = "gray40") +
  geom_linerange(data = hsc_est, aes(ymin = parmlow, ymax = parmhigh)) +
  geom_point(data = hsc_est, aes(y = parmest), size = 2) +
  facet_grid(parm ~ species, scales = "free_y", switch = "y") +
  xlab("Relative horseshoe crab egg abundance") +
  ylab("") +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 14, margin = margin(1,1,2,1)),
        strip.placement = "outside",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)
fig3

save_plot(fig3, file = "fig3_phi-gamma-hsc.tiff",
          base_width = 10, base_height = 6)

