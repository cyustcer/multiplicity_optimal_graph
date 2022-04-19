source(here::here('src/utils/', 'tools.R'))

# dir.sim = "outputs/grids/mp-0.5-0.8-0.9_weight-0.333333333333333-0.333333333333333-0.333333333333333_corr-0_n.sim-6_w.step-0.01_w.range-0-1_G.step-0.2/"
dir.sim = "outputs/grids/mp-0.5-0.8-0.9_weight-0.333333333333333-0.333333333333333-0.333333333333333_corr-0_n.sim-6_w.specified-0.09965126-0.36192224-0.53842651_G.step-0.01/"
sim.grid = load.sim.grid(dir.sim)

sim.grid %>% filter(disjunctive.power==max(disjunctive.power))
sim.grid %>% filter(weighted.power==max(weighted.power))
sim.grid %>% filter(conjunctive.power==max(conjunctive.power))


sim.grid %>% filter(disjunctive.power==max(disjunctive.power)) %>%
  filter(weighted.power==max(weighted.power))

sim.grid %>%
  rowwise() %>% mutate(is_connected=is.connected(G1_2, G1_3, G2_1, G2_3, G3_1, G3_2)) %>%
  filter(is_connected)
