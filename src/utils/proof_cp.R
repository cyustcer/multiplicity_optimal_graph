source('src/optimization/conjunctive_power.R')

w = seq(0, 1, by=0.01)
g_2w_corr <- function(w, rho, mp=c(0.8, 0.8), alpha=0.025) {
  conjunctive_power_corr(w=c(w, 1-w), alpha=alpha, mp=mp, rho=rho)
}


for (rho in seq(0., 1, by=0.05)) {
  warpper <- function(w) {
    cpc = c()
    for (weight in w) {
      cpc = c(cpc, g_2w_corr(weight, rho))
    }
    cpc
  }
  flag_add = rho != 0.
  curve(warpper, add=flag_add, ylim=c(-0.05, .75))
}
