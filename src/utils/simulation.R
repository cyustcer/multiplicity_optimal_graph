library(MASS)
source(here::here('src/utils/power.R'))

#####################################
### Functions for data simulation ###
#####################################


### Marginal treatment effect
# marginal treatment effect vector calculated based on one-sided type I error
margin.trt.eff <- function(type.i.error, margin.power) {
  # type.i.error: float, one-sided type I error
  # margin.power: vector, marginal power of each hypothese
  trt.eff = qnorm(1-type.i.error) - qnorm(margin.power, lower.tail=FALSE)
  return(trt.eff)
}


simulate.pval <- function(type.i.error, margin.power, sigma, n.sim) {
  # Calculate marginal treatment effect
  mte = margin.trt.eff(type.i.error=type.i.error, margin.power=margin.power)
  # Simulate treatment effects
  trt.eff.sim = t(mvrnorm(n=n.sim, mte, Sigma=sigma))
  # Calculate p-value from simulated treatment effect
  pval.sim = t(pnorm(trt.eff.sim, lower.tail=FALSE))
  return(pval.sim)
}


simulate.grid.power <- function(n.sim, margin.power,
                                ws, Gs, sigma,
                                type.i.error, weight,
                                seed=2022, n.cluster=1) {
  set.seed(seed)
  pvals = simulate.pval(type.i.error, margin.power, sigma, n.sim)
  n.graph = dim(ws)[1]
  n.hypo = dim(ws)[2]
  data.sim = cbind(ws, matrix(aperm(Gs, c(3, 2, 1)), nrow=n.graph))
  data.sim = data.frame(data.sim)
  colnames(data.sim) = c(paste0('w', 1:n.hypo),
                         paste0('G', as.vector(sapply(1:n.hypo, 
                                                      function(x){paste0(x,"_", 1:n.hypo)}))))
  data.sim <- data.sim %>%
    add_column("disjunctive.power"=0, "conjunctive.power"=0, "weighted.power"=0)
  for (i.graph in 1:n.graph) {
    out.success = graph.power(as.vector(ws[i.graph, ]),
                              as.matrix(Gs[, , i.graph]),
                              type.i.error, pvals)
    power.disjunctive = disjunctive.power(out.success)
    power.conjunctive = conjunctive.power(out.success)
    power.weighted = weighted.power(out.success, weight=weight)
    data.sim[i.graph, "disjunctive.power"] = power.disjunctive
    data.sim[i.graph, "conjunctive.power"] = power.conjunctive
    data.sim[i.graph, "weighted.power"] = power.weighted
  }
  return(data.sim)
}
