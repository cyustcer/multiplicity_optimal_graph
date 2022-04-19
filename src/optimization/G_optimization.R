source(here::here("src/optimization/disjunctive_power.R"))

optim_G <- function(alpha, mp, w) {
  n.hypo = length(mp)
  optim.G = diag(0, n.hypo)
  for (i in 1:n.hypo) {
    w.prime = optim_w(alpha, mp[-i], lb=w[-i])$w
    optim.G[i, -i] = (w.prime - w[-i])/w[i]
  }
  return(optim.G)
}


optim_G_constrain <- function(alpha, mp, w, constraint.G=NULL) {
  n.hypo = length(mp)
  optim.G = diag(0, n.hypo)
  for (i in 1:n.hypo) {
    lb = w[-i]
    ub = rep(1, n.hypo-1)
    if (!is.null(constraint.G)) {
      constraint.w = constraint.G[i, -i]
      if (sum(lb == 0) > 0) {
        if (sum(lb) == 1) {
          constraint.w[lb == 0] = 0
        }
        else {
          lb[lb == 0] = 1e-6
        }
        lb = lb[is.na(constraint.w) | constraint.w != 0]
        ub = ub[is.na(constraint.w) | constraint.w != 0]
      }
    }
    else {
      constraint.w = rep(NA, n.hypo-1)
    }
    print(lb)
    print(ub)
    print(constraint.w)
    w.prime = optim_w(alpha, mp[-i], lb=lb, ub=ub, constraint.w=constraint.w)$w
    print(w.prime)
    optim.G[i, -i] = (w.prime - w[-i])/w[i]
  }
  return(optim.G)
}
