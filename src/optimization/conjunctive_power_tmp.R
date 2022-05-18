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


optim_G_constrain <- function(alpha, mp, w, constraint.G=NULL, tolerance=1e-6) {
  n.hypo = length(mp)
  optim.G = diag(0, n.hypo)
  if (is.null(constraint.G)) {
    optim.G = optim_G(alpha, mp, w)
  }
  else {
    stopifnot(dim(constraint.G)==c(n.hypo, n.hypo))
    for (i in 1:n.hypo) {
      stopifnot(sum(constraint.G[i, ], na.rm=T) <= 1)
      # Check if row sum up to 1 already
      if (sum(constraint.G[i, ], na.rm=T) == 1) {
        optim.G[i, ] = constraint.G[i, ]
        optim.G[i, is.na(optim.G[i, ])] = 0
        next
      }
      # Check if row is well defined
      if (sum(is.na(constraint.G[i, ])) == 0) {
        stopifnot(sum(constraint.G[i, ]) == 1)
        optim.G[i, ] = constraint.G[i, ]
        next
      }
      if (sum(is.na(constraint.G[i, ])) == 1) {
        optim.G[i, ] = constraint.G[i, ]
        constraint.G[i, is.na(constraint.G[i, ])] = 1 - sum(constraint.G[i, ], na.rm=T)
        next
      }
      # Optimize w constraint
      lb = w[-i]
      ub = rep(1, n.hypo-1)
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
      else {
        constraint.w = rep(NA, n.hypo-1)
      }
      if (w[i] == 0) {
        w[i] = tolerance
      }
      # print(lb)
      # print(ub)
      # print(constraint.w)
      w.prime = optim_w(alpha, mp[-i], lb=lb, ub=ub, constraint.w=constraint.w)$w
      # print(w.prime)
      optim.G[i, -i] = (w.prime - w[-i])/w[i]
    }
  }
  return(optim.G)
}
