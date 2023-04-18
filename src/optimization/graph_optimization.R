source("src/optimization/w_optimization.R")


optim_row_dp <- function(i, alpha, w, mp, rho=NULL, constraint.G=NULL, eps=1e-3) {
  # Check precision of w
  if (sum(w) != 1) {
    w = w + w / sum(w) * (1 - sum(w.optim))
  }
  n.hypo = length(mp)
  if (length(rho) != 1) {
    rho = rho[-i, -i]
  }
  if (w[i] == 0) {
    w[i] = eps
    for (j in 1:n.hypo) {
      if (j != i) {
        w[j] = w[j] - w[j] * (eps)
      }
    }
  }
  lb = as.numeric(w[-i])
  ub = rep(1, n.hypo -1)
  if (!is.null(constraint.G)) {
    constraint.w = constraint.G[i, -i] * w[i] + w[-i]
    if (sum(lb == 0) > 0) {  # if there exists lower bound equals 0
      if (sum(lb) == 1) {
        constraint.w[lb == 0] = 0
      }
    }
  }
  else {
    constraint.w = rep(NA, n.hypo-1)
  }
  res = go_optim_w_dp(alpha, mp[-i],
                      rho=rho, lb=lb, ub=ub,
                      constraint.w=constraint.w)
  w.prime = res[1, ] %>%
    select(starts_with("w")) %>%
    as.numeric()
  Gi = (w.prime - w[-i]) / w[i]
  return(Gi)
}


# n = 8
# for (w3 in rep(1, n) /  10^seq(n)) {
#   w = optim_w_dp(alpha=0.025, mp=c(0.9, 0.8, 0.7), rho=0,
#                  constraint.w=c(NA, NA, w3))$w
#   print(optim_row_dp(3, alpha=0.025, w=w, mp=c(0.9, 0.8, 0.7)))
# }

optim_G_dp <- function(alpha, w, mp, rho=NULL, constraint.G=NULL, eps=1e-4) {
  n.hypo = length(mp)
  optim.G = diag(0, n.hypo)
  for (i in 1:n.hypo) {
    Gi = optim_row_dp(i, alpha=alpha, w=w, mp=mp,
                      rho=rho, constraint.G=constraint.G, eps=eps)
    optim.G[i, -i] = Gi
  }
  return(optim.G)
}


# optim_G_dp_2 <- function(alpha, w, mp, rho=NULL, constraint.G=NULL) {
#   n.hypo = length(mp)
#   W.prime = diag(0, n.hypo)
#   optim.G = diag(0, n.hypo)
#   for (i in 1:n.hypo) {
#     w.prime = optim_row_dp(i, alpha=alpha, w=w, mp=mp,
#                            rho=rho, constraint.G=constraint.G)
#     W.prime[i, -i] = w.prime
#   }
#   print(W.prime)
#   for (i in 1:n.hypo) {
#     for (j in 1:n.hypo) {
#       if (i != j) {
#         optim.G[i, j] = (W.prime[j, i] - w[i]) / W.prime[i, j]
#       }
#     }
#   }
#   return(optim.G)
# }
