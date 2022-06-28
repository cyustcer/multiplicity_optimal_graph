# Import libraries
library(mvtnorm)
library(tidyverse)
library(nloptr)

#######################
### Correlated Case ###
#######################

power_corr <- function(mus, Sigma, alphas) {
  zs = c()
  for (alpha in alphas) {
    zs = c(zs, qnorm(1-alpha))
  }
  set.seed(2022)
  pmvnorm(lower=zs, mean=mus, sigma=Sigma)
}

solve_mu <- function(alpha, p) {
  qnorm(1-alpha) - qnorm(1-p)
}

sigma_matrix <- function(rho, n) {
  sigma = matrix(rho, n, n)
  diag(sigma) = 1
  sigma
}

conjunctive_power_corr <- function(w=c(0.5, 0.5), alpha=0.025, mp=rep(0.9, 2), rho=0) {
  n = length(mp)
  mus = c()
  for (p in mp) {
    mus = c(mus, solve_mu(alpha=alpha, p=p))
  }
  if (length(rho) == 1) {
    Sigma = sigma_matrix(rho, n)
  }
  else {
    Sigma = rho
  }
  power_corr(mus=mus, Sigma=Sigma, alphas=alpha*w)
}

loss_corr <- function(w, alpha=ALPHA, mp=MP, rho=RHO, min.w=1e-8) {
  -conjunctive_power_corr(w, alpha=alpha, mp=mp, rho=rho)
}

conditional_distribution <- function(z, i, mus, Sigma) {
  mus_i = mus[-i] + Sigma[-i, i] * (z[i]-mus[i])
  Sigma_i = Sigma[-i, -i] - Sigma[-i, i] %*% t(Sigma[i, -i])
  z_i = z[-i]
  set.seed(2022)
  pmvnorm(lower=z_i, mean=as.vector(mus_i), sigma=Sigma_i)
}

marginal_distribution <- function(z, i, mus) {
  mus_i = mus[i]
  z_i = z[i]
  dnorm(z_i, mean=as.vector(mus_i))
}

inverse_derivative <- function(w, alpha, i) {
  -alpha / dnorm(qnorm(1-w[i]*alpha))
}

loss_d_corr <- function(w, k, alpha=ALPHA, mp=MP, rho=RHO, min.w=1e-8) {
  n = length(mp)
  stopifnot(length(w) == n)
  stopifnot(k <= n)
  if (w[k] <= min.w) {
    w[k] = min.w
  }
  mus = c()
  for (p in mp) {
    mus = c(mus, solve_mu(alpha=alpha, p=p))
  }
  zs = c()
  for (t in w) {
    zs = c(zs, qnorm(1-t*alpha))
  }
  Sigma = sigma_matrix(rho, n)
  derivative = conditional_distribution(z=zs, i=k, mus=mus, Sigma=Sigma) * inverse_derivative(w=w, alpha=alpha, i=k) * marginal_distribution(z=zs, i=k, mus=mus)
  derivative
}

loss_grad_corr <- function(w, alpha=ALPHA, mp=MP, rho=RHO, min.w=1e-8) {
  n = length(mp)
  stopifnot(length(w) == n)
  grad = c()
  for (k in 1:n) {
    grad = c(grad, loss_d_corr(w, k, alpha=alpha, mp=mp, rho=rho, min.w=min.w))
  }
  grad
}

eqn_corr <- function(w, alpha, mp, rho, min.w=1e-8) {
  return(sum(w)-1.)
}

eqn_grad_corr <- function(w, alpha, mp, rho, min.w=1e-8) {
  return(rep(1, length(w)))
}

#########################
### Optimization of w ###
#########################

opts <- list("algorithm" = "NLOPT_LD_SLSQP",
             "xtol_rel" = 0,
             "maxeval" = 1e4)

optim_w <- function(alpha, mp, rho=NULL,
                    initial_w=NULL, constraint.w=NULL,
                    lb=NULL, ub=NULL, min.w=1e-8,
                    optim_opts=opts) {
  # Filter non-zeros marginal power
  if (is.null(constraint.w)) {
    constraint.w = rep(NA, length(mp))
  }
  mp.p = mp[constraint.w != 0 | is.na(constraint.w)]
  constraint.p = constraint.w[constraint.w != 0 | is.na(constraint.w)]
  # Optimize non zeros w's
  n.hypo.p = length(mp.p)
  if (is.null(lb)) {
    # lb = rep(1e-8, n.hypo.p)
    lb = rep(0, n.hypo.p)
  }
  if (is.null(ub)) {
    ub = rep(1, n.hypo.p)
  }
  lb[!is.na(constraint.p)] = constraint.p[!is.na(constraint.p)]
  ub[!is.na(constraint.p)] = constraint.p[!is.na(constraint.p)]
  if (is.null(initial_w)) {
    remain_weight = 1 - sum(lb)
    initial_w = (ub - lb) / sum(ub - lb) * remain_weight + lb
  }
  if (is.null(rho)) {
    rho = 0
  }
  res <- nloptr(x0=initial_w,
                eval_f=loss_corr,
                eval_grad_f = loss_grad_corr,
                alpha=alpha, mp=mp.p, rho=rho,
                lb=lb, ub=ub, min.w=min.w,
                eval_g_eq = eqn_corr,
                eval_jac_g_eq = eqn_grad_corr,
                opts=opts)
  w.p = res$solution
  w = rep(0, length(mp))
  w[constraint.w != 0 | is.na(constraint.w)] = w.p
  return(list(w = w,
              optima = res))
}
