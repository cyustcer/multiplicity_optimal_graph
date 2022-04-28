# Import libraries
library(mvtnorm)
library(tidyverse)
library(nloptr)

########################
### Independent Case ###
########################

power <- function(alpha, eta=2.5) {
  1 - pnorm(qnorm(1-alpha) - eta)
}


# solve_eta <- function(alpha, p) {
#   equation <- function(alpha, eta) {
#     power(alpha, eta) - p
#   }
#   root = uniroot(equation, alpha=alpha, lower=-1000, upper=1000)$root
#   root
# }

solve_eta <- function(alpha, p) {
  qnorm(1-alpha) - qnorm(1-p)
}

disjunctive_power <- function(w=c(0.5, 0.5), alpha=0.025, mp=c(0.9, 0.9)) {
  n = length(mp)
  stopifnot(length(w) == n)
  etas = c()
  for (p in mp) {
    etas = c(etas, solve_eta(alpha=alpha, p=p))
  }
  fail = 1
  for (i in 1:n) {
    fail = fail * (1 - power(w[i]*alpha, eta=etas[i]))
  }
  1 - fail
}

loss <- function(w, alpha=ALPHA, mp=MP, min.w=1e-8) {
  -disjunctive_power(w=w, alpha=alpha, mp=mp)
}

loss_d <- function(w, k, alpha=ALPHA, mp=MP, min.w=1e-8) {
  n = length(mp)
  stopifnot(length(w) == n)
  stopifnot(k <= n)
  if (w[k] <= min.w) {
    w[k] = min.w
  }
  etas = c()
  for (p in mp) {
    etas = c(etas, solve_eta(alpha=alpha, p=p))
  }
  derivative = 1
  for (i in 1:n) {
    if (i != k) {
      derivative = derivative * (1 - power(w[i]*alpha, eta=etas[i]))
    }
    else {
      derivative = -derivative * alpha * dnorm(qnorm(1-w[i]*alpha)-etas[i]) / dnorm(qnorm(1-w[i]*alpha))
    }
  }
  derivative
}

loss_grad <- function(w, alpha=ALPHA, mp=MP, min.w=1e-8) {
  n = length(mp)
  stopifnot(length(w) == n)
  grad = c()
  for (k in 1:n) {
    grad = c(grad, loss_d(w, k, alpha=alpha, mp=mp, min.w=min.w))
  }
  grad
}

eqn <- function(w, alpha, mp) {
  return(c(sum(w)-1))
}

eqn_grad <- function(w, alpha, mp) {
  return(rep(1, length(w)))
}

#######################
### Correlated Case ###
#######################

power_corr <- function(mus, Sigma, alphas) {
  zs = c()
  for (alpha in alphas) {
    zs = c(zs, qnorm(1-alpha))
  }
  1 - pmvnorm(upper=zs, mean=mus, sigma=Sigma)
}

# solve_mu <- function(alpha, p) {
#   equation <- function(alpha, mu) {
#     z = qnorm(1-alpha)
#     1 - pnorm(z, mean=mu) - p
#   }
#   root = uniroot(equation, alpha=alpha, lower=-1000, upper=1000)$root
#   root
# }

solve_mu <- function(alpha, p) {
  qnorm(1-alpha) - qnorm(1-p)
}

sigma_matrix <- function(rho, n) {
  sigma = matrix(rho, n, n)
  diag(sigma) = 1
  sigma
}

disjunctive_power_corr <- function(w=c(0.5, 0.5), alpha=0.025, mp=c(0.9, 0.9), rho=0) {
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
  -disjunctive_power_corr(w, alpha=alpha, mp=mp, rho=rho)
}

conditional_distribution <- function(z, i, mus, Sigma) {
  mus_i = mus[-i] + Sigma[-i, i] * (z[i]-mus[i])
  Sigma_i = Sigma[-i, -i] - Sigma[-i, i] %*% t(Sigma[i, -i])
  z_i = z[-i]
  pmvnorm(upper=z_i, mean=as.vector(mus_i), sigma=Sigma_i)
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
  return(c(sum(w)-1))
}

eqn_grad_corr <- function(w, alpha, mp, rho, min.w=1e-8) {
  return(rep(1, length(w)))
}


# Optimization of w
opts <- list("algorithm" = "NLOPT_LD_SLSQP",
             "xtol_rel" = 0,
             "maxeval" = 1e5)

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
    res <- nloptr(x0=initial_w,
                  eval_f=loss,
                  eval_grad_f = loss_grad,
                  alpha=alpha, mp=mp.p,
                  lb=lb, ub=ub, min.w=min.w,
                  eval_g_eq = eqn,
                  eval_jac_g_eq = eqn_grad,
                  opts=optim_opts)
  }
  else {
    res <- nloptr(x0=initial_w,
                  eval_f=loss_corr,
                  eval_grad_f = loss_grad_corr,
                  alpha=alpha, mp=mp.p, rho=rho,
                  lb=lb, ub=ub, min.w=min.w,
                  eval_g_eq = eqn_corr,
                  eval_jac_g_eq = eqn_grad_corr,
                  opts=opts)
  }
  w.p = res$solution
  w = rep(0, length(mp))
  w[constraint.w != 0 | is.na(constraint.w)] = w.p
  return(list(w = w,
              optima = res))
}

create_initial_ws <- function(n) {
  cand.ws <- NULL
  for (i in 1:n) {
    cand = list(c(0, 1))
    names(cand) = paste0("w", i)
    cand.ws <- c(cand.ws, cand)
  }
  initial.ws <- cand.ws %>% expand.grid()
  initial.ws = initial.ws[rowSums(initial.ws) != 0, ]
  initial.ws = initial.ws / rowSums(initial.ws)
  row.names(initial.ws) <- NULL
  initial.ws
}

go_optim_w <- function(alpha, mp, rho=NULL,
                       constraint.w=NULL, lb=NULL, ub=NULL,
                       min.w=1e-8, optim_opts=opts) {
  n = length(mp)
  initial.ws = create_initial_ws(n)
  N = nrow(initial.ws)
  optimas <- data.frame(matrix(ncol = n+1, nrow = 0))
  colnames(optimas) <- c(colnames(initial.ws), 'optimal_value')
  for (i in 1:N) {
    initial.w = as.matrix(initial.ws)[i, ]
    # print(initial.w)
    res = optim_w(alpha=alpha, mp=mp, rho=rho,
                  initial_w=initial.w, constraint.w=constraint.w,
                  lb=lb, ub=ub, min.w=min.w,
                  optim_opts=optim_opts)
    # optimal.w = round(res$w, digits=pcs.digits)
    optimal.w = res$w
    optimal.value = -res$optima$objective
    optimas[i, ] <- c(optimal.w, optimal.value)
  }
  optimas %>% arrange(desc(optimal_value))
}

# approximate view
optimas %>% mutate(across(w1:w6, ~round(.x, digits=3)))
