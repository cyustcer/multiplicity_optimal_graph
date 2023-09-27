source(here::here('src/utils/misc.R'))

cpc <- function(mus, Sigma, alphas) {
  zs = c()
  for (alpha in alphas) {
    zs = c(zs, qnorm(1-alpha))
  }
  set.seed(2022)
  pmvnorm(lower=zs, mean=mus, sigma=Sigma,
          algorithm = mvtnorm::GenzBretz(maxpts = 1e6,
                                         abseps = 1e-6,
                                         releps = 1e-6))
}

conjunctive_power_corr <- function(w=c(0.5, 0.5), alpha=0.025, mp=rep(0.9, 2), rho=0) {
  n = length(mp)
  mus = c()
  for (p in mp) {
    mus = c(mus, solve_eta(alpha=alpha, p=p))
  }
  if (length(rho) == 1) {
    Sigma = sigma_matrix(rho, n)
  }
  else {
    Sigma = rho
  }
  cpc(mus=mus, Sigma=Sigma, alphas=alpha*w)
}


loss_cpc <- function(w, alpha=ALPHA, mp=MP, rho=RHO, min.w=1e-8) {
  -conjunctive_power_corr(w, alpha=alpha, mp=mp, rho=rho)
}


loss_cpc_d <- function(w, k, alpha, mp, rho, min.w=1e-8) {
  n <- length(mp)
  stopifnot(length(w) == n)
  stopifnot(k <= n)
  if (w[k] <= min.w) {
    w[k] <- min.w
  }
  mus = c()
  for (p in mp) {
    mus = c(mus, solve_eta(alpha=alpha, p=p))
  }
  zs = c()
  for (t in w) {
    zs = c(zs, qnorm(1-t*alpha))
  }
  Sigma = sigma_matrix(rho, n)
  cd = conditional_distribution(zs, i = k, mus = mus,
                                Sigma = Sigma)
  derivative = pmvnorm(lower = zs[-k],
                       mean = cd$mean,
                       sigma = cd$covar,
                       algorithm = mvtnorm::GenzBretz(maxpts = 1e6,
                                                      abseps = 1e-6,
                                                      releps = 1e-6))
  derivative = derivative * inverse_derivative(w = w,
                                               alpha = alpha,
                                               i = k)
  derivative = derivative * dnorm(zs[k], mean=mus[k])
  derivative
}


#' A function to calculate partial derivative of loss function
#'
#' @param w vector, weight assigned to each hypotheses
#' @param alpha float, type I error
#' @param mp vector of float, marginal power of each  hypotheses
#' @param rho float or matrix, correlation
#'
#' @return vector of partial derivatives
#' @export
loss_cpc_grad <- function(w, alpha, mp, rho, min.w=1e-8) {
  n <- length(mp)
  stopifnot(length(w) == n)
  grad <- c()
  for (k in 1:n) {
    grad <- c(grad, loss_cpc_d(w, k, alpha=alpha,
                               mp=mp, rho=rho,
                               min.w=min.w))
  }
  grad
}


#' A function to calculate constraint
#'
#' @param w vector, weight assigned to each hypotheses
#' @param alpha float, type I error
#' @param mp vector of float, marginal power of each  hypotheses
#'
#' @return constraints
#' @export
eqn_corr <- function(w, alpha, mp, rho, min.w) {
  return(sum(w) - 1.)
}


#' A function to calculate gradient of constraint
#'
#' @param w vector, weight assigned to each hypotheses
#' @param alpha float, type I error
#' @param mp vector of float, marginal power of each  hypotheses
#'
#' @return gradient of constraints
#' @export
eqn_corr_grad <- function(w, alpha, mp, rho, min.w) {
  return(rep(1, length(w)))
}
