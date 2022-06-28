source('src/utils/misc.R')

#' A function to transfer single correlation to matrix
#'
#' @param rho float, single correlation
#' @param n int, number of hypotheses
#'
#' @return correlation matrix
#' @export
sigma_matrix <- function(rho, n) {
  sigma <- matrix(rho, n, n)
  diag(sigma) <- 1  # make diagnal to 1
  sigma
}

dpc <- function(mus, Sigma, alphas) {
  zs = c()
  for (alpha in alphas) {
    zs = c(zs, qnorm(1 - alpha))
  }
  set.seed(2022)
  1 - pmvnorm(upper = zs, mean = mus, sigma = Sigma,
              algorithm = GenzBretz(maxpts = 1e6, 
                                    abseps = 1e-6, 
                                    releps = 1e-6))
}

#' A function to calculate correlated disjunctive power of
#' a Bonferroni test
#'
#' @param w vector, weight assigned to each hypotheses
#' @param alpha float, type I error
#' @param mp vector of float, marginal power of each  hypotheses
#' @param rho float or matrix, correlation or correlation matrix
#'
#' @return The disjunctive power
#' @export
disjunctive_power_corr <- function(w, alpha, mp, rho=0) {
  n <- length(mp)
  mus <- c()
  for (p in mp) {
    mus <- c(mus, solve_eta(alpha = alpha, p = p))
  }
  if (length(rho) == 1) {
    # If single correlation is provided
    Sigma <- sigma_matrix(rho, n)
  }
  else {
    Sigma <- rho
  }
  dpc(mus = mus, Sigma = Sigma, alphas = alpha*w)
}


#' A function to calculate loss for disjunctive power, i.e.
#' negative of disjunctive power
#'
#' @param w vector, weight assigned to each hypotheses
#' @param alpha float, type I error
#' @param mp vector of float, marginal power of each  hypotheses
#' @param rho float or matrix, correlation
#'
#' @return loss
#' @export
loss_dpc <- function(w, alpha, mp, rho, min.w=1e-8) {
  -disjunctive_power_corr(w, alpha=alpha, mp=mp, rho=rho)
}


#' A function to calculate the kth partial derivative of
#' loss function
#'
#' @param w vector, weight assigned to each hypotheses
#' @param k int, the kth variable, i.e. kth weight
#' @param alpha float, type I error
#' @param mp vector of float, marginal power of each  hypotheses
#' @param rho float or matrix, correlation
#'
#' @return kth partial derivative of loss function
#' @export
loss_dpc_d <- function(w, k, alpha, mp, rho, min.w=1e-8) {
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
  derivative = pmvnorm(upper = zs[-k],
                       mean = cd$mean,
                       sigma = cd$covar,
                       algorithm = GenzBretz(maxpts = 1e6,
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
loss_dpc_grad <- function(w, alpha, mp, rho, min.w=1e-8) {
  n <- length(mp)
  stopifnot(length(w) == n)
  grad <- c()
  for (k in 1:n) {
    grad <- c(grad, loss_dpc_d(w, k, alpha=alpha,
                                mp=mp, rho=rho,
                                min.w=min.w))
  }
  grad
}
