library(mvtnorm)
library(tidyverse)

#' A function
#'
#' @param alpha type I error
#' @param eta 
#'
#' @return The marginal power
#' @export
solve_eta <- function(alpha, p) {
  qnorm(1 - alpha) - qnorm(1 - p)
}


#' A function to calculate conditional distribution
#'
#' @param alpha type I error
#' @param eta 
#'
#' @return The marginal power
#' @export
conditional_distribution <- function(z, i, mus, Sigma) {
  mus_i <- mus[-i] + Sigma[-i, i] * (z[i]-mus[i])
  Sigma_i <- Sigma[-i, -i] - Sigma[-i, i] %*% t(Sigma[i, -i])
  return(list(mean = as.vector(mus_i),
              covar = Sigma_i))
}


#' A function to calculate conditional distribution
#'
#' @param alpha type I error
#' @param eta 
#'
#' @return The marginal power
#' @export
inverse_derivative <- function(w, alpha, i) {
  -alpha / dnorm(qnorm(1 - w[i] * alpha))
}


#' A function to calculate constraint
#'
#' @param w vector, weight assigned to each hypotheses
#' @param alpha float, type I error
#' @param mp vector of float, marginal power of each  hypotheses
#'
#' @return constraints
#' @export
eqn <- function(w, alpha, mp, min.w) {
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
eqn_grad <- function(w, alpha, mp, min.w) {
  return(rep(1, length(w)))
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


#' A function to calculate inequality constraint w.o. correlation
#'
#' @param w vector, weight assigned to each hypotheses
#' @param alpha float, type I error
#' @param mp vector of float, marginal power of each  hypotheses
#'
#' @return constraints
#' @export
ineqn <- function(w, alpha, mp, min.w) {
  h <- numeric(2)
  h[1] <- sum(w) - 1.
  h[2] <- 1. - sum(w)
  return(h)
}


#' A function to calculate inequality constraint w. correlation
#'
#' @param w vector, weight assigned to each hypotheses
#' @param alpha float, type I error
#' @param mp vector of float, marginal power of each  hypotheses
#'
#' @return constraints
#' @export
ineqn_corr <- function(w, alpha, mp, rho, min.w) {
  h <- numeric(2)
  h[1] <- sum(w) - 1.
  h[2] <- 1. - sum(w)
  return(h)
}
