library(mvtnorm)
library(tidyverse)
library(nloptr)
source("src/utils/misc.R")


#' A function
#'
#' @param alpha type I error
#' @param eta 
#'
#' @return The marginal power
#' @export
marginal_power <- function(alpha, eta=2.5) {
  1 - pnorm(qnorm(1 - alpha) - eta)
}


#' A function to calculate disjunctive power of a Bonferroni test
#'
#' @param w vector, weight assigned to each hypotheses
#' @param alpha float, type I error
#' @param mp vector of float, marginal power of each  hypotheses
#'
#' @return The disjunctive power
#' @export
disjunctive_power <- function(w, alpha, mp) {
  n <- length(mp)
  stopifnot(length(w) == n)
  etas <- c()
  for (p in mp) {
    etas <- c(etas, solve_eta(alpha = alpha, p = p))
  }
  fail <- 1
  for (i in 1:n) {
    fail <- fail * (1 - marginal_power(w[i] * alpha, eta = etas[i]))
  }
  1 - fail
}


#' A function to calculate loss for disjunctive power, i.e.
#' negative of disjunctive power
#'
#' @param w vector, weight assigned to each hypotheses
#' @param alpha float, type I error
#' @param mp vector of float, marginal power of each  hypotheses
#'
#' @return loss
#' @export
loss_dp <- function(w, alpha, mp, min.w = 1e-8) {
  -disjunctive_power(w = w, alpha = alpha, mp = mp)
}


#' A function to calculate the kth partial derivative of
#' loss function
#'
#' @param w vector, weight assigned to each hypotheses
#' @param k int, the kth variable, i.e. kth weight
#' @param alpha float, type I error
#' @param mp vector of float, marginal power of each  hypotheses
#'
#' @return kth partial derivative of loss function
#' @export
loss_dp_d <- function(w, k, alpha, mp, min.w = 1e-8) {
  n <- length(mp)
  stopifnot(length(w) == n)
  stopifnot(k <= n)
  if (w[k] <= min.w) {
    w[k] <- min.w
  }
  etas <- c()
  for (p in mp) {
    etas <- c(etas, solve_eta(alpha = alpha, p = p))
  }
  derivative <- 1
  for (i in 1:n) {
    if (i != k) {
      derivative <- derivative *
        (1 - marginal_power(w[i] * alpha, eta = etas[i]))
    }
    else {
      derivative <- -derivative * alpha *
        dnorm(qnorm(1 - w[i] * alpha) - etas[i]) /
        dnorm(qnorm(1 - w[i] * alpha))
    }
  }
  derivative
}

#' A function to calculate partial derivative of loss function
#'
#' @param w vector, weight assigned to each hypotheses
#' @param alpha float, type I error
#' @param mp vector of float, marginal power of each  hypotheses
#'
#' @return vector of partial derivatives
#' @export
loss_dp_grad <- function(w, alpha, mp, min.w = 1e-8) {
  n <- length(mp)
  stopifnot(length(w) == n)
  grad <- c()
  for (k in 1:n) {
    grad <- c(grad, loss_dp_d(w, k, alpha = alpha, mp = mp, min.w = min.w))
  }
  grad
}
