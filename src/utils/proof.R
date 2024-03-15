inverse_Mills <- function(x) {
  return(dnorm(x)/pnorm(x))
}

weight_ratio_dp <- function(w, p, alpha) {
  eta = solve_eta(alpha, p)
  mills_ratio = inverse_Mills(qnorm(1-w*alpha)-eta)
  ratio = 1/mills_ratio*dnorm(qnorm(1-w*alpha))
  return(ratio)
}

w = seq(0, 1, by=0.001)
plot(w, weight_ratio(w, p=0.14, alpha=0.1))
points(w, weight_ratio(w, p=0.8, alpha=0.025))

weights_ratio_dp <- function(w, mp, alpha) {
  ratios = c()
  n = length(w)
  for (i in 1:n) {
    ratio = weight_ratio_dp(w[i], mp[i], alpha)
    ratios = c(ratios, ratio)
  }
  return(ratios)
}

weight_ratio_cp <- function(w, p, alpha) {
  eta = solve_eta(alpha, p)
  di = 1 - pnorm(qnorm(1 - w * alpha) - eta)
  di_prime = dnorm(qnorm(1 - w * alpha)) /
    dnorm(qnorm(1 - w * alpha) - eta)
  di * di_prime
}

weights_ratio_cp <- function(w, mp, alpha) {
  ratios = c()
  n = length(w)
  for (i in 1:n) {
    ratio = weight_ratio_cp(w[i], mp[i], alpha)
    ratios = c(ratios, ratio)
  }
  return(ratios)
}

w = seq(0, 1, by=0.001)
plot(w, weight_ratio_cp(w, p=0.4, alpha=0.025))
# plot(w, weight_ratio_cp(w, p=0.8, alpha=0.025))

weights_ratio <- function(w, mp, alpha) {
  ratios = c()
  n = length(w)
  for (i in 1:n) {
    ratio = weight_ratio(w, mp[i], alpha)
    ratios = c(ratios, ratio)
  }
  return(ratios)
}


test <- function(w, p, alpha) {
  # C = solve_eta(alpha, p)
  C = qnorm(1-alpha) - qnorm(1-p)
  return(dnorm(qnorm(1-w*alpha) - C) - C * pnorm(qnorm(1-w*alpha) - C))
}

plot(w, test(w, p=0.18, alpha=0.1))

power_corr_split <- function(rho, p=0.9, n=2, alpha=0.025) {
  disjunctive_power_corr(w=rep(1/n, n), mp=rep(p, n), alpha=alpha, rho=rho)
}

n = 2
rhos = seq(0.1, 0.95, by=0.01)
dp = c()
for (rho in rhos) {
  dp = c(dp, power_corr_split(rho, n=n))
}
plot(rhos, dp)
points(rhos, dp)
abline(h=0.9)

w = seq(0, 1, by=0.01)
g_2w_corr <- function(w, mp=c(0.8, 0.8), alpha=0.025, rho=0.9) {
  disjunctive_power_corr(w=c(w, 1-w), alpha = alpha, mp = mp, rho = rho)
}

dpc = c()
for (weight in w) {
  dpc = c(dpc, g_2w_corr(weight, rho=0.79))
}
plot(w, dpc)


p_rho_2 <- function(p, rhos=seq(0, 1, by=0.01), n=2, alpha=0.025) {
  C.alpha = qnorm(1-1/n*alpha) - qnorm(1-alpha)
  p_rhos =c()
  for (rho in rhos) {
    p_rhos = c(p_rhos, 1/(1+rho)*qnorm(p*(1+rho)*sqrt(2*pi*(1-rho^2))) - qnorm(1-p))
  }
  plot(rhos, p_rhos)
  abline(h=C.alpha)
}
