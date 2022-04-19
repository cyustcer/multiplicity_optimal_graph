########################
### Independent Case ###
########################

power <- function(alpha, eta=2.5) {
  1 - pnorm(qnorm(1-alpha) - eta)
}


solve_eta <- function(alpha, p) {
  equation <- function(alpha, eta) {
    power(alpha, eta) - p
  }
  root = uniroot(equation, alpha=alpha, lower=-1000, upper=1000)$root
  root
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
