# Time assessment


#' Check if it is necessary to keep independent scenario or
#' dependent case with rho = 0
#' Conclusion: independent case is much faster than dependent
#' case with rho = 0
time.idp = c()
iter.idp = c()
for (i in 1:100) {
  start = Sys.time()
  res = optim_w_dp(alpha = 0.025, mp = c(0.9, 0.8, 0.7))
  end = Sys.time()
  time.idp = c(time.idp, end - start)
  iter.idp = c(iter.idp, res$optima$iterations)
}

time.dp = c()
iter.dp = c()
for (i in 1:100) {
  start = Sys.time()
  optim_w_dp(alpha = 0.025, mp = c(0.9, 0.8, 0.7), rho = 0)
  end = Sys.time()
  time.dp = c(time.dp, end - start)
  iter.dp = c(iter.dp, res$optima$iterations)
}
