source('src/objectives/independent_disjunctive_power.R')
source('src/objectives/correlated_disjunctive_power.R')
source('src/objectives/correlated_conjunctive_power.R')


create_initial_ws <- function(constraint.w) {
  n = length(constraint.w)
  remain.ws = 1 - sum(constraint.w, na.rm=T)
  cand.ws <- NULL
  for (i in 1:n) {
    if (is.na(constraint.w[i])) {
      cand = list(c(0, 1))
      names(cand) = paste0("w", i)
      cand.ws <- c(cand.ws, cand)
    }
    else {
      cand = list(c(constraint.w[i]))
    }
  }
  initial.ws <- cand.ws %>% expand.grid()
  initial.ws = initial.ws[rowSums(initial.ws) != 0, , drop=F]
  initial.ws = initial.ws / rowSums(initial.ws) * remain.ws
  row.names(initial.ws) <- NULL
  for (i in 1:n) {
    if (!is.na(constraint.w[i])) {
      initial.ws[paste0("w", i)] = constraint.w[i]
    }
  }
  initial.ws %>% select(sort(names(initial.ws)))
}


opts <- list("algorithm" = "NLOPT_LD_SLSQP",
             "xtol_rel" = 1e-3,
             "maxeval" = 1e3)


optim_w_dp <- function(alpha, mp, rho=NULL,
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
  if (optim_opts$algorithm == "NLOPT_LN_COBYLA") {
    if (is.null(rho)) {
      res <- nloptr(x0 = initial_w,
                    eval_f = loss_dp,
                    eval_grad_f = loss_dp_grad,
                    alpha = alpha, mp = mp.p,
                    lb = lb, ub = ub, min.w = min.w,
                    eval_g_ineq = ineqn,
                    opts=optim_opts)
    }
    else {
      res <- nloptr(x0=initial_w,
                    eval_f=loss_dpc,
                    eval_grad_f = loss_dpc_grad,
                    alpha=alpha, mp=mp.p, rho=rho,
                    lb=lb, ub=ub, min.w=min.w,
                    eval_g_ineq = ineqn_corr,
                    opts=optim_opts)
    }
  }
  else {
    if (is.null(rho)) {
      res <- nloptr(x0 = initial_w,
                    eval_f = loss_dp,
                    eval_grad_f = loss_dp_grad,
                    alpha = alpha, mp = mp.p,
                    lb = lb, ub = ub, min.w = min.w,
                    eval_g_eq = eqn,
                    eval_jac_g_eq = eqn_grad,
                    opts=optim_opts)
    }
    else {
      res <- nloptr(x0=initial_w,
                    eval_f=loss_dpc,
                    eval_grad_f = loss_dpc_grad,
                    alpha=alpha, mp=mp.p, rho=rho,
                    lb=lb, ub=ub, min.w=min.w,
                    eval_g_eq = eqn_corr,
                    eval_jac_g_eq = eqn_corr_grad,
                    opts=optim_opts)
    }
  }
  w.p = res$solution
  w = rep(0, length(mp))
  w[constraint.w != 0 | is.na(constraint.w)] = w.p
  return(list(w = w,
              optima = res))
}

# Example
#' Miwa always return equal split
#' Miwa does not work for more than 5-dimension
#' GenzBretz time variance is large
#' GenzBretz does not return equal split
# res <- optim_w_dp(alpha=0.025, mp=rep(0.9, 5), rho=0.6)

go_optim_w_dp <- function(alpha, mp, rho=NULL,
                          constraint.w=NULL, lb=NULL, ub=NULL,
                          min.w=1e-12, optim_opts=opts) {
  n = length(mp)
  if (is.null(constraint.w)) {
    constraint.w = rep(NA, n)
  }
  initial.ws = create_initial_ws(constraint.w)
  N = nrow(initial.ws)
  optimas <- data.frame(matrix(ncol = 2*n+1, nrow = 0))
  colnames(optimas) <- c(paste0('initial ', colnames(initial.ws)), 
                         colnames(initial.ws), 'optimal_value')
  for (i in 1:N) {
    initial.w = as.matrix(initial.ws)[i, ]
    # print(initial.w)
    res = optim_w_dp(alpha = alpha, mp = mp, rho = rho,
                     initial_w = initial.w, 
                     constraint.w = constraint.w,
                     lb = lb, ub = ub, min.w = min.w,
                     optim_opts = optim_opts)
    # optimal.w = round(res$w, digits=pcs.digits)
    optimal.w = res$w
    optimal.value = -res$optima$objective
    optimas[i, ] <- c(initial.w, optimal.w, optimal.value)
  }
  optimas %>% arrange(desc(optimal_value))
}


optim_w_cp <- function(alpha, mp, rho=NULL,
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
                eval_f=loss_cpc,
                eval_grad_f = loss_cpc_grad,
                alpha=alpha, mp=mp.p, rho=rho,
                lb=lb, ub=ub, min.w=min.w,
                eval_g_eq = eqn_corr,
                eval_jac_g_eq = eqn_corr_grad,
                opts=opts)
  w.p = res$solution
  w = rep(0, length(mp))
  w[constraint.w != 0 | is.na(constraint.w)] = w.p
  return(list(w = w,
              optima = res))
}
