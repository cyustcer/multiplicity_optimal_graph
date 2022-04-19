# w optimization
source(here::here("src/optimization/disjunctive_power.R"))

# =================================================
# =============== Independent cases =============== 
# =================================================

# Optimization set up
opts <- list("algorithm" = "NLOPT_LD_SLSQP",
             "xtol_rel" = 0,
             "maxeval" = 10000)

# ========= Example One =========
ALPHA = 0.025
MP = c(0.9, 0.9)
initial_w = c(0.999, 0.001)
LB = c(0, 0)
UB = c(1, 1)

res <- nloptr(x0=initial_w,
              eval_f=loss,
              eval_grad_f = loss_grad,
              alpha=ALPHA,
              mp=MP,
              lb=LB, ub=UB,
              eval_g_eq = eqn,
              eval_jac_g_eq = eqn_grad,
              opts=opts)


# ========= Example Two =========
ALPHA = 0.025
MP = c(0.8, 0.8, 0.8)
# initial_w = c(0.998, 0.001, 0.001)
initial_w = rep(1, length(MP))/length(MP)
LB = c(0, 0, 0)
UB = c(1, 1, 1)


res <- nloptr(x0=initial_w,
              eval_f=loss,
              eval_grad_f = loss_grad,
              alpha=ALPHA,
              mp=MP,
              lb=LB, ub=UB,
              eval_g_eq = eqn,
              eval_jac_g_eq = eqn_grad,
              opts=opts)


# Other optimization try
# res <- solnp(initial_w, loss, eqfun=eqn, eqB=c(1), LB=LB, UB=UB)

# ========= Example Three =========
ALPHA = 0.025
MP = c(0.9, 0.8, 0.1)
# initial_w = c(0.998, 0.001, 0.001)
initial_w = rep(1, length(MP))/length(MP)
LB = c(0, 0, 0)
UB = c(1, 1, 1)


res <- nloptr(x0=initial_w,
              eval_f=loss,
              eval_grad_f = loss_grad,
              alpha=ALPHA,
              mp=MP,
              lb=LB, ub=UB,
              eval_g_eq = eqn,
              eval_jac_g_eq = eqn_grad,
              opts=opts)


# ================================================
# =============== Correlated cases =============== 
# ================================================


# Optimization set up
opts <- list("algorithm" = "NLOPT_LD_SLSQP",
             "xtol_rel" = 0,
             "maxeval" = 10000)

# ========= Example One =========
ALPHA = 0.025
MP = c(0.9, 0.9)
initial_w = rep(1, length(MP))/length(MP)
RHO = 0.5
LB = c(0, 0)
UB = c(1, 1)

res <- nloptr(x0=initial_w,
              eval_f=loss_corr,
              eval_grad_f = loss_grad_corr,
              alpha=ALPHA,
              mp=MP,
              rho=RHO,
              lb=LB, ub=UB,
              eval_g_eq = eqn_corr,
              eval_jac_g_eq = eqn_grad_corr,
              opts=opts)

# ========= Example Two =========
ALPHA = 0.025
MP = c(0.8, 0.8, 0.8)
initial_w = rep(1, length(MP))/length(MP)
RHO = 2/3
LB = c(0, 0, 0)
UB = c(1, 1, 1)


res <- nloptr(x0=initial_w,
              eval_f=loss_corr,
              eval_grad_f = loss_grad_corr,
              alpha=ALPHA,
              mp=MP,
              rho=RHO,
              lb=LB, ub=UB,
              eval_g_eq = eqn_corr,
              eval_jac_g_eq = eqn_grad_corr,
              opts=opts)

# ========= Example Three =========
ALPHA = 0.025
MP = c(0.9, 0.5, 0.3)
initial_w = c(0.998, 0.001, 0.001)
# initial_w = rep(1, length(MP))/length(MP)
# initial_w = c(1/3, 1/3, 1/3)
RHO = 0.5
LB = c(1e-6, 1e-6, 1e-6)
UB = c(1, 1, 1)


res <- nloptr(x0=initial_w,
              eval_f=loss_corr,
              eval_grad_f = loss_grad_corr,
              alpha=ALPHA,
              mp=MP,
              rho=RHO,
              lb=LB, ub=UB,
              eval_g_eq = eqn_corr,
              eval_jac_g_eq = eqn_grad_corr,
              opts=opts)

