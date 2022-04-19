source(here::here('src/optimization/disjunctive_power.R'))
source(here::here("src/optimization/conjunctive_power.R"))

####################
### Optimization ###
####################

alpha = 0.025

# Example 1
mp.1 = c(0.3, 0.7, 0.9)
(w.optim.1 = optim_w(alpha, mp.1)$w)
(G.optim.1 = optim_G(alpha, mp.1, w.optim.1))

# Example 2
mp.2 = c(0.7, 0.8, 0.9, 0.9)
(w.optim.2 = optim_w(alpha, mp.2)$w)
(G.optim.2 = optim_G(alpha, mp.2, w.optim.2))

# Example 3
mp.3 = rep(0.9, 5)
(w.optim.3 = optim_w(alpha, mp.3)$w)
(G.optim.3 = optim_G(alpha, mp.3, w.optim.3))


###################################
### Optimization w. constraints ###
###################################

alpha = 0.025

# Example 1 w. constraints
mp.1 = c(0.3, 0.7, 0.9)
w.constraint.1 = c(0, NA, NA)
G.constraint.1 = matrix(NA, length(mp.1), length(mp.1))
diag(G.constraint.1) = 0  # diagonal needs to be 0
G.constraint.1[1, 3] = 1  # first hypothesis only pass to 3rd hypothesis
(w.optim.constraint.1 = optim_w(alpha, mp.1, constraint.w=w.constraint.1)$w)
(G.optim.constraint.1 = optim_G_constrain(alpha, mp.1, w.optim.constraint.1, constraint.G=G.constraint.1))

# Example 4
mp.4 = c(0.9, 0.9, 0.9, 0.9)
w.constraint.4 = c(NA, NA, 0, 0)
G.constraint.4 = matrix(NA, length(mp.4), length(mp.4))
diag(G.constraint.4) = 0
G.constraint.4[1, 4] = 0
G.constraint.4[2, 3] = 0
G.constraint.4[3, 2] = 1
G.constraint.4[4, 1] = 1
(w.optim.constraint.4 = optim_w(alpha, mp.4, constraint.w=w.constraint.4)$w)
(G.optim.constraint.4 = optim_G_constrain(alpha, mp.4, w.optim.constraint.4, constraint.G=G.constraint.4))

# Example 5
mp.5 = c(0.9, 0.8, 0.6, 0.5)
w.constraint.5 = c(NA, NA, 0, 0)
G.constraint.5 = matrix(NA, length(mp.5), length(mp.5))
diag(G.constraint.5) = 0
G.constraint.5[1, 4] = 0
G.constraint.5[2, 3] = 0
G.constraint.5[3, 2] = 1
G.constraint.5[4, 1] = 1
(w.optim.constraint.5 = optim_w(alpha, mp.5, constraint.w=w.constraint.5)$w)
(G.optim.constraint.5 = optim_G_constrain(alpha, mp.5, w.optim.constraint.5, constraint.G=G.constraint.5))

# Example 6
mp.6 = c(0.9, 0.8, 0.6, 0.3)
w.constraint.6 = c(NA, NA, 0, 0)
G.constraint.6 = matrix(NA, length(mp.6), length(mp.6))
diag(G.constraint.6) = 0
G.constraint.6[1, 4] = 0
G.constraint.6[2, 3] = 0
G.constraint.6[3, 2] = 1
G.constraint.6[4, 1] = 1
(w.optim.constraint.6 = optim_w(alpha, mp.6, constraint.w=w.constraint.6)$w)
(G.optim.constraint.6 = optim_G_constrain(alpha, mp.6, w.optim.constraint.6, constraint.G=G.constraint.6))
