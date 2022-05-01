library(mvtnorm)
library(ggplot2)
library(future.apply)

################################### 2 hypotheses
w1 <- seq(0, 1, 0.001)
w2 <- seq(0, 1, 0.001)
w <- expand.grid(w2, w1)
w <- subset(w, Var1 + Var2 == 1)
rho <- seq(0, 1, 0.01)
mp <- seq(0.3, 0.9, 0.1)
arg <- expand.grid(mp, rho)
n <- nrow(w)
arg <- do.call("rbind", replicate(n, arg, simplify = FALSE))
arg <- arg[order(arg$Var1, arg$Var2),]
w_list <- list()
for(i in 1:(length(rho) * length(mp))) {
  w_list[[i]] <- w
}
scen <- cbind(arg, do.call(rbind, w_list))
colnames(scen) <- c("mp", "rho", "w1", "w2")
rownames(scen) <- NULL
alpha <- 0.025
scen <- list(scen, alpha)

search_2_eq <- function(x, scen) {
  alpha <- scen[[2]]
  mp <- rep(scen[[1]][x, 1], 2)
  sigma <- matrix(scen[[1]][x, 2], nrow = length(mp), ncol = length(mp))
  diag(sigma) <- 1
  temp <- as.numeric(scen[[1]][x, 3:4])
  y <- 1 - pmvnorm(upper = qnorm(1 - alpha * temp),
                   mean = qnorm(1 - alpha) - qnorm(1 - mp),
                   sigma = sigma,
                   algorithm = GenzBretz(maxpts = 1e7, abseps = 1e-7, releps = 1e-7)
  )
  y <- matrix(c(mp, scen[[1]][x, 2], temp, as.numeric(y)), nrow = 1)
  colnames(y) <- c("mp1", "mp2", "rho", "w1", "w2", "dp")
  return(y)
}

search_2_uneq <- function(x, scen) {
  alpha <- scen[[2]]
  mp <- c(0.9, scen[[1]][x, 1])
  sigma <- matrix(scen[[1]][x, 2], nrow = length(mp), ncol = length(mp))
  diag(sigma) <- 1
  temp <- as.numeric(scen[[1]][x, 3:4])
  y <- 1 - pmvnorm(upper = qnorm(1 - alpha * temp),
                   mean = qnorm(1 - alpha) - qnorm(1 - mp),
                   sigma = sigma,
                   algorithm = GenzBretz(maxpts = 1e7, abseps = 1e-7, releps = 1e-7)
  )
  y <- matrix(c(mp, scen[[1]][x, 2], temp, as.numeric(y)), nrow = 1)
  colnames(y) <- c("mp1", "mp2", "rho", "w1", "w2", "dp")
  return(y)
}

# Parallel
library(future.apply)
plan(multisession)
# 2 hypotheses with equal marginal power
start <- proc.time()
result <- future_lapply(1:2, FUN = search_2_eq, future.seed = TRUE,
                        future.packages = c("mvtnorm"), scen = scen)
(proc.time() - start)[3]
data <- as.data.frame(do.call(rbind, result))
write.csv(data, file="dp_2_eq.csv")

# 2 hypotheses with unequal marginal power
start <- proc.time()
result <- future_lapply(1:2, FUN = search_2_eq, future.seed = TRUE,
                        future.packages = c("mvtnorm"), scen = scen)
(proc.time() - start)[3]
data <- as.data.frame(do.call(rbind, result))
write.csv(data, file="dp_2_uneq.csv")

# Example call
setwd("C:/Users/dxi1/OneDrive - Gilead Sciences/Paper/Optimal graph/code")
data <- read.csv("dp_2_eq.csv")
mp_1 <- 0.9
mp_2 <- 0.9
cr <- 0.8
temp <- subset(data, mp1 == mp_1 & mp2 == mp_2 & rho == cr)
temp[order(-temp$dp), ][1, ]
subset(temp, w1 %in% c(0, 0.5, 1))


# Figure 1
setwd("C:/Users/dxi1/OneDrive - Gilead Sciences/Paper/Optimal graph/code")
data <- read.csv("dp_2_eq.csv")
temp <- subset(data, mp1 == 0.9 & mp2 == 0.9 & rho %in% c(0.3, 0.6, 0.82, 0.9))
ggplot(data = temp, aes(x = w1, y = dp)) +
  geom_line(size = 1) +
  facet_wrap(~ rho, labeller = label_both) +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(breaks = seq(0.85, 1, 0.05), limits = c(0.875, 0.975)) +
  xlab("w1") +
  ylab("Disjunctive power") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        strip.text.x = element_text(size = 15))

# Figure 2
setwd("C:/Users/dxi1/OneDrive - Gilead Sciences/Paper/Optimal graph/code")
data <- read.csv("dp_2_eq.csv")
mp <- seq(0.5, 0.9, 0.1)
data_plot <- expand.grid(rho, mp)
data_plot <- data.frame(data_plot, maxdp = NA)
colnames(data_plot) <- c("rho", "mp", "maxdp")
for (i in 1:nrow(data_plot)) {
  temp <- subset(data, mp1 == data_plot$mp[i] & mp2 == data_plot$mp[i] &
                   abs(rho - data_plot$rho[i]) < 1e-7)
  data_plot[i, 3] <- max(temp$dp)
}
ggplot(data = data_plot, aes(x = rho, y = maxdp, group = as.factor(mp))) +
  geom_line(aes(linetype = as.factor(mp)), size = 1) +
  scale_linetype_manual(values=c("solid", "longdash", "dashed", "dotdash", "dotted"))+
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(breaks = seq(0.5, 1, 0.1), limits = c(0.5, 1)) +
  xlab("Correlation") +
  ylab("Optimal disjunctive power") +
  labs(linetype = "Marginal power") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = "top")

# Table 3
setwd("C:/Users/dxi1/OneDrive - Gilead Sciences/Paper/Optimal graph/code")
data <- read.csv("dp_2_eq.csv")
mp <- c(0.9, 0.8, 0.7)
sigma <- c(0.9, 0.8, 0.7, 0.5)
tab <- NULL
for (i in 1:length(mp)) {
  for (j in 1:length(sigma)) {
    temp <- subset(data, mp1 == mp[i] & mp2 == mp[i] & abs(rho - sigma[j]) < 1e-7)
    x <- temp[order(-temp$dp), ][1, ]
    x[length(x)] <- round(100 * x[length(x)], 3)
    tab <- rbind(tab, x)
    if (x$w1 != 0.5) {
      y <- temp[order(-temp$dp), ][2, ]
      if (temp[order(-temp$dp), ][1, ]$dp == y$dp) {
        y[length(y)] <- round(100 * y[length(y)], 3)
        tab <- rbind(tab, y)
      }
      z <- subset(temp, w1 == 0.5)
      z[length(z)] <- round(100 * z[length(z)], 3)
      tab <- rbind(tab, z)
    }
  }
}
tab

# Table 4
setwd("C:/Users/dxi1/OneDrive - Gilead Sciences/Paper/Optimal graph/code")
data <- read.csv("dp_2_uneq.csv")
mp <- c(0.8, 0.6, 0.3)
sigma <- c(0.9, 0.8, 0.7, 0.5)
tab <- NULL
for (i in 1:length(mp)) {
  for (j in 1:length(sigma)) {
    temp <- subset(data, mp1 == 0.9 & mp2 == mp[i] & abs(rho - sigma[j]) < 1e-7)
    x <- temp[order(-temp$dp), ][1, ]
    x[length(x)] <- round(100 * x[length(x)], 3)
    tab <- rbind(tab, x)
  }
}
tab

bench <- NULL
alpha <- 0.025
for (i in 1:length(sigma)) {
  cr <- matrix(sigma[i], nrow = 2, ncol = 2)
  diag(cr) <- 1
  y <- 1 - pmvnorm(upper = qnorm(1 - alpha * c(0.597, 0.403)),
                   mean = qnorm(1 - alpha) - qnorm(1 - c(0.9, 0.8)),
                   sigma = cr,
                   algorithm = GenzBretz(maxpts = 1e7, abseps = 1e-7, releps = 1e-7)
  )
  bench <- rbind(bench, c(0.8, sigma[i], 0.597, 0.403, round(100 * as.numeric(y), 3)))
}
for (i in 1:length(sigma)) {
  cr <- matrix(sigma[i], nrow = 2, ncol = 2)
  diag(cr) <- 1
  y <- 1 - pmvnorm(upper = qnorm(1 - alpha * c(0.755, 0.245)),
                   mean = qnorm(1 - alpha) - qnorm(1 - c(0.9, 0.6)),
                   sigma = cr,
                   algorithm = GenzBretz(maxpts = 1e7, abseps = 1e-7, releps = 1e-7)
  )
  bench <- rbind(bench, c(0.6, sigma[i], 0.755, 0.245, round(100 * as.numeric(y), 3)))
}
for (i in 1:length(sigma)) {
  cr <- matrix(sigma[i], nrow = 2, ncol = 2)
  diag(cr) <- 1
  y <- 1 - pmvnorm(upper = qnorm(1 - alpha * c(0.957, 0.043)),
                   mean = qnorm(1 - alpha) - qnorm(1 - c(0.9, 0.3)),
                   sigma = cr,
                   algorithm = GenzBretz(maxpts = 1e7, abseps = 1e-7, releps = 1e-7)
  )
  bench <- rbind(bench, c(0.3, sigma[i], 0.957, 0.043, round(100 * as.numeric(y), 3)))
}
bench

################################### 3 hypotheses

search_3_eq <- function(x, scen) {
  alpha <- scen[[2]]
  mp <- rep(scen[[3]], 3)
  sigma <- matrix(scen[[4]], nrow = length(mp), ncol = length(mp))
  diag(sigma) <- 1
  temp <- as.numeric(scen[[1]][x, ])
  y <- 1 - pmvnorm(upper = qnorm(1 - alpha * temp),
                   mean = qnorm(1 - alpha) - qnorm(1 - mp),
                   sigma = sigma,
                   algorithm = GenzBretz(maxpts = 1e7, abseps = 1e-7, releps = 1e-7)
  )
  y <- matrix(c(mp, scen[[4]], temp, as.numeric(y)), nrow = 1)
  colnames(y) <- c("mp1", "mp2", "mp3", "rho", "w1", "w2", "w3", "dp")
  return(y)
}

w1 <- seq(0, 1, 0.001)
w2 <- seq(0, 1, 0.001)
w3 <- seq(0, 1, 0.001)
w <- expand.grid(w2, w1)
w <- subset(w, Var1 + Var2 <= 1)
arg <- rep(NA, 3)
for (i in 1:nrow(w)) {
  temp <- as.numeric(w[i, ])
  id <- as.vector(which(sum(temp) + w3 == 1))
  if (length(id) > 0) {
    for (j in 1:length(id)) {
      arg <- rbind(arg, c(temp, w3[id[j]]))
    }
  }
}
arg <- arg[-1, ]
write.csv(arg, file="w_3.csv")
alpha <- 0.025
mp <- 0.9
sigma <- 0.9
scen <- list(arg, alpha, mp, sigma)

plan(multisession)
# 3 hypotheses with equal marginal power
start <- proc.time()
result <- future_lapply(1:nrow(scen[[1]]), FUN = search_3_eq, future.seed = TRUE,
                        future.packages = c("mvtnorm"), scen = scen)
(proc.time() - start)[3]
data <- as.data.frame(do.call(rbind, result))

# Table 5
alpha <- 0.025
setwd("C:/Users/dxi1/OneDrive - Gilead Sciences/Paper/Optimal graph/code")
data <- read.csv("dp_3_eq_0.9.csv")[, -1]
sigma <- c(0.9, 0.8, 0.7, 0.5)
tab_0.9 <- NULL
for (j in 1:length(sigma)) {
  temp <- subset(data, abs(rho - sigma[j]) < 1e-7)
  x <- temp[order(-temp$dp), ][1, ]
  x[length(x)] <- round(100 * x[length(x)], 3)
  if (x$w1 %in% c(0, 0.5, 1) |
      x$w2 %in% c(0, 0.5, 1) |
      x$w3 %in% c(0, 0.5, 1)) {
    tab_0.9 <- rbind(tab_0.9, x)
    y <- temp[order(-temp$dp), ][2, ]
    if (temp[order(-temp$dp), ][1, ]$dp == y$dp) {
      y[length(y)] <- round(100 * y[length(y)], 3)
      tab_0.9 <- rbind(tab_0.9, y)
    }
    y <- temp[order(-temp$dp), ][3, ]
    if (temp[order(-temp$dp), ][1, ]$dp == y$dp) {
      y[length(y)] <- round(100 * y[length(y)], 3)
      tab_0.9 <- rbind(tab_0.9, y)
    }
    cr <- matrix(sigma[j], nrow = 3, ncol = 3)
    diag(cr) <- 1
    y <- 1 - pmvnorm(upper = qnorm(1 - alpha * rep(1/3, 3)),
                     mean = qnorm(1 - alpha) - qnorm(1 - rep(0.9, 3)),
                     sigma = cr,
                     algorithm = GenzBretz(maxpts = 1e7, abseps = 1e-7, releps = 1e-7)
    )
    z <- c(rep(0.9, 3), sigma[j], round(rep(1/3, 3), 3), round(100 * y, 3))
    tab_0.9 <- rbind(tab_0.9, z)
  } else if (x$w1 %in% c(0.333, 0.334) &
             x$w2 %in% c(0.333, 0.334) &
             x$w3 %in% c(0.333, 0.334)) {
    cr <- matrix(sigma[j], nrow = 3, ncol = 3)
    diag(cr) <- 1
    y <- 1 - pmvnorm(upper = qnorm(1 - alpha * rep(1/3, 3)),
                     mean = qnorm(1 - alpha) - qnorm(1 - rep(0.9, 3)),
                     sigma = cr,
                     algorithm = GenzBretz(maxpts = 1e7, abseps = 1e-7, releps = 1e-7)
    )
    z <- c(rep(0.9, 3), sigma[j], round(rep(1/3, 3), 3), round(100 * y, 3))
    tab_0.9 <- rbind(tab_0.9, z)
  }
}
tab_0.9

data <- read.csv("dp_3_eq_0.8.csv")[, -1]
sigma <- c(0.9, 0.79, 0.7, 0.5)
tab_0.8 <- NULL
for (j in 1:length(sigma)) {
  temp <- subset(data, abs(rho - sigma[j]) < 1e-7)
  x <- temp[order(-temp$dp), ][1, ]
  x[length(x)] <- round(100 * x[length(x)], 3)
  if (x$w1 %in% c(0, 0.5, 1) |
      x$w2 %in% c(0, 0.5, 1) |
      x$w3 %in% c(0, 0.5, 1)) {
    tab_0.8 <- rbind(tab_0.8, x)
    y <- temp[order(-temp$dp), ][2, ]
    if (temp[order(-temp$dp), ][1, ]$dp == y$dp) {
      y[length(y)] <- round(100 * y[length(y)], 3)
      tab_0.8 <- rbind(tab_0.8, y)
    }
    y <- temp[order(-temp$dp), ][3, ]
    if (temp[order(-temp$dp), ][1, ]$dp == y$dp) {
      y[length(y)] <- round(100 * y[length(y)], 3)
      tab_0.8 <- rbind(tab_0.8, y)
    }
    cr <- matrix(sigma[j], nrow = 3, ncol = 3)
    diag(cr) <- 1
    y <- 1 - pmvnorm(upper = qnorm(1 - alpha * rep(1/3, 3)),
                     mean = qnorm(1 - alpha) - qnorm(1 - rep(0.8, 3)),
                     sigma = cr,
                     algorithm = GenzBretz(maxpts = 1e7, abseps = 1e-7, releps = 1e-7)
    )
    z <- c(rep(0.8, 3), sigma[j], round(rep(1/3, 3), 3), round(100 * y, 3))
    tab_0.8 <- rbind(tab_0.8, z)
  } else if (x$w1 %in% c(0.333, 0.334) &
             x$w2 %in% c(0.333, 0.334) &
             x$w3 %in% c(0.333, 0.334)) {
    cr <- matrix(sigma[j], nrow = 3, ncol = 3)
    diag(cr) <- 1
    y <- 1 - pmvnorm(upper = qnorm(1 - alpha * rep(1/3, 3)),
                     mean = qnorm(1 - alpha) - qnorm(1 - rep(0.8, 3)),
                     sigma = cr,
                     algorithm = GenzBretz(maxpts = 1e7, abseps = 1e-7, releps = 1e-7)
    )
    z <- c(rep(0.8, 3), sigma[j], round(rep(1/3, 3), 3), round(100 * y, 3))
    tab_0.8 <- rbind(tab_0.8, z)
  }
}
tab_0.8


data <- read.csv("dp_3_eq_0.7.csv")[, -1]
sigma <- c(0.9, 0.77, 0.7, 0.5)
tab_0.7 <- NULL
for (j in 1:length(sigma)) {
  temp <- subset(data, abs(rho - sigma[j]) < 1e-7)
  x <- temp[order(-temp$dp), ][1, ]
  x[length(x)] <- round(100 * x[length(x)], 3)
  if (x$w1 %in% c(0, 0.5, 1) |
      x$w2 %in% c(0, 0.5, 1) |
      x$w3 %in% c(0, 0.5, 1)) {
    tab_0.7 <- rbind(tab_0.7, x)
    y <- temp[order(-temp$dp), ][2, ]
    if (temp[order(-temp$dp), ][1, ]$dp == y$dp) {
      y[length(y)] <- round(100 * y[length(y)], 3)
      tab_0.7 <- rbind(tab_0.7, y)
    }
    y <- temp[order(-temp$dp), ][3, ]
    if (temp[order(-temp$dp), ][1, ]$dp == y$dp) {
      y[length(y)] <- round(100 * y[length(y)], 3)
      tab_0.7 <- rbind(tab_0.7, y)
    }
    cr <- matrix(sigma[j], nrow = 3, ncol = 3)
    diag(cr) <- 1
    y <- 1 - pmvnorm(upper = qnorm(1 - alpha * rep(1/3, 3)),
                     mean = qnorm(1 - alpha) - qnorm(1 - rep(0.7, 3)),
                     sigma = cr,
                     algorithm = GenzBretz(maxpts = 1e7, abseps = 1e-7, releps = 1e-7)
    )
    z <- c(rep(0.7, 3), sigma[j], round(rep(1/3, 3), 3), round(100 * y, 3))
    tab_0.7 <- rbind(tab_0.7, z)
  } else if (x$w1 %in% c(0.333, 0.334) &
             x$w2 %in% c(0.333, 0.334) &
             x$w3 %in% c(0.333, 0.334)) {
    cr <- matrix(sigma[j], nrow = 3, ncol = 3)
    diag(cr) <- 1
    y <- 1 - pmvnorm(upper = qnorm(1 - alpha * rep(1/3, 3)),
                     mean = qnorm(1 - alpha) - qnorm(1 - rep(0.7, 3)),
                     sigma = cr,
                     algorithm = GenzBretz(maxpts = 1e7, abseps = 1e-7, releps = 1e-7)
    )
    z <- c(rep(0.7, 3), sigma[j], round(rep(1/3, 3), 3), round(100 * y, 3))
    tab_0.7 <- rbind(tab_0.7, z)
  }
}
tab_0.7
