
is.connected <- function(G1_2, G1_3, G2_1, G2_3, G3_1, G3_2) {
  if ((G1_2 > 0 | (G1_3 > 0 & G3_2 > 0)) &
      (G1_3 > 0 | (G1_2 > 0 & G2_3 > 0)) &
      (G2_3 > 0 | (G2_1 > 0 & G1_3 > 0)) &
      (G2_1 > 0 | (G2_3 > 0 & G3_1 > 0)) &
      (G3_1 > 0 | (G3_2 > 0 & G2_1 > 0)) &
      (G3_2 > 0 | (G3_1 > 0 & G1_2 > 0))) {
    return(TRUE)
  }
  return(FALSE)
}


load.sim.grid <- function(filedir){
  sim.grid.complete <- NULL
  files = list.files(filedir)
  for (file in files) {
    filepath = paste0(filedir, file)
    load(filepath)
    sim.grid.complete <- rbind(sim.grid.complete, data.sim)
  }
  
  return(sim.grid.complete)
}


is.valid.bound <- function(bound, is.upper=TRUE) {
  if (is.upper) {
    return(max(bound) <= 1)
  }
  else {
    return(sum(bound) <= 1)
  }
}


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


powerset <- function(n){
  l <- matrix(0, nrow = 2^n, ncol = n)
  counter = 1
  for(x in n:1){
    for(subset in 1:counter){
      counter <- counter + 1
      l[counter, ] <- l[subset, ]
      l[counter, x] <- 1
    }
  }
  return(l)
}


create_initial_ws_dx <- function(constraint.w) {
  n <- length(constraint.w)
  remain.ws <- 1 - sum(constraint.w, na.rm = T)
  loc <- is.na(constraint.w)
  n1 <- sum(loc)
  initial.ws <- powerset(n1)[-1, ]
  initial.ws <- initial.ws / rowSums(initial.ws) * remain.ws
  out <- matrix(NA, nrow = nrow(initial.ws), ncol = n)
  out[, loc] <- initial.ws
  out[, !loc] <- matrix(rep(constraint.w[!loc], nrow(initial.ws)), nrow = nrow(initial.ws))
  colnames(out) <- paste0("w", 1:n)
  return(as.data.frame(out))
}


adjust_initial_ws <- function(initial.ws, lb, ub) {
  n = nrow(initial.ws)
  to_remove = c()
  for (i in 1:n) {
    if (!(all(initial.ws[i, ] <= ub))) {
      initial.ws[i, initial.ws[i, ] > ub] = ub[initial.ws[i, ] > ub]
      remain.ws <- 1 - sum(initial.ws[i, initial.ws[i, ] == ub])
      to_modify = initial.ws[i, ] != ub & initial.ws[i, ] != 0
      initial.ws[i, to_modify] = remain.ws / sum(to_modify)
    }
    if (!(all(initial.ws[i, ] >= lb))) {
      if (any(initial.ws[i, initial.ws[i, ] < lb] == 0)) {
        to_remove = c(to_remove, i)
      }
      else {
        initial.ws[i, initial.ws[i, ] < lb] = lb[initial.ws[i, ] < lb]
        remain.ws <- 1 - sum(initial.ws[i, initial.ws[i, ] == lb])
        to_modify = initial.ws[i, ] != lb & initial.ws[i, ] != 0
        initial.ws[i, to_modify] = remain.ws / sum(to_modify)
      }
    }
  }
  if (length(to_remove) > 0) {
    initial.ws = initial.ws[-to_remove, ]
  }
  initial.ws = initial.ws[rowSums(initial.ws) == 1, ]
  return(initial.ws)
}

