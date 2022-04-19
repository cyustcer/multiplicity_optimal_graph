library(tidyverse)
library(data.table)


gen.w.grid <- function(n.hypo, w.step=0.1, w.range=c(0.33, 0.34)) {
  w.cand <- NULL
  for (i in 1:n.hypo) {
    tmp = list(seq(w.range[1], w.range[2], w.step))
    names(tmp) = paste0("w", i)
    w.cand <- c(w.cand, tmp)
  }
  w.grid = w.cand %>% expand.grid()
  w.grid = w.grid[rowSums(w.grid) == 1, ]
  w.grid = as.matrix(w.grid)
  row.names(w.grid) <- NULL
  return(w.grid)
}


gen.G.grid <- function(n.hypo, G.step=0.25) {
  G.grid <- c(0)  # Placeholder for merge
  for (i in 1:n.hypo) {
    G.cand <- NULL
    for (j in 1:n.hypo) {
      if (i != j) {
        tmp = list(seq(0, 1, G.step))
        names(tmp) = paste0("G", i, "_", j)
        G.cand <- c(G.cand, tmp)
      }
      if (i == j) {
        tmp = list(c(0))
        names(tmp) = paste0("G", i, "_", j)
        G.cand <- c(G.cand, tmp)
      }
    }
    G.tmp <- G.cand %>% expand.grid()
    G.tmp <- G.tmp[rowSums(G.tmp) == 1, ]
    G.grid <- merge(G.grid, G.tmp)
  }
  G.grid = G.grid[, -1]
  G.grid = as.matrix(G.grid)
  return(G.grid)
}


gen.data.grid <- function(n.hypo, w.step=0.1, w.range=c(0.33, 0.34), G.step=0.25, save.dir=here::here('data/grids')) {
  grid.name.w = paste0("n.hypo_", n.hypo, "_w.step_", w.step, "_w.range_", paste(w.range, collapse="_"), ".rds")
  grid.name.G = paste0("n.hypo_", n.hypo, "_G.step_", G.step, ".rds")
  grid.name.data = paste0("n.hypo_", n.hypo, "_w.step_", w.step,
                          "_w.range_", paste(w.range, collapse="_"),
                          "_G.step_", G.step, ".rds")
  path.grid.data <- here::here(save.dir, grid.name.data)
  path.grid.w <- here::here(save.dir, grid.name.w)
  path.grid.G <- here::here(save.dir, grid.name.G)
  if (!file.exists(path.grid.data)) {
    # generate or load w grids
    if (!file.exists(path.grid.w)) {
      w.grid = gen.w.grid(n.hypo, w.step=w.step, w.range=w.range)
      saveRDS(w.grid, path.grid.w)
    } else {
      w.grid = readRDS(path.grid.w)
    }
    # generate or load G grids
    if (!file.exists(path.grid.G)) {
      G.grid = gen.G.grid(n.hypo, G.step=G.step)
      saveRDS(G.grid, path.grid.G)
    } else {
      G.grid = readRDS(path.grid.G)
    }
    # generate data grids
    data.grid = merge(w.grid, G.grid)
    saveRDS(data.grid, path.grid.data)
  } else {
    # load data grids
    data.grid = readRDS(path.grid.data)
  }
  return(data.grid)
}

load.grid <- function(data.grid) {
  w.grid <- as.matrix(data.grid[, grepl('w', colnames(data.grid))])
  G.grid <- as.matrix(data.grid[, grepl('G', colnames(data.grid))])
  n_hypo = dim(w.grid)[2]
  n_G <- dim(G.grid)[1]
  dim(G.grid)  = c(n_G, n_hypo, n_hypo)
  G.grid <- aperm(G.grid, c(3, 2, 1))
  return(list('ws'=w.grid, 'Gs'=G.grid))
}
