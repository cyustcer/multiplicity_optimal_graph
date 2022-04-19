
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
