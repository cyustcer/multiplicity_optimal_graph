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
