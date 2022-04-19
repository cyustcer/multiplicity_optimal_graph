library(gMCP)

graph.power <- function(w, G, type.i.error, pvals, weights) {
  graph = matrix2graph(G, w)
  out.success = graphTest(pvalues=pvals, graph=graph, alpha=type.i.error)
  return(out.success)
}


weighted.power <- function(out.success, weight) {
  out.power = apply(out.success, 2, mean) %*% weight / sum(weight)
  return(out.power)
}

disjunctive.power <- function(out.success) {
  out.power = mean(apply(out.success, 1, sum) > 0)
  return(out.power)
}

conjunctive.power <- function(out.success) {
  n.hypo = dim(out.success)[2]
  out.power = mean(apply(out.success, 1, sum) == n.hypo)
  return(out.power)
}
