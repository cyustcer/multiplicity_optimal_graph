# Recursively calculates weights for a given graph and intersection hypothesis
#
# Recursively calculates weights for a given graph and intersection hypothesis
#
# @param h A numeric vector with only binary entries (0,1).
# @param g Graph represented as transition matrix.
# @param w A numeric vector of weights.
# @return A weight vector.
# @author Florian Klinglmueller \email{float@@lefant.net}
# @keywords graphs
# @examples
#
# g <- BonferroniHolm(4)@m
# w <- rep(1/4, 4)
# h <- c(0,1,0,1)
# gMCP:::mtp.weights(h,g,w)
# gMCP:::mtp.edges(h,g,w)
#
mtp.weights <- function(h, g, w){
  ## recursively compute weights for a given graph and intersection hypothesis
  if(sum(h) == length(h)){
    return(w)
  } else {
    j <- which(h == 0)[1]
    h[j] <- 1
    wu <- mtp.weights(h, g, w)
    gu <- mtp.edges(h, g, w)
    guj <- gu[j, ]
    wt <- wu + wu[j] * guj
    wt[j] <- 0
    return(wt)
  }
}
