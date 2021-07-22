#' Parameters Test
#'
#' @description Compute relative change in x
#'
#' @param x n times 1 vector
#' @param dx n times 1 vector, change in x
#' @param typx n times 1 vector, typical elements of x
#'
#' @return crit scalar
#' @export

ParTest <- function(x,dx,typx){
  n <- length(x)
  crit <- matlab::zeros(n,1)
  i <- 1
  while (i <= n){
    crit[i] <- abs(dx[i])/max(c(abs(x[i]),abs(typx[i])))
    i <- i+1
  }

  return(max(crit))

}





