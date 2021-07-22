#' Min Step
#'
#' @description Compute the minimal step size so that x1=x0+mstep*dx=0
#' (in terms of the parameter tolerance criterium pTol)
#'
#' @param x n times 1 vector
#' @param dx n times 1 vector, change in x
#' @param typx n times 1 vector, typical elements of x
#' @param pTol scalar
#'
#' @return mstep scalar

MinStep <- function(x,dx,typx,pTol){
  n <- length(x)
  i <- 1
  converge <- 0

  while (i <= n){
    temp <- abs(dx[i])/max(c(abs(x[i]),abs(typx[i])))
    if (temp > converge){
      converge <- temp
    }
    i=i+1;
  }
  return(pTol/converge)
}





