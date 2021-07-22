#' CheckBounds2
#'
#' @description Finds smallest s so that x+s*dx is within bounds
#'
#' @param x  k times 1 vector, the given value of x0
#' @param dx  k times 1 vector, the Newton step from x0 to x1
#' @param bounds k times 2 matrix, the first (second) colum store the lower (upper) bounds of x
#'
#' @return s a scalar
#' @export


CheckBounds2 <- function(x, dx, bounds){
  # browser()
  l1 <- (bounds[,1]-x)/dx
  s1 <- min(l1[l1 > 0])
  l2 <- (bounds[,2]-x)/dx
  s2 <- min(l2[l2 > 0])
  s1 <- min(c(0.98*s1,0.98*s2,1))
  return(s1)
}





