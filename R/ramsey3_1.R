#' Nonlinear equations (3.2) from section 3.1 in HM(2009)
#'
#' @description we don't write a equation set, and get \code{k} through \code{c},
#' and it reduce the number of equation solved.
#'
#' @param k a vector which length denotes the number of periods from \code{k0} to \code{kT}
#' @param k0 initial value for capital stock
#' @param kT final value for capital stock
#' @param alpha a parameter for production function
#' @param delta rate of depreciation
#' @param betax discount rate for unitity function
#' @param eta a parameter for unitity funciton with flexible labor supply
#'
#' @importFrom magrittr `%>%`
#'
#' @export

ramsey <- function(k,k0,kT,alpha = 0.27, delta = 0.011,betax = 0.994,eta = 2){
# browser()
  fx <- matlab::zeros(length(k),1) %>% as.numeric()
  x <- c(k0,k,kT)

  c0 <- (x[1]^alpha)+(1-delta)*x[1] - x[2]
  if (c0 < 0){
    fx[1] <- NA
    # return(fx)
  }

  for (tt in 1:length(k)){
    c1 <- (x[tt+1]^alpha)+(1-delta)*x[tt+1] - x[tt+2];
    if (c1 < 0){
      # print(tt)
      fx[tt] <- NA
      # return(fx)
    }
    fx[tt] <- (c0^(-eta)) - betax*(c1^(-eta))*(1-delta+alpha*(x[tt+1]^(alpha-1)))
    c0 <- c1
  }
  return(fx)
}
