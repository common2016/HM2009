#' FixvMN2
#' @description  Solves a system of non-linear equations using a modified Newton Method
#' @param f Pointer the vector valued function F(x), whose
#'                  zero, F(x1)=0, is to be computed.
#' @param bounds k timex 2 vector of upper and lower bounds on x
#' @param x0  k times 1 vector of starting values
#' @param stopc stopping criterium
#' @param pTol parameter convergence criterium
#' @param MNR_CDJac wheather using \code{CDJac} computing Jacobian matrix, with fault value is TRUE.
#' @param MNR_Global wheather using line search, with fault value is FALSE.
#' @param ... any additional arguments passed to \code{f}.
#'
#' @importFrom matlab ones
#' @return a list including x1 and crit. x1 is k times 1 vector, the approximate solution to F(x1)=0.
#'
#'  \code{crit} is 1 vector, where
#'  \itemize{
#'       \item \itemize{
#'     \item \code{crit[1]=0} : normal termination
#'    \item      =1 : function evaluation failed,
#'    \item      =2 : no further decrease in function value possible,
#'     \item =3 : maximum number of iterations exceedes
#'              }
#'       \item \code{crit[2]}   : termination criterion: max(abs(f(x)))
#'       \item \code{crit[3]}   : the maximum relative change of x between the last two iterations
#'       \item \code{crit[4]}   : f(x)'f(x)/2
#'       \item \code{crit[5]}   : number of iterations
#'       }
#' @export
#'
#' @examples
#' # para setting
#' library(matlab)
#' alpha <- 0.27
#' betax <- 0.994
#' tmax <- 60
#'
#' x0 <- ones(tmax,1)
#' bounds <- data.frame(dw = ones(tmax,1)*0.01, up = ones(tmax,1)*30)
#' FixVMN2(x0 = x0,bounds = bounds, f = ramsey, k0 = 4.4, kT = 0, delta = 0.011)
#'
#' # another function to solve it
#' rootSolve::multiroot(f = ramsey,start = x0,k0 = 4.4, kT = 0, delta = 0.011)





FixVMN2 <- function(x0,bounds,f, stopc = 1e-7, pTol = .Machine$double.eps^(2/3),
                    MNR_CDJac = TRUE, MNR_Global = FALSE, ...){

  # Initialize
    maxit <- 500      # stop after 500 Iterations
    x1 <- x0
    typf <- f(x0, ...)

    # Start Iterations
    crit <- matlab::ones(5,1)
    crit[1] <- 0
    critold <- 2

   #  if _MNR_Print; DosWin; cls; endif; # clear screen if output is printed to the screen
    while (crit[5] <= maxit) {# start iterations

      # method of computing Jacbian
      if (MNR_CDJac){
        df <- CDJac(f,x1,...)
      } else {
        df <- numDeriv::jacobian(f,x1,...)
      }
# browser()
      # initial value evaluating
      if (any(is.na(df))){
        crit[1] <- 1
        return(list(x1 = x1,crit = crit))
      }
      fx <- f(x1,...)
      if (MNR_Global)  dg <- matrix(fx,nrow = 1) %*% df
      dx <- optR::optR(df,-fx, method = 'LU')$beta[[1]] # use the LU factorization
      x2 <- x1 + dx
      step1 <- 1
      mstep <- MinStep(x1,dx,matlab::ones(length(x1),1),pTol)
      step1 <- CheckBounds2(x1,dx,bounds) # scale down by step if array bounds are exceeded
    #   if MNR_Print;
    # locate 4,2;
    # ?"mStep= " mStep;
    # locate 5,2;
    # ?"Step1= " step1;
    # endif;
    if ((mstep <= 1) & (step1 < mstep)){
      crit[1] <- 2
      return(list(x1,crit = crit))
    }
    dx <- step1*dx
    if (MNR_Global){
      step2 <- MNRStep(x1,dx,dg,pTol,f,...)
    } else step2 <- 1

    # if _MNR_Print;
    # locate 6,2;
    # ?"Step2= " step2;
    # if ismiss(step2); crit[1]=1; retp(x1,crit); endif;
    # endif;
    x2 <- x1 + step2*dx
    crit[2] <- max(abs(f(x2,...)))
    crit[3] <- ParTest(x1,step2*dx,matlab::ones(length(x1),1))
    if (crit[2] < stopc){
      crit[1] <- 0
      return(list(x1 = x2,crit  = crit))
    } else {
      if (step1*step2 < mstep){
        crit[1] <- 2
        return(list(x1 = x2,crit = crit))
      }
    }


    if (MNR_Global){
      critold <- crit[4]
      crit[4] <- sum(f(x2,...) * f(x2,...))/2
    }
    x1 <- x2
    crit[5] <- crit[5]+1
    }
    if (crit[5] >= maxit) crit[1] <- 3
    return(list(x1 = x1,crit = crit))
}



