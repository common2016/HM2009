#' MNRStep
#'
#' @description   Find the step size s so that the Newton-Raphson algorithm
#'      always moves in the direction of a (local) minimum of (1/2)(f(x)'f(x))
#'
#' @param x0 n times 1 vector, the initial point
#' @param dx0 n times 1 vector, the Newton direction
#' @param dg  1 times n vector, the gradient of (1/2)(f'f) at x_0
#' @param pTol scalar, parameter convergence criterium
#' @param f pointer to the function whose zero x solves the system of equations
#' @param ... any additional arguments passed to \code{f}.
#'
#' @return s  admissible stepsize
#'
#' @export


MNRStep <- function(x0,dx0,dg,pTol,f,...){
  # Fixed parameters of the algorithm
  smult <- 1.0e-4
  smin <- 0.1
  smax <- 0.5

  # Initialize
  s1 <- 1.0
  amat <- matlab::zeros(2,2)
  bvec <- matlab::zeros(2,1)
  g0 <- (1/2)*sum(f(x0,...)*f(x0,...))
  dgdx <- dg*dx0                        #dg(x0)*dx0
  g1 <- (1/2)*sum(f(x0+dx0,...) * f(x0+dx0,...))

  # Try the full Newton step s=1
  if (g1 <= g0 + smult*dgdx){
    return(s1)
  } else {
    s <- -dgdx/(2*(g1-g0-dgdx))
    if (s < smin) s <- smin
    if (s > smax) s <- smax
    x1 <- x0+s*dx0
    g2=(1/2)*sum(f(x1, ...) * f(x1, ...))
    if (is.na(g2)) return(g2)
  }
  s2 <- s

 # Reduce s2 further unless g2 < g0 + s2*smult*dgdx */
 while (g2 > (g0 + smult*s2*dgdx) ){
   amat[1,1] <- 1/(s2^2)
   amat[1,2] <- -1/(s1^2)
   amat[2,1] <- -s1/(s2^2)
   amat[2,2] <- s2/(s1^2)
   bvec[1] <- g2-s2*dgdx-g0
   bvec[2] <- g1-s1*dgdx-g0
   ab <- (amat*bvec)/(s2-s1)

   if (ab[1,1] == 0.0){
     s <- -dgdx/(2*ab[2,1])
   } else {
     disc <- (ab[2,1]^2)-3*ab[1,1]*dgdx
     if (disc < 0.0){
       s=s2*smax;
     } else if (ab[2,1] <= 0.0){
       s=(-ab[2,1]+sqrt(disc))/(3*ab[1,1]);
     } else {
       s <- -dgdx/(ab[2,1]+sqrt(disc))
     }
   }


   if (s < s2*smin) s <- s2*smin
   if (s > s2*smax) s <- s2*smax

   #tol=sqrt((s*dx0)'(s*dx0))/(1+sqrt(x0'x0))
  if (s < MinStep(x0,s*dx0,matlab::ones(length(x0),1),pTol)) return(s)
   # if tol < pTol; retp(s); endif
  s1 <- s2
  s2 <- s
  g1 <- g2
  x1 <- x0+s2*dx0
  g2 <- (1/2)*sum(f(x1, ...) * f(x1, ...))
  if (is.na(g2)) return(g2)
 }
  return(s2)
}



