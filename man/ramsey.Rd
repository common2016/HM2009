% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ramsey3_1.R
\name{ramsey}
\alias{ramsey}
\title{Nonlinear equations (3.2) from section 3.1 in HM(2009)}
\usage{
ramsey(k, k0, kT, alpha = 0.27, delta = 0.011, betax = 0.994, eta = 2)
}
\arguments{
\item{k}{a vector which length denotes the number of periods from \code{k0} to \code{kT}}

\item{k0}{initial value for capital stock}

\item{kT}{final value for capital stock}

\item{alpha}{a parameter for production function}

\item{delta}{rate of depreciation}

\item{betax}{discount rate for unitity function}

\item{eta}{a parameter for unitity funciton with flexible labor supply}
}
\description{
we don't write a equation set, and get \code{k} through \code{c},
and it reduce the number of equation solved.
}
