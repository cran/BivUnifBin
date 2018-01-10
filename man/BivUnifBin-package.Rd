\name{BivUnifBin-package}
\alias{BivUnifBin-package}
\alias{BivUnifBin}
\docType{package}
\title{
Generation of Bivariate Uniform Data and its Relation to Bivariate Binary Data
}
\description{
This package implements the algorithms described in Demirtas (2014) for generation of bivariate uniform data with coverage on the whole feasible correlation range and its relationship to bivariate binary data. Two functions for simulation of bivariate uniform data are provided: \code{genbivunif.a} simulates bivariate uniform data via g(u)~Beta(a,1) and a specified correlation (the tetrachoric correlation, which is referred to as rho), and \code{genbivunif.t} simulates bivariate uniform data via g(u)~Beta(1-t,1+t) and a specified correlation rho.

The computations based on the relationship between bivariate uniform data and binary data are implemented in \code{rho2delta} and \code{delta2rho}. In \code{rho2delta}, bivariate uniform data are simulated given a sample size and correlation rho. The data are then dichotomized based on the given expected values for binary data, and an empirical and algorithmic correlation of the binary data (the phi coefficient, which is referred to as delta) are provided. In \code{delta2rho}, algorithmic rho is calculated given a specified delta. Bivariate uniform data are simulated based on the calculated algorithmic rho. The uniform data are dichotomized based on the given expected values for binary data, and an empirical and algorithmic correlation of the binary data are provided. Before any computations, the specified delta is compared to the lower and upper correlation bounds for the two binary variables using the generate, sort and correlate (GSC) algorithm in Demirtas and Hedeker (2011).
}
\details{
\tabular{ll}{
Package: \tab BivUnifBin\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2018-01-10\cr
License: \tab GPL-2 | GPL-3\cr
}}
\author{
Hakan Demirtas, Rawan Allozi

Maintainer: Rawan Allozi <ralloz2@uic.edu>
}
\references{
Demirtas, H. and Hedeker, D. (2011). A practical way for computing approximate lower and upper correlation bounds. \emph{The American Statistician}, \bold{65(2)}, 104-109.

Demirtas, H. (2014). Generating bivariate uniform data with a full range of correlations and connections to bivariate binary data. \emph{Communications in Statistics-Theory and Methods}, \bold{43(17)}, 3574-3579.
}
