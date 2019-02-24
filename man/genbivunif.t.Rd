\name{genbivunif.t}
\alias{genbivunif.t}
\title{Simulating Bivariate Uniform Data via Beta(1-t,1+t)}
\description{This function simulates bivariate uniform data via Beta(1-t,1+t) and a specified correlation (rho) as described in Demirtas (2014).}
\usage{
genbivunif.t(N=10000, rho, print.cor=TRUE)
}
\arguments{
  \item{N}{The sample size for the bivariate uniform data to be simulated. Default is 10,000.}
  \item{rho}{Theoretical correlation for the uniform data to be simulated.}
  \item{print.cor}{Option to print correlation results. Default is TRUE.}
}
\value{A list of length 3 containing the simulated bivariate uniform data, the theoretical correlation specified by the user, and the empirical correlation of the simulated data titled unif.dat, specified.rho, and empirical.rho, respectively.}
\examples{
set.seed(567)
res.gent<-genbivunif.t(N=10000, rho=0.9)
#"Specified rho is 0.9 and empirical rho is 0.897883."
}
