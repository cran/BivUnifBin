\name{delta2rho}
\alias{delta2rho}
\title{Bivariate Binary Data and its Connection to Bivariate Uniform Data}
\description{This function calculates algorithmic rho given delta, simulates bivariate uniform data with the calculated rho using \code{\link{genbivunif.a}}, dichotomizes the data based on p1 and p2 to generate bivariate binary data, and calculates the algorithmic and empirical correlation (delta) of the bivariate binary data as described in Demirtas (2014).}
\usage{
delta2rho(N=10000, p1, p2, delta, print.cor=TRUE)
}
\arguments{
  \item{N}{The sample size for the data to be simulated. Default is 10,000.}
  \item{p1}{The expected value of the first binary variable w1.}
  \item{p2}{The expected value of the second binary variable w2.}
  \item{delta}{Theoretical correlation for the bivariate binary data to be simulated.}
  \item{print.cor}{Option to print correlation results. Default is TRUE.}
}
\value{A list of length 5 containing simulated bivariate uniform data, binary data derived from dichotomizing the uniform data, correlation of the binary data specified by the user, empirical correlation of the binary data, and algorithmic correlation of the uniform data, titled unif.dat, bin.dat, specified.delta, empirical.delta, and algorithmic.rho, respectively.}
\seealso{\code{\link{genbivunif.a}}}
\examples{
set.seed(147)
res.d2r<-delta2rho(N=10000, p1=0.4, p2=0.7, delta=0.5)
#"Specified delta is 0.5, empirical delta is 0.488197, and algorithmic rho is 0.909093."
}
