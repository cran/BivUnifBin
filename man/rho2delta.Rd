\name{rho2delta}
\alias{rho2delta}
\title{Bivariate Uniform Data and its Connection to Bivariate Binary Data}
\description{This function simulates bivariate uniform data with the given correlation (rho) using \code{\link{genbivunif.a}}, dichotomizes the data based on p1 and p2 to generate bivariate binary data, and calculates the algorithmic and empirical correlation (delta) of the bivariate binary data as described in Demirtas (2014).}
\usage{
rho2delta(N=10000, p1, p2, rho, print.cor=TRUE)
}
\arguments{
  \item{N}{The sample size for the data to be simulated. Default is 10,000.}
  \item{p1}{The expected value of the first variable w1 to be dichotomized.}
  \item{p2}{The expected value of the second variable w2 to be dichotomized.}
  \item{rho}{Theoretical correlation for the bivariate uniform data to be simulated.}
  \item{print.cor}{Option to print correlation results. Default is TRUE.}
}
\value{A list of length 5 containing simulated bivariate uniform data, binary data derived from dichotomizing the uniform data, theoretical correlation specified by the user, algorithmic correlation of the binary data, and the empirical correlation of the binary data, titled unif.dat, bin.dat, specified.rho, algorithmic.delta, and empirical.delta, respectively.}
\seealso{\code{\link{genbivunif.a}}}
\examples{
set.seed(234)
res.r2d<-rho2delta(N=10000, p1=0.4, p2=0.7, rho=0.7)
#"Algorithmic delta is 0.393887 and empirical delta is 0.395547."
}
