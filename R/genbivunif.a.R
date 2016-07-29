genbivunif.a<-function(N=10000, rho, print.cor=TRUE) {
  tol<-.Machine$double.eps^0.5
  if(N<=0 | abs(N - round(N)) > tol) {
    stop("N must be a positive integer.")
  }
  if(rho<=-1 | rho>=1) {
    stop("rho can take values between -1 and 1.")
  }
  a<-(-5/2)+(1/2)*sqrt((rho+49)/(rho+1))
  
  x<-runif(N)
  v2<-runif(N)
  u<-rbeta(n=N, shape1=a, shape2=1)
  
  y<-ifelse(v2<0.5, abs(u-x), 1-abs(1-u-x))
 
  e.rho<-round(cor(x,y), 6)
  if(print.cor==TRUE) {
    print(paste("Specified rho is ", round(rho, 6), " and empirical rho is ", e.rho, ".", sep=""))
  }
  return(list(unif.dat=data.frame(x,y), specified.rho=round(rho,6), empirical.rho=e.rho)) 
}

