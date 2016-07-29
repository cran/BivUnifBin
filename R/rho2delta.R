rho2delta<-function(N=10000, p1, p2, rho, print.cor=TRUE) {
  if(p1<=0 | p1>=1) {
    stop("p for variable 1 must be between 0 and 1.")
  }
  if(p2<=0 | p2>=1) {
    stop("p for variable 2 must be between 0 and 1.")
  }
  
  #solve for a
  a<-(-5/2)+(1/2)*sqrt((rho+49)/(rho+1))
  
  #generate uniform data
  unif.dat<-genbivunif.a(N=N, rho=rho, print.cor=FALSE)$unif.dat
  
  #generate binary data and calculate empirical delta
  q1<-1-p1
  c1<-quantile(unif.dat$x, prob=q1)
  w1<-ifelse(unif.dat$x>c1, 1, 0)
  
  q2<-1-p2
  c2<-quantile(unif.dat$y, prob=q2)
  w2<-ifelse(unif.dat$y>c2, 1, 0)
  
  bin.dat<-data.frame(w1=w1, w2=w2)
  emp.delta<-round(cor(w1,w2), 6)
  
  #calculate algorithmic delta
  if(p1<p2 & p1+p2<=1) { #region1
    Fp1p2<-((1)/(2*(a+1)))*((p1+p2)^(a+1)-(p2-p1)^(a+1))
  } else if(p1>=p2 & p1+p2<=1) { #region2
    Fp1p2<-((1)/(2*(a+1)))*((p1+p2)^(a+1)-(p1-p2)^(a+1))
  } else if(p1<=p2 & p1+p2>=1) { #region3
    Fp1p2<-p2-1+p1+((1)/(2*(a+1)))*((2-p1-p2)^(a+1)-(p2-p1)^(a+1))
  } else if(p1>p2 & p1+p2>1) { #region4
    Fp1p2<-p2-1+p1+((1)/(2*(a+1)))*((2-p1-p2)^(a+1)-(p1-p2)^(a+1))
  }
  
  alg.delta<-round((Fp1p2-p1*p2)/sqrt(p1*q1*p2*q2), 6)
  
  if(print.cor==TRUE) {
    print(paste("Algorithmic delta is ", alg.delta, " and empirical delta is ", emp.delta, ".", sep=""))
  }
  return(list(unif.dat=unif.dat, 
              bin.dat=bin.dat,
              specified.rho=rho, 
              algorithmic.delta=alg.delta,
              empirical.delta=emp.delta)) 
  
}