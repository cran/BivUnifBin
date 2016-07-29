delta2rho<-function(N=10000, p1, p2, delta, print.cor=TRUE) {
  if(p1<=0 | p1>=1) {
    stop("p for variable 1 must be between 0 and 1.")
  }
  if(p2<=0 | p2>=1) {
    stop("p for variable 2 must be between 0 and 1.")
  }
  
  #ensure that phi coefficient is within feasible range
  corr.limits<-valid.limits.BinOrdNN(plist=list(p1, p2), no.bin=2, no.ord=0, no.NN=0)
  
  if(delta<corr.limits$lower[2,1] | delta>corr.limits$upper[2,1]) {
    stop(paste('Specified delta is not within the feasible correlation range of [',
               corr.limits$lower[2,1],', ', corr.limits$upper[2,1], '] for the given distributional characteristics.', sep=''))
  }  
  
  #find theoretical rho
  q1<-1-p1
  q2<-1-p2
  b<-delta*(p1*q1*p2*q2)^(1/2)+(p1*p2)
  
  if(p1<p2 & p1+p2<=1) { #region1
    find.a<-function(a) {
      (((1)/(2*(a[1]+1)))*((p1+p2)^(a[1]+1)-(p2-p1)^(a[1]+1)))-b
    }
  } else if(p1>=p2 & p1+p2<=1) { #region2
    find.a<-function(a) {
      (((1)/(2*(a[1]+1)))*((p1+p2)^(a[1]+1)-(p1-p2)^(a[1]+1)))-b
    }
  } else if(p1<=p2 & p1+p2>=1) { #region3
    find.a<-function(a) {
      (p2-1+p1+((1)/(2*(a[1]+1)))*((2-p1-p2)^(a[1]+1)-(p2-p1)^(a[1]+1)))-b
    }
  } else if(p1>p2 & p1+p2>1) { #region4
    find.a<-function(a) {
      (p2-1+p1+((1)/(2*(a[1]+1)))*((2-p1-p2)^(a[1]+1)-(p1-p2)^(a[1]+1)))-b
    }
  }
  
  a<-multiroot(f=find.a, start=1)$root
  alg.rho<-round(((1-a)*(6-a))/((2+a)*(3+a)), 6)

  #generate uniform data
  unif.dat<-genbivunif.a(N=N, rho=alg.rho, print.cor=FALSE)$unif.dat
  
  #generate binary data and calculate empirical delta
  q1<-1-p1
  c1<-quantile(unif.dat$x, prob=q1)
  w1<-ifelse(unif.dat$x>c1, 1, 0)
  
  q2<-1-p2
  c2<-quantile(unif.dat$y, prob=q2)
  w2<-ifelse(unif.dat$y>c2, 1, 0)
  
  bin.dat<-data.frame(w1=w1, w2=w2)
  emp.delta<-round(cor(w1,w2), 6)
  
  if(print.cor==TRUE) {
    print(paste("Specified delta is ", delta, ", empirical delta is ", emp.delta,", and algorithmic rho is ", alg.rho,".", sep=''))
  }
  return(list(unif.dat=unif.dat,
              bin.dat=bin.dat,
              specified.delta=delta,
              empirical.delta=emp.delta,
              algorithmic.rho=alg.rho)) 
}