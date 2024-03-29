#' line 6-282: basic functions
#' line 284-330: Type I error rates
#' line 332-391: Powers
#' line 393-455: Comparison on computation time

library(MASS)
library(PearsonDS)

#' Ep calculates the expectation for tr(AB)
Ep <- function(A, B){
    j = sum(diag(A))
    T = sum(diag(B))
    n = nrow(A)
    Expectation = j*T/(n-1)
    return(Expectation)
}

#' Vp calculates the variance for tr(AB)
Vp <- function(A, B){
    n = nrow(A)
    j = sum(diag(A))
    T = sum(diag(B))
    t2 = sum(diag(A %*% A))
    T2 = sum(diag(B %*% B))
    s2 = sum(diag(A)*diag(A))
    S2 = sum(diag(B)*diag(B))
    Variance = 
    (s2*S2*(n-1)^2*(n-2)*(n-3)+((j^2-s2)*(T^2-S2)+2*(t2-s2)*(T2-S2)+4*s2*S2)*(n-1)*(n-2)*(n-3)+(4*(2*s2-t2)*(2*S2-T2)+
    2*(2*s2-j^2)*(2*S2-T^2))*(n-1)*(n-3)+(2*t2-6*s2+j^2)*(2*T2-6*S2+T^2)*(n-1)- (j^2)*(T^2)*(n-2)*(n-3)*n)/(n*((n-1)^2)*(n-2)*(n-3))
    #s2*S2/n + ((j^2-s2)*(T^2-S2)+2*(t2-s2)*(T2-S2)+4*s2*S2)/(n*(n-1)) + (4*(2*s2-t2)*(2*S2-T2) +
    #2*(2*s2-j^2)*(2*S2-T^2))/(n*(n-1)*(n-2)) + (2*t2-6*s2+j^2)*(2*T2-6*S2+T^2)/(n*(n-1)*(n-2)*(n-3)) - (j^2)*(T^2)/((n-1)^2)
    return(Variance)
}

#' Tp calculates the third moment for tr(AB)
Tp <- function(A, B){
    n = nrow(A)
    j = sum(diag(A))
    T = sum(diag(B))
    t2 = sum(diag(A %*% A))
    T2 = sum(diag(B %*% B))
    s2 = sum(diag(A)*diag(A))
    S2 = sum(diag(B)*diag(B))
    t3 = sum(diag(A %*% A %*% A))
    T3 = sum(diag(B %*% B %*% B))
    s3 = sum(diag(A)*diag(A)*diag(A))
    S3 = sum(diag(B)*diag(B)*diag(B))
    u = sum(A*A*A)
    U = sum(B*B*B)
    d = sum(t(diag(A)) %*% diag(A %*% A))
    D = sum(t(diag(B)) %*% diag(B %*% B))
    b = sum(t(diag(A)) %*% A %*% diag(A))
    B = sum(t(diag(B)) %*% B %*% diag(B))
    ThirdMoment = ((n^2)*(n+1)*(n^2+15*n-4)*s3*S3+4*(n^4-8*n^3+19*n^2-4*n-16)*u*U+24*(n^2-n-4)*(u*B+b*U) 
    +6*(n^4-8*n^3+21*n^2-6*n-24)*b*B+12*(n^4-n^3-8*n^2+36*n-48)*d*D
    +12*(n^3-2*n^2+9*n-12)*(j*s2*D+d*T*S2)+3*(n^4-4*n^3-2*n^2+9*n-12)*j*T*s2*S2
    +24*((n^3-3*n^2-2*n+8)*(d*U+u*D)+(n^3-2*n^2-3*n+12)*(d*B+b*D))
    +12*(n^2-n+4)*(j*s2*U+u*T*S2)+6*(2*n^3-7*n^2-3*n+12)*(j*s2*B+b*T*S2)
    -2*n*(n-1)*(n^2-n+4)*((2*u+3*b)*S3+(2*U+3*B)*s3)
    -3*n*((n-1)^2)*(n+4)*((j*s2+4*d)*S3+(T*S2+4*D)*s3)
    +2*n*(n-1)*(n-2)*((j^3+6*j*t2+8*t3)*S3+(T^3+6*T*T2+8*T3)*s3)
    +(j^3)*((n^3-9*(n^2)+23*n-14)*T^3+6*(n-4)*T*T2+8*T3)
    +6*j*t2*((n-4)*T^3+(n^3-9*n^2+24*n-14)*T*T2+4*(n-3)*T3)
    +8*t3*(T^3+3*(n-3)*T*T2+(n^3-9*n^2+26*n-22)*T3)
    -16*((j^3)*U+u*(T^3))-12*(n^2-5*n+8)*(j*t2*U+u*T*T2)-8*(3*n^2-15*n+16)*(t3*U+u*T3)
    -6*(n^2-5*n+4)*((j^3)*B+b*(T^3))-24*((n^2-5*n+8)*(t3*B+b*T3)+(n^2-5*n+6)*(j*t2*B+b*T*T2))
    -3*(n-2)*(8*((j^3)*D+d*(T^3))+4*(n^2-5*n+12)*(j*t2*D+d*T*T2)
    +8*(n^2-5*n+8)*(t3*D+d*T3)+(n^2-5*n+2)*((j^3)*T*S2+(T^3)*j*s2)
    +2*(n^2-5*n+6)*(j*t2*T*S2+j*s2*T*T2)+16*(t3*T*S2+j*s2*T3)))/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))
    return(ThirdMoment)
}

#' Trans transform random data into genotype data
Trans <- function(x){
  n = nrow(x)
  m = ncol(x)
  for (i in 1:n) {
    for (j in 1:m) {
      q = 0.05+0.45*(j-1)/(m-1)
      if(x[i,j] < qnorm(q^2,0,1)) x[i,j] = 0
      else if (x[i,j] > qnorm(1-(1-q)^2,0,1)) x[i,j] = 2
      else x[i,j] = 1
    }
  }
  return(x)
}

#' f calculates the similarity matrix for genotype data x using IBS kernel
f <- function(x){
  n = nrow(x)
  m = ncol(x)  
  G = matrix(NA, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      G[i,j]=(sum(2-abs(x[i,]-x[j,])))/(2*m)
    }
  }
  return(G)
}

#' g calculates the similarity matrix for phenotype data x using Gaussian kernel
g <- function(x){
  n = nrow(x)
  G = matrix(NA, nrow = n, ncol = n)
  E = matrix(NA,1, n*(n-1)/2)
  for(i in 1:n){
    for(j in 1:n){
      G[i,j]=sum((x[i,]-x[j,])^2)
    }
  }
  k = 1 
  for(i in 2:n){
    for(j in 1:i-1){
      E[k]=sum((x[i,]-x[j,])^2)
      k = k+1  
    }
  }
  h = median(E)
  G1 = exp(-G/h)
  return(G1)
}

#' PVD calculates the p value when similarity matrix G and E are given
PVD <- function(G,E){
  n = ncol(G)
  H = diag(n) - (1/n)*rep(1,n)%*%t(rep(1,n))

  A = H %*% G %*% H
  B = H %*% E %*% H

  eigA = eigen(A)
  valA0 = eigA$values
  kA = sum(Re(valA0)>0.01)
  valA = valA0[1:kA]
  vecA = eigA$vectors[,1:kA]

  eigB = eigen(B)  
  valB0 = eigB$values
  kB = sum(Re(valB0)>0.01)
  valB = valB0[1:kB]
  vecB = eigB$vectors[,1:kB]

  X = c(1/16,1/4,1,4)
  P = rep(NA,4)
  
  for(i in 1:4){
    valA1 = valA^X[i]
    A1 = vecA %*% diag(valA1) %*% t(vecA)
    valB1 = valB^X[i]
    B1 = vecB %*% diag(valB1) %*% t(vecB)

    s = sum(A1*B1)
    
    x = Ep(A1,B1)
    y = Vp(A1,B1)
    z = Tp(A1,B1)

    a = (z-3*x*y-x^3)/(y^(3/2))
    b = 4/(a^2)

    tt = (s-x)/(sqrt(y))

    if(Re(a)>0){
      P[i] = ppearsonIII(q=Re(tt), shape=Re(b), location=Re(-sqrt(b)), scale=Re(1/sqrt(b)), lower.tail=F)
    }
    else{
      P[i] = ppearsonIII(q=Re(tt), shape=Re(b), location=Re(sqrt(b)), scale=Re(-1/sqrt(b)), lower.tail=F)
    }   
       
  }

  t0 = sum((1/4)*tan((1/2-P)*pi))
  p0 = 0.5 - atan(t0)/pi
 
  return(p0)
   
}

#' permCA calculates the p value when similarity matrix G and E and the times of permutation perN are given
permCA <- function(G,E,perN){
  n = ncol(G)
  H = diag(n) - (1/n)*rep(1,n) %*% t(rep(1,n))
  
  A = H %*% G %*% H
  B = H %*% E %*% H

  eigA = eigen(A)
  valA0 = eigA$values
  kA = sum(valA0>0.01)
  valA = valA0[1:kA]
  vecA = eigA$vectors[,1:kA]

  eigB = eigen(B)  
  valB0 = eigB$values
  kB = sum(valB0>0.01)
  valB = valB0[1:kB]
  vecB = eigB$vectors[,1:kB]

  X = c(1/16,1/4,1,4)
  P = rep(NA,4)

  for(i in 1:4){
    valA1 = valA^X[i]
    A1 = vecA %*% diag(valA1) %*% t(vecA)
    valB1 = valB^X[i]
    B1 = vecB %*% diag(valB1) %*% t(vecB)

    s = sum(A1*B1)
    P1 = rep(NA, perN)

    for(j in 1:perN){
        ind1 = sample(n,n)
        vecB.per = vecB[ind1,]
        B2 = vecB.per %*% diag(valB1) %*% t(vecB.per)
        P1[j] = sum(B2*A1)
    }

    P[i] = sum(s<P1)/perN
  }

  t0 = sum((1/4)*tan((1/2-P)*pi))
  p0 = 0.5 - atan(t0)/pi
  return(p0)
}

#' approCAP calculates the p value when similarity matrix G and E and the times of permutation perN are given
approCAP <- function(G,E,perN){
  n = ncol(G)
  H = diag(n) - (1/n)*rep(1,n) %*% t(rep(1,n))
  
  A = H %*% G %*% H
  B = H %*% E %*% H

  eigA = eigen(A)
  valA0 = eigA$values
  kA = sum(valA0>0.01)
  valA = valA0[1:kA]
  vecA = eigA$vectors[,1:kA]

  eigB = eigen(B)  
  valB0 = eigB$values
  kB = sum(valB0>0.01)
  valB = valB0[1:kB]
  vecB = eigB$vectors[,1:kB]

  X = c(1/16,1/4,1,4)
  P = rep(NA,4)

  for(i in 1:4){
    valA1 = valA^X[i]
    A1 = vecA %*% diag(valA1) %*% t(vecA)
    valB1 = valB^X[i]
    B1 = vecB %*% diag(valB1) %*% t(vecB)

    s = sum(A1*B1)
    P1 = rep(NA, perN)

    for(j in 1:perN){
        ind1 = sample(n,n)
        vecB.per = vecB[ind1,]
        B2 = vecB.per %*% diag(valB1) %*% t(vecB.per)
        P1[j] = sum(B2*A1)
    }
    
    a = calc.skew(P1)
    b = 4/(a^2)

    tt = (s-mean(P1))/(sqrt(var(P1)))

    if(Re(a)>0){
      P[i] = ppearsonIII(q=Re(tt), shape=Re(b), location=Re(-sqrt(b)), scale=Re(1/sqrt(b)), lower.tail=F)
    }
    else{
      P[i] = ppearsonIII(q=Re(tt), shape=Re(b), location=Re(sqrt(b)), scale=Re(-1/sqrt(b)), lower.tail=F)
    }   
       
  }

  t0 = sum((1/4)*tan((1/2-P)*pi))
  p0 = 0.5 - atan(t0)/pi
  return(p0)
}

#' Error calculates type 1 error for the test 
#' ( p is the parameter in the covariance matrix; m1 is dim(geno); m2 is dim(pheno) )
Error <- function(p, m1, m2){
  Dx = matrix(NA, m1, m1)
  for(i in 1:m1){
    for(j in 1:m1){
    Dx[i,j]=0.5^(abs(i-j))
    }
  }
  Dy = matrix(NA, m2, m2)
  for(i in 1:m2){
    for(j in 1:m2){
    Dy[i,j]=p^(abs(i-j))
    }
  }
  #' or Dy = (1-p)*diag(m2)+p*rep(1, m2) %*% t(rep(1, m2))
  Z1 = rep(NA, 1000)

  for(i in 1:1000){
    X = mvrnorm(200, rep(0,m1), Dx)
    X = Trans(X)
    Y = mvrnorm(200, rep(0,m2), Dy)
    Z1[i] = PVD(f2(X), g1(Y))
  }
  return(Z1) 
}

#' Theo calculates theoratical points of p value
#' (N times)
Theo <- function(N){
  Z2 = rep(NA, N)
  for(i in 1:N){
    Z2[i] = (i-0.5)/N
  }
  Z2 = -log(Z2)/log(10)

  return(Z2)
}

m1 = 30 #' or 60
m2 = 30 #' or 60
p = 0.2 #' or 0.5 or 0.8
X<- Error(p, m1, m2)
X1 = sort(X)
X2 = -log(X1)/log(10)
Y<- Theo(1000)
plot(X2,Y,ylab=expression("log"[0.1]*"(Theoretical value)"),xlab=expression("log"[0.1]*"(P-value)"),type = "p",main =expression("m"[1]*"=30, "*"m"[2]*"=30, "*symbol(r)*"=0.2"));abline(a=0,b=1)

#' Power calculates powers for the test 
#' ( p is the parameter in the covariance matrix; m1 is dim(geno); m2 is dim(pheno) )
Power <- function(p){
  Po = rep(NA,7)
  for(I in 1:7){
    M = (I+1)*50
    m1 = 60
    m2 = 60
    Dx = matrix(NA, m1, m1)
    for(i in 1:m1){
      for(j in 1:m1){
      Dx[i,j]=0.5^(abs(i-j))
      }
    }
    Dy = matrix(NA, m2, m2)
    for(i in 1:m2){
      for(j in 1:m2){
      Dy[i,j]=p^(abs(i-j))
      }
    }
#' or Dy = (1-p)*diag(m2)+p*rep(1, m2) %*% t(rep(1, m2))

    Z1 = rep(NA, 1000)

    betax = matrix(NA, m1, m2)
      for(j in 1:m2){
        if(j < 11 ){
          betax[,j] = 0.005 + 0.005*(j-1)/9
        }
        else{
          betax[,j] = 0
        }
      }

    for(i in 1:1000){
      X = mvrnorm(M, rep(0,m1), Dx)
      X = Trans(X)
      Y = exp(X) %*% betax + mvrnorm(M, rep(0, m2), Dy)
    
      Z1[i] = PVD(f2(X), g1(Y))
    }

    Po[I] = sum(Z1<0.05)/1000
  }

  return(Po)
}

#' Theo calculates sample sizes
Theor <- function(N){
  Z2 = rep(NA, N)
  for(i in 1:N){
    Z2[i] = (i+1)*50
  }
  return(Z2)
}

Y<- Power(0.2) #' or 0.5 or 0.8
X<- Theor(7)
plot(X,Y,ylab="Power",xlab="Sample size",type = "b",main =expression("pho"*"=0.2"))

m1=30 #' or 60
m2=30 #' or 60
Dx = matrix(NA, m1, m1)
for(i in 1:m1){
  for(j in 1:m1){
  Dx[i,j]=0.5^(abs(i-j))
  }
}
Dy = matrix(NA, m2, m2)
for(i in 1:m2){
  for(j in 1:m2){
  Dy[i,j]=0.5^(abs(i-j))
  }
}

X = mvrnorm(200, rep(0,m1), Dx)
X = Trans(X)
Y = mvrnorm(200, rep(0,m2), Dy)

G=f2(X)
E=g1(Y)

time_start<-Sys.time()# # Record the time when the code starts executing
PVD(G,E)
PVD(G,E)
PVD(G,E)
PVD(G,E)
PVD(G,E)
PVD(G,E)
PVD(G,E)
PVD(G,E)
PVD(G,E)
PVD(G,E)
exc_time<-difftime(Sys.time(),time_start,units = 'mins')# Record the time difference between the current time and time_start
print(paste0('Execution time of the code',round(exc_time,2),'mins'))# Output: Keep the time difference of two decimal points

time_start<-Sys.time()# # Record the time when the code starts executing
permCA(G,E,100) #' or permCA(G,E,1000) or permCA(G,E,10000)
permCA(G,E,100) #' or permCA(G,E,1000) or permCA(G,E,10000)
permCA(G,E,100) #' or permCA(G,E,1000) or permCA(G,E,10000)
permCA(G,E,100) #' or permCA(G,E,1000) or permCA(G,E,10000)
permCA(G,E,100) #' or permCA(G,E,1000) or permCA(G,E,10000)
permCA(G,E,100) #' or permCA(G,E,1000) or permCA(G,E,10000)
permCA(G,E,100) #' or permCA(G,E,1000) or permCA(G,E,10000)
permCA(G,E,100) #' or permCA(G,E,1000) or permCA(G,E,10000)
permCA(G,E,100) #' or permCA(G,E,1000) or permCA(G,E,10000)
permCA(G,E,100) #' or permCA(G,E,1000) or permCA(G,E,10000)
exc_time<-difftime(Sys.time(),time_start,units = 'mins')# Record the time difference between the current time and time_start
print(paste0('Execution time of the code',round(exc_time,2),'mins'))# Output: Keep the time difference of two decimal points

time_start<-Sys.time()# # Record the time when the code starts executing
approCAP(G,E,100) #' or approCAP(G,E,1000) or approCAP(G,E,10000)
approCAP(G,E,100) #' or approCAP(G,E,1000) or approCAP(G,E,10000)
approCAP(G,E,100) #' or approCAP(G,E,1000) or approCAP(G,E,10000)
approCAP(G,E,100) #' or approCAP(G,E,1000) or approCAP(G,E,10000)
approCAP(G,E,100) #' or approCAP(G,E,1000) or approCAP(G,E,10000)
approCAP(G,E,100) #' or approCAP(G,E,1000) or approCAP(G,E,10000)
approCAP(G,E,100) #' or approCAP(G,E,1000) or approCAP(G,E,10000)
approCAP(G,E,100) #' or approCAP(G,E,1000) or approCAP(G,E,10000)
approCAP(G,E,100) #' or approCAP(G,E,1000) or approCAP(G,E,10000)
approCAP(G,E,100) #' or approCAP(G,E,1000) or approCAP(G,E,10000)
exc_time<-difftime(Sys.time(),time_start,units = 'mins')# Record the time difference between the current time and time_start
print(paste0('Execution time of the code',round(exc_time,2),'mins'))# Output: Keep the time difference of two decimal points