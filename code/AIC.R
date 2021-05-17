
times <- dat$sample_times
ph <- dat$pheno_H
pd <- dat$pheno_D
allpheno<-cbind(ph,pd)
mpheno <- as.numeric(colMeans(allpheno))
mmh <- max(mpheno[1:24])
mmd <- max(mpheno[-c(1:24)])
n<-dim(allpheno)[1]


double.curve <- function(par,times){
  
  a <- par[1]
  b <- par[2]
  k <- par[3]
  k*exp(-exp(a-b*times)) 
}


sum.curve <-function(par,yy,times) { 
  sum( (yy - double.curve(par,times) )^2 )
}

rss1<-0
for (i in 1:n) {
  ny <- ph[i,]
  ny<-c(t(ny))
  a <- max(ny)
  parin <- c(1,0.1,a)
  r <- optim(parin,sum.curve,yy=ny,times=times,method="BFGS",control=list(maxit=1500))
  rss1<-rss1+r$value
}

parin <- c(1,0.1,mmh)
m0 <- optim(parin,sum.curve,yy=mpheno[1:24],times=times,method="BFGS",control=list(maxit=1500))

tss1<-sum( (ph - mpheno[1:24] )^2 )
R2 <- 1-rss1/tss1
L1<-(1+log(2*pi)+log(rss1/n))*(-n)/2
AIC<-(6-2*L1)/n
#AIC <- 6+n*log(rss1/n)
adj.R2 <- 1-((n-1)/(n-3-1))*(1-R2)
RSD <- sqrt(rss1/n)
#________________________________________________________
rss2<-0
for (i in 1:n) {
  ny <- pd[i,]
  ny<-c(t(ny))
  a <- max(ny)
  parin <- c(1,0.1,a)
  r <- optim(parin,sum.curve,yy=ny,times=times,method="BFGS",control=list(maxit=1500))
  rss2<-rss2+r$value
}

parin <- c(1,0.1,mmd)
m0 <- optim(parin,sum.curve,yy=mpheno[-c(1:24)],times=times,method="BFGS",control=list(maxit=1500))

tss2<-sum( (pd - mpheno[-c(1:24)] )^2 )
R2 <- 1-rss2/tss2
L2<-(1+log(2*pi)+log(rss2/n))*(-n)/2
AIC<-(6-2*L2)/n
adj.R2 <- 1-((n-1)/(n-3-1))*(1-R2)
RSD <- sqrt(rss2/n)

R2 <- 1-(rss1+rss2)/(tss1+tss2)
adj.R2 <- 1-((n-1)/(n-6-1))*(1-R2)
RSD <- sqrt((rss1+rss2)/n)
L<-L1+L2
AIC<-(12-2*L)/n
SC<--2*L/n+(6*log(n))/n
HQ<--2*L/n+2*log(log(n))*6/n
####################################################################
double.curve <- function(par,times){
  
  a1 <- par[1]
  b1 <- par[2]
  k1 <- par[3]
  a2 <- par[4]
  b2 <- par[5]
  k2 <- par[6]
  
  k1*exp(-exp(a1-b1*times)) + k2*exp(-exp(a2-b2*times))
}


sum.curve <-function(par,yy,times) { 
  sum( (yy - double.curve(par,times) )^2 )
}

rss1<-0
for (i in 1:n) {
  ny <- ph[i,]
  #names(ny) <- NULL
  ny<-c(t(ny))
  a <- max(ny)/2
  parin <- c(1,0.1,a,1,0.1,a)
  r <- optim(parin,sum.curve,yy=ny,times=times,method="BFGS",control=list(maxit=1500))
  rss1<-rss1+r$value
}
parin <- c(1,0.1,mmh/2,1,0.1,mmh/2)
m0 <- optim(parin,sum.curve,yy=mpheno[1:24],times=times,method="BFGS",control=list(maxit=1500))

tss1<-sum( (ph - mpheno[1:24] )^2 )
R2 <- 1-rss1/tss1
L1<-(1+log(2*pi)+log(rss1/n))*(-n)/2
AIC<-(12-2*L1)/n
#AIC <- 6+n*log(rss1/n)
adj.R2 <- 1-((n-1)/(n-6-1))*(1-R2)
RSD <- sqrt(rss1/n)
#___________________________________________________________________________
rss2<-0
for (i in 1:n) {
  ny <- pd[i,]
  ny<-c(t(ny))
  a <- max(ny)/2
  parin <- c(1,0.1,a,1,0.1,a)
  r <- optim(parin,sum.curve,yy=ny,times=times,method="BFGS",control=list(maxit=1500))
  rss2<-rss2+r$value
}

parin <- c(1,0.1,mmd/2,1,0.1,mmd/2)
m0 <- optim(parin,sum.curve,yy=mpheno[-c(1:24)],times=times,method="BFGS",control=list(maxit=1500))

tss2<-sum( (pd - mpheno[-c(1:24)] )^2 )
R2 <- 1-rss2/tss2
L2<-(1+log(2*pi)+log(rss2/n))*(-n)/2
AIC<-(12-2*L2)/n
adj.R2 <- 1-((n-1)/(n-6-1))*(1-R2)
RSD <- sqrt(rss2/n)

R2 <- 1-(rss1+rss2)/(tss1+tss2)
adj.R2 <- 1-((n-1)/(n-12-1))*(1-R2)
RSD <- sqrt((rss1+rss2)/n)
L<-L1+L2
AIC<-(24-2*L)/n
SC<--2*L/n+(12*log(n))/n
HQ<--2*L/n+2*log(log(n))*12/n
#########################################################################
trible.G.curve <- function(par,times){
  
  a1 <- par[1]
  b1 <- par[2]
  k1 <- par[3]
  a2 <- par[4]
  b2 <- par[5]
  k2 <- par[6]
  a3 <- par[7]
  b3 <- par[8]
  k3 <- par[9]
  
  k1*exp(-exp(a1-b1*times)) + k2*exp(-exp(a2-b2*times)) + k3*exp(-exp(a3-b3*times))
}

sum.curve <-function(par,yy,times) { 
  sum( (yy - trible.G.curve(par,times) )^2 )
}

rss1<-0
for (i in 1:n) {
  ny <- ph[i,]
  ny<-c(t(ny))
  a <- max(ny)/3
  parin <- c(1,0.1,a,1,0.1,a,1,0.1,a)
  r <- optim(parin,sum.curve,yy=ny,times=times,method="BFGS",control=list(maxit=1500))
  rss1<-rss1+r$value
}

tss1<-sum( (ph - mpheno[1:24] )^2 )
R2 <- 1-rss1/tss1
L1<-(1+log(2*pi)+log(rss1/n))*(-n)/2
AIC<-(18-2*L1)/n
SC<--2*L1/n+(9*log(n))/n
HQ<--2*L1/n+2*log(log(n))*9/n
adj.R2 <- 1-((n-1)/(n-9-1))*(1-R2)
RSD <- sqrt(rss1/n)

#__________________________
rss2<-0
for (i in 1:n) {
  ny <- pd[i,]
  ny<-c(t(ny))
  a <- max(ny)/3
  parin <- c(1,0.1,a,1,0.1,a,1,0.1,a)
  r <- optim(parin,sum.curve,yy=ny,times=times,method="BFGS",control=list(maxit=1500))
  rss2<-rss2+r$value
}

tss2<-sum( (pd - mpheno[-c(1:24)] )^2 )
R2 <- 1-rss2/tss2
L2<-(1+log(2*pi)+log(rss2/n))*(-n)/2
AIC<-(18-2*L2)/n
SC<--2*L2/n+(9*log(n))/n
HQ<--2*L2/n+2*log(log(n))*9/n
adj.R2 <- 1-((n-1)/(n-9-1))*(1-R2)
RSD <- sqrt(rss2/n)

R2 <- 1-(rss1+rss2)/(tss1+tss2)
adj.R2 <- 1-((n-1)/(n-18-1))*(1-R2)
RSD <- sqrt((rss1+rss2)/n)
L<-L1+L2
AIC<-(36-2*L)/n
SC<--2*L/n+(18*log(n))/n
HQ<--2*L/n+2*log(log(n))*18/n
############################
source("ode2.gompertz.R")
rss<-0
for (i in 1:n) {
  ny <- allpheno[i,]
  ny<-c(t(ny))
  a1 <- max(ny[1:24])/2
  a2 <- max(ny[-c(1:24)])/2
  parin <- c(1,0.1,a1,1,0.1,a2,0.5,a1,0,0.5,a2,0)
  r <- optim(parin,s.mle,s.y=ny,s.t=times,x1=ny[1],x2=ny[25],method="BFGS",control=list(maxit=1500))
  rss<-rss+r$value
}
parin <- c(1,0.1,mmh/2,1,0.1,mmd/2,0.5,mmh/2,0,0.5,mmd/2,0)
m0 <- optim(parin,s.mle,s.y=mpheno,s.t=times,x1=mpheno[1],x2=mpheno[25],method="BFGS",control=list(maxit=1500))

tss<-sum( (allpheno - mpheno )^2 )
R2 <- 1-rss/tss
L<-(1+log(2*pi)+log(rss/n))*(-n)/2
AIC<-(24-2*L)/n
#AIC <- 6+n*log(rss1/n)
adj.R2 <- 1-((n-1)/(n-12-1))*(1-R2)
RSD <- sqrt(rss/n)
SC<--2*L/n+(12*log(n))/n
HQ<--2*L/n+2*log(log(n))*12/n
