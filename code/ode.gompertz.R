require("deSolve")
phasic1 <- function(par,times){
  
  a <- par[1]
  b <- par[2]
  k <- par[3]
  k*exp(-exp(a-b*times))
}


phasic2 <- function(par, times,x1,x2, options=list())
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      r1 = par[1],
      k1 = par[2],
      a1 = par[3],
      r2 = par[4],
      k2 = par[5],
      a2 = par[6]);
  }
  
  state0 <- c(X=x1, Y=x2);
  y <- COMP.f( par0, state0, times );
  
  #if (class(y)=="try-error" )
  #  return(rep(NaN, length(times)*2));
  index.time <- 1:length(times)
  return ( c(y[index.time, 2:3] ) );
}

COMP.f <-function( parameters, state, times ){
  Lorenz<-function(t, state, parameters) 
  {
    with( as.list(c(state, parameters)),
          {
            dX <- r1*X*(1-X/k1) + a1*Y*X
            dY <- r2*Y*(1-Y/k2) + a2*X*Y
            
            list(c(dX, dY))
          }
          
    ) # end with(as.list ...
  }
  
  out <- try(rk(y = state, times = times, func = Lorenz, parms = parameters,method="rk4") );
  out;
}

yy<-function(par,times,x1,x2){
  H<-phasic1(par[1:3],times)+phasic2(par[-c(1:6)], times,x1,x2)[1:length(times)]
  D<-phasic1(par[4:6],times)+phasic2(par[-c(1:6)], times,x1,x2)[(1+length(times)):(2*length(times))]
  return(c(H,D))
}

s.mle <- function(s.par,s.y,s.t,x1,x2){
  fit.y<-yy(s.par,s.t,x1,x2)
  A <- sum((s.y - fit.y)^2 )
  A
}


com.get_mu1 <- function(par, times,x1,x2, options=list())
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      r1 = par[1],
      k1 = par[2],
      r2 = par[3],
      k2 = par[4]);
  }
  
  state0 <- c(X=x1, Y=x2);
  y <- COMP.f1( par0, state0, times );
  
  #if (class(y)=="try-error" )
  #  return(rep(NaN, length(times)*2));
  index.time <- 1:length(times)
  return ( c(y[index.time, 2:3] ) );
}


COMP.f1 <-function( parameters, state, times ){
  Lorenz<-function(t, state, parameters) 
  {
    with( as.list(c(state, parameters)),
          {
            dX <- r1*X*(1-(X/k1))
            dY <- r2*Y*(1-(Y/k2))
            
            list(c(dX, dY))
          }
          
    ) # end with(as.list ...
  }
  
  out <- try(rk(y = state, times = times, func = Lorenz, parms = parameters,method="rk4") );
  out;
}


com.get_mu2 <- function(par, times,x1,x2, options=list())
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      a1 = par[1],
      a2 = par[2]);
  }
  
  state0 <- c(X=x1, Y=x2);
  y <- COMP.f2( par0, state0, times );
  
  #if (class(y)=="try-error" )
  #  return(rep(NaN, length(times)*2));
  index.time <- 1:length(times)
  return ( c(y[index.time, 2:3] ) );
}

COMP.f2 <-function( parameters, state, times ){
  Lorenz<-function(t, state, parameters) 
  {
    with( as.list(c(state, parameters)),
          {
            dX <- a1*Y*X
            dY <- a2*X*Y
            
            list(c(dX, dY))
          }
          
    ) # end with(as.list ...
  }
  
  out <- try(rk(y = state, times = times, func = Lorenz, parms = parameters,method="rk4") );
  out;
}

R2<-function(par,H1,D1,nt){
  
  r1 = par[1]
  k1 = par[2]
  a1 = par[3]
  r2 = par[4]
  k2 = par[5]
  a2 = par[6]
  Z<-phasic2(par,nt,x1=H1,x2=D1)
  X<-Z[1:length(nt)]
  Y<-Z[(length(nt)+1):(2*length(nt))]
  RH<-r1*X*(1-X/k1) + a1*Y*X
  RD<-r2*Y*(1-Y/k2) + a2*X*Y
  return(list(RH=RH,RD=RD))
}

R1 <- function(par,times){
  a <- par[1]
  b <- par[2]
  k <- par[3]
  
  k*b*exp(-exp(a-b*times))*exp(a-b*times)
}

R2.ind<-function(par,H1,D1,nt){
  
  r1 = par[1]
  k1 = par[2]
  r2 = par[3]
  k2 = par[4]
  Z<-com.get_mu1(par,nt,x1=H1,x2=D1)
  X<-Z[1:length(nt)]
  Y<-Z[(length(nt)+1):(2*length(nt))]
  RH<-r1*X*(1-X/k1)
  RD<-r2*Y*(1-Y/k2)
  return(list(RH=RH,RD=RD))
}
