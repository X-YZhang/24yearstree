com.DH.H0 <- function(dat){
  
  allpheno <- cbind(dat$pheno_H,dat$pheno_D)
  times <- dat$sample_times
  
  parin2 <- c( 0.700027831,0.919989750,0.775633812,1.548642740,0.310478718,
               1.914807661,0.177739581,14.708957769,2.976438785,0.162247023,14.953668322,
               0.560697795,15.580939795,-0.007372783,0.521424779,18.601865700,0.009362770)
  mpheno <- as.numeric(colMeans(allpheno))
  loop_k <- 1
  max_iter <- 100
  epsi <- 10^-4
  max_err <- 1
  
  while(loop_k<max_iter && max_err>epsi){
    
    
    oldpar <-c(parin2);
    mle.covar1 <- function(npar){
      
      nnpar <- c(npar,parin2[6:17])
      AA <- curve.mle(nnpar,y=allpheno,time.std =times,x1=mpheno[1],x2=mpheno[25])
      AA
    }
    r1.covar <- optim(parin2[1:5],mle.covar1,method = "BFGS",control=list(maxit=32000))
    new.covar1 <- r1.covar$par
    #cat("new.coavr1:",unlist( new.covar1), "\n");
    
    mle.1 <- function(npar){
      
      nnpar <- c(new.covar1,npar)
      AA <- curve.mle(nnpar,y=allpheno,time.std =times,x1=mpheno[1],x2=mpheno[25])
      AA
    }
    r1 <- optim(c(parin2[6:17]),mle.1,method = "BFGS",control=list(maxit=32000))    
    new1 <- r1$par
    
    #r2<- optim(c(new.covar1,new1),curve.mle,y=allpheno,time.std =times,x1=mpheno[1],x2=mpheno[15],
    #method = "BFGS",control=list(maxit=32000,trace=T))    
    #cat("new1:",unlist( new1), "\n");
    
    nparin <- c(new.covar1,new1)
    
    newpar <- c(nparin)
    #cat("newpar:", newpar, "\n");
    
    max_err <- max(abs( oldpar - newpar) );
    
    parin2 <- nparin
    #cat(loop_k, "max.err=", max_err, "allpar", newpar,"\n");
    loop_k <- loop_k+1; 
  }
  LL <- curve.mle(parin2,y=allpheno,time.std =times,x1=mpheno[1],x2=mpheno[25])
  return(c(LL,parin2))
}



curve.mle <-function( par,y,time.std,x1,x2)
{
  len.cov <- 5
  par.covar <- par[1:len.cov]
  n  <- length(y[,1])
  #sig.inv3 <- SAD3.get_inv_mat(par.covar,time.std, 2)
  sigma <- SAD3.get_mat(par.covar,time.std, 2)#solve( sig.inv3 )
  
  curve.par <- par[(len.cov+1):(len.cov+ 12)]
  mu <- yy(curve.par,time.std,x1,x2)
  
  yy <- y
  fy <- dmvnorm(yy,mu,sigma)
  #fy[which(fy<=.Machine$double.eps)] <- .Machine$double.eps
  A <- -sum(log(fy))
  #cat("LL=",A,"\n")
  return(A)
}


double.H0.debug <- function(dat){
  
  ny <- cbind(dat$pheno_H,dat$pheno_D)
  times <- dat$sample_times
  mpheno <- as.numeric(colMeans(ny))
  mmh <- max(mpheno[1:24])/2
  mmd <- max(mpheno[-c(1:24)])/2
  init.par <- c(1.9,0.1,mmh,2.5,0.1,mmd,0.5,mmh,0,0.5,mmd,0)
  
  c0 <- optim(init.par,s.mle,s.y=mpheno,s.t=times,x1=mpheno[1],x2=mpheno[25],
              method="BFGS",control=list(maxit=1500))
  covar.par <- c(0.5,sd(dat$pheno_H[,12]),0.5,sd(dat$pheno_D[,12]),0.5)
  parin <- c(covar.par,c0$par)
  time.std <- times
  
  loop_k <- 1
  max_iter <- 100
  epsi <- 10^-4
  max_err <- 1
  
  while(loop_k<max_iter && max_err>epsi){
    
    
    oldpar <-c(parin);
    mle.covar1 <- function(npar){
      
      nnpar <- c(npar,parin[6:17])
      AA <- curve.mle(nnpar,y=ny,time.std =times,x1=mpheno[1],x2=mpheno[25])
      AA
    }
    r1.covar <- optim(parin[1:5],mle.covar1,method = "BFGS",control=list(maxit=32000))
    new.covar1 <- r1.covar$par
    #cat("new.coavr1:",unlist( new.covar1), "\n");
    
    mle.1 <- function(npar){
      
      nnpar <- c(new.covar1,npar)
      AA <- curve.mle(nnpar,y=ny,time.std =times,x1=mpheno[1],x2=mpheno[25])
      AA
    }
    r1 <- optim(c(parin[6:17]),mle.1,method = "BFGS",control=list(maxit=32000))    
    new1 <- r1$par
    #cat("new1:",unlist( new1), "\n");
    
    nparin <- c(new.covar1,new1)
    
    newpar <- c(nparin)
    #cat("newpar:", newpar, "\n");
    
    max_err <- max(abs( oldpar - newpar) );
    
    parin <- nparin
    cat(loop_k, "max.err=", max_err, "allpar", newpar,"\n");
    loop_k <- loop_k+1; 
  }
  LL <- curve.mle(parin,y=ny,time.std =times,x1=mpheno[1],x2=mpheno[25])
  return(c(LL,parin))
}
