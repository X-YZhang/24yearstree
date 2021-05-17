com.DH.load <- function(datt,pheno,pheno1,snpinfo,nmiss,nstart,nlen,times){
  
  
  mergedat <- fg_read_data(datt,pheno,nstart, nlen, times)
  mergedat1 <- fg_read_data(datt,pheno1,nstart, nlen, times)
  
  info <- read.csv(snpinfo)
  
  final <- get_geno_code(mergedat$snp.info, mergedat$gen, mergedat$phe, mergedat$ids)
  
  final1 <- get_geno_code(mergedat1$snp.info, mergedat1$gen, mergedat1$phe, mergedat1$ids)
  
  nobs <- mergedat1$n.obs
  gene_table <- matrix(final$gen0,nrow=nobs,byrow=T) 
  colnames(gene_table) <- mergedat$snp.info[,2]
  
  snp <- t(gene_table)
  new_D <- final$phe0
  new_H <- final1$phe0
  dattt <- list(pheno_D=new_D[,1:24],pheno_H=new_H[,1:24],snps=snp,sample_times=1:24,snpinfo=info)
  return(dattt)
  
}


fg_read_phe<-function( phe.csv )
{
  tb.phe<- read.csv(phe.csv);
  y.dim <- dim(tb.phe);
  
  return(list(dat=tb.phe[,-1], id=tb.phe[,1], n.obs=dim(tb.phe)[1], times=c(1:(y.dim[2]-1)), phe.csv=phe.csv ))
}

#scaffoldId, Loci, RefBase, AltBase, P1.A, P1.B, P2.A, P2.B....
fg_read_data<-function(geno.dat, phe.csv, nstart, nlen, times=NA)
{
  tb.gen<-read.table(geno.dat, header=T);
  tb.phe<-read.csv(phe.csv);
  
  cat("gen.table[", dim(tb.gen), "]\n");
  
  tb.ids<-colnames(tb.gen);
  tb.ids<-tb.ids[-c(1:4)]
  tb2.ids<-tb.phe$X
  tb2.idx <- c();
  for(i in 1:length(tb.ids))
  {
    tb2.idx<-c(tb2.idx, which(tb2.ids==tb.ids[i]) );
  }
  
  if (length(tb2.idx)>0)
    tb.phe<- tb.phe[tb2.idx,]
  
  nstop<-(nstart+nlen-1);
  if (nstop > dim(tb.gen)[1])
    nstop <- dim(tb.gen)[1];
  
  if (nstart==0 && nlen==0)
  {
    nstop <- dim(tb.gen)[1];
    nstart<- 1;
  }		
  
  n.obs <- dim(tb.phe)[1];
  times  <- 1:(dim(tb.phe)[2]-1); 
  
  snp.info <- tb.gen[c(nstart:nstop),1:4];
  gen   <- tb.gen[c(nstart:nstop), c(5:(n.obs+4))]
  n.snp <- dim(gen)[1];
  
  return(list(ids=tb.phe[,1], n.obs=n.obs, n.snp=n.snp, times=times, 
              snp.info=snp.info, gen=gen, phe=tb.phe[,-1], geno.dat=geno.dat, phe.csv=phe.csv ))
}

get_geno_code<-function(d.snp, d.gen, d.phe, d.ids)
{
  d.gen2 <- as.character(unlist(d.gen));
  
  snpB <- as.character(unlist(d.snp[4]));
  snpA <- as.character(unlist(d.snp[3]));
  
  QQ2<- paste(snpB, snpB, sep=""); 
  qq0<- paste(snpA, snpA, sep="");
  Qq1<- list(paste(snpA, snpB, sep=""), paste( snpB, snpA, sep="") ) ;
  
  d.g <- array( 9, length(d.ids) );
  d.g[which(d.gen2==QQ2)]<-2;
  d.g[which(d.gen2==qq0)]<-0;
  d.g[which(d.gen2==c(Qq1[[1]]))]<-1;
  d.g[which(d.gen2==c(Qq1[[2]]))]<-1;
  
  #if( fg.sys$log >= LOG_DEBUG ) show(d.gen);
  
  d.umiss <- d.g;
  if (length(which(d.g==9))>0)
  {
    d.ids <- d.ids[-which(d.g==9)];
    d.phe <- d.phe[-which(d.g==9),];
    d.g <- d.g[-which(d.g==9)];
  }
  
  return(list(gen0=d.g, phe0=d.phe, ids0=d.ids, umiss=d.umiss))
}




double.parin <- function(dat,ret=res,test.file="testcross.csv",
                         inter.file="intercross.csv",n.t=1,n.i=1){
  
  test <- read.csv(test.file)[,-1]
  inter <- read.csv(inter.file)[,-1]
  
  testloci <- test[order(test[,5]),1][n.t]
  interloci <- inter[order(inter[,5]),1][n.i]
  
  
  testloci.chr <- test[order(test[,5]),3][n.t]
  interloci.chr <- inter[order(inter[,5]),3][n.i]
  
  snpnames <- dat$snpinfo[,3]
  
  testindex <- snpnames[testloci]
  
  interindex <- snpnames[interloci]
  
  parin <- ret[c(testloci,interloci),1:43]
  parin <- cbind(parin,c(testloci.chr,interloci.chr),c(testloci,interloci))
  rownames(parin) <- c(testindex,interindex)
  #parin[order(parin[,24]),]
  parin
}


Impute<-function(Z, impute.method){
  
  p<-dim(Z)[2]
  
  if(impute.method =="random"){
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-rbinom(length(IDX),2,maf1)
      }
    }
  } else if(impute.method =="fixed"){
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-2 * maf1
      }
    }
  } else if(impute.method =="bestguess") {
    
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-round(2 * maf1)
      }
    }
    
  } else {
    stop("Error: Imputation method shoud be \"fixed\", \"random\" or \"bestguess\" ")
  }
  
  return(Z)
}


com.parin <- function(dat,ret=ret,test.file="testcross.csv",
                      inter.file="intercross.csv",n.t=c(1,2,3,4,6),n.i=1:5){
  
  test <- read.csv(test.file)[,-1]
  inter <- read.csv(inter.file)[,-1]
  
  testloci <- test[order(test[,5]),1][n.t]
  interloci <- inter[order(inter[,5]),1][n.i]
  
  
  testloci.chr <- test[order(test[,5]),3][n.t]
  interloci.chr <- inter[order(inter[,5]),3][n.i]
  
  snpnames <- dat$snpinfo[,3]
  testindex <- snpnames[testloci]
  interindex <- snpnames[interloci]
  
  test.lineno<-dat$snpinfo[testloci,1]
  inter.lineno<-dat$snpinfo[interloci,1]
  
  parin <- ret[c(testloci,interloci),1:25]
  parin <- cbind(parin,c(testloci.chr,interloci.chr),c(testloci,interloci))
  #rownames(parin) <- c(testindex,interindex)
  #parin[order(parin[,24]),]
  parin
}