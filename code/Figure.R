Figure1<-function(dat){
  ph <- dat$pheno_H
  pd <- dat$pheno_D
  times<-dat$sample_times
  nt<-seq(1,24,0.1)
  allpheno<-cbind(ph,pd)
  mpheno <- as.numeric(colMeans(allpheno))
  mmh <- max(mpheno[1:24])
  mmd <- max(mpheno[-c(1:24)])
  n<-dim(allpheno)[1]
  
  sigle.G.curve <- function(par,times){
    
    a <- par[1]
    b <- par[2]
    k <- par[3]
    k*exp(-exp(a-b*times)) 
  }
  
  double.G.curve <- function(par,times){
    
    a1 <- par[1]
    b1 <- par[2]
    k1 <- par[3]
    a2 <- par[4]
    b2 <- par[5]
    k2 <- par[6]
    
    k1*exp(-exp(a1-b1*times)) + k2*exp(-exp(a2-b2*times))
  }
  
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
  
  ######################### G n=1
  sum.curve <-function(par,yy,times) { 
    sum( (yy - sigle.G.curve(par,times) )^2 )
  }
  
  parin <- c(1,0.1,mmh)
  m0 <- optim(parin,sum.curve,yy=mpheno[1:24],times=times,method="BFGS",control=list(maxit=1500))
  h1<-sigle.G.curve(m0$par,nt)
  
  parin <- c(1,0.1,mmd)
  m1 <- optim(parin,sum.curve,yy=mpheno[-c(1:24)],times=times,method="BFGS",control=list(maxit=1500))
  d1<-sigle.G.curve(m1$par,nt)
  
  ######################### G n=2
  sum.curve <-function(par,yy,times) { 
    sum( (yy - double.G.curve(par,times) )^2 )
  }
  
  parin <- c(1,0.1,mmh/2,1,0.1,mmh/2)
  m0 <- optim(parin,sum.curve,yy=mpheno[1:24],times=times,method="BFGS",control=list(maxit=1500))
  h2<-double.G.curve(m0$par,nt)
  
  parin <- c(1,0.1,mmd/2,1,0.1,mmd/2)
  m1 <- optim(parin,sum.curve,yy=mpheno[-c(1:24)],times=times,method="BFGS",control=list(maxit=1500))
  d2<-double.G.curve(m1$par,nt)
  
  ######################### G n=3
  sum.curve <-function(par,yy,times) { 
    sum( (yy - trible.G.curve(par,times) )^2 )
  }
  
  parin <- c(1,0.1,mmh/3,1,0.1,mmh/3,1,0.1,mmh/3)
  m0 <- optim(parin,sum.curve,yy=mpheno[1:24],times=times,method="BFGS",control=list(maxit=1500))
  h3<-trible.G.curve(m0$par,nt)
  
  parin <- c(1,0.1,mmd/3,1,0.1,mmd/3,1,0.1,mmd/3)
  m1 <- optim(parin,sum.curve,yy=mpheno[-c(1:24)],times=times,method="BFGS",control=list(maxit=1500))
  d3<-trible.G.curve(m1$par,nt)
  
  parin <- c(1,0.1,mmh/2,1,0.1,mmd/2,0.5,mmh/2,0,0.5,mmd/2,0)
  m3 <- optim(parin,s.mle,s.y=mpheno,s.t=times,x1=mpheno[1],x2=mpheno[25],method="BFGS",control=list(maxit=1500))
  y3<-yy(m3$par,nt,x1=mpheno[1],x2=mpheno[25])[1:length(nt)]
  y31<-phasic2(m3$par[-c(1:6)], nt,x1=mpheno[1],x2=mpheno[25])[1:length(nt)]
  ind.fit<-com.get_mu1(m3$par[-c(1:6)][c(1,2,4,5)],nt,x1=mpheno[1],x2=mpheno[25])[1:length(nt)]
  #inter.fit<-com.get_mu2(m3$par[-c(1:6)][c(3,6)],nt,x1=mpheno[1],x2=mpheno[25])[1:length(nt)]
  inter.fit<-y31-ind.fit
  y32 <- phasic1(m3$par[1:3],nt)
  
  R.lv<-R2(par=m3$par[-c(1:6)],H1=mpheno[1],D1=mpheno[25],nt)
  R.lv.ind<-R2.ind(par=m3$par[-c(1:6)][c(1,2,4,5)],H1=mpheno[1],D1=mpheno[25],nt)
  R.Logistic.D<-R1(m3$par[4:6],nt)
  R.Logistic.H<-R1(m3$par[1:3],nt)
  
  mt2<-nt[which(R.lv$RH==max(R.lv$RH))]
  my<-phasic2(m3$par[-c(1:6)], seq(1,mt2,0.1),x1=mpheno[1],x2=mpheno[25])[length(seq(1,mt2,0.1))]
  t22<-2*(max(y31)-my)/max(R.lv$RH)
  #t21<-2*(my)/max(R.lv$RH)
  a2<-seq(1,(mt2+t22),0.1)
  a2p<-mt2+t22 
  #b2p<-mt2-t22 
  
  mt2.ind<-nt[which(R.lv.ind$RH==max(R.lv.ind$RH))]
  my.ind<-com.get_mu1(m3$par[-c(1:6)][c(1,2,4,5)], seq(1,mt2.ind,0.1),x1=mpheno[1],x2=mpheno[25])[length(seq(1,mt2.ind,0.1))]
  a2p.ind<-mt2.ind+2*(max(ind.fit)-my.ind)/max(R.lv.ind$RH)
  a2.ind<-seq(1,a2p.ind,0.1)
  
  mt1<-m3$par[1]/m3$par[2]
  
  #col2rgb("mistyrose")
  ymax <- apply(dat$pheno_H,2,sort)[66,]
  dat1 <- data.frame(x=1:24,y=c(ymax))
  ymaxl <- loess(y ~ x, dat1)
  ymin <- apply(dat$pheno_H,2,sort)[1,]
  dat2 <- data.frame(x=1:24,y=c(ymin))
  yminl <- loess(y ~ x, dat2)
  tiff("./Figures/figure1.tiff",width=20,height=10,units = "cm",res=300)
  par(mfrow=c(1,2))
  par(mar=c(4,4,2,1))
  plot(c(0,0), c(2,0),xlim=c(0,25), ylim=c(-6,31), xlab="Age (year)",
       ylab = "Stem Height Growth (m)",type="n", xaxt="n", yaxt="n",
       xaxs="i", yaxs="i",cex.lab=1.2,mgp=c(2.5,1,0));
  polygon(c(seq(1,24,0.1),seq(24,1,-0.1)),c(predict(ymaxl,seq(1,24,0.1)),
                                            predict(yminl,seq(24,1,-0.1))),
          col = rgb(255,165,79,alpha=100,max=255),lwd = 0.01, border = NA)
  lines(nt,h1,col=rgb(0,205,0,alpha=220,max=255),lwd=2)
  lines(nt,h2,col="black",lwd=2)
  lines(nt,h3,col=rgb(0,0,255,alpha=220,max=255),lwd=2)
  lines(nt,y3,col=rgb(255,0,0,alpha=220,max=255),lwd=2)
  lines(nt,y32,col="orchid",lwd=2)
  # lines(a2,y31[1:length(a2)],col="royalblue",lwd=2)
  # lines(a2.ind,ind.fit[1:length(a2.ind)],col="royalblue",lwd=2,lty=2)
  # lines(a2.ind,inter.fit[1:length(a2.ind)],col="royalblue",lwd=2,lty=3)
  lines(nt,y31,col="orangered",lwd=2)
  lines(nt,ind.fit,col="orangered",lwd=2,lty=2)
  lines(nt,inter.fit,col="orangered",lwd=2,lty=3)
  axis(1,seq(1,24,3),seq(1,24,3),cex.axis=1)
  axis(2,seq(-5,30,5),seq(-5,30,5),las = 1,cex.axis=1)
  abline(h=0)
  mtext("A",side=3, cex=1.5,adj=-0.17,padj =-0.2)
  
  
  y3<-yy(m3$par,nt,x1=mpheno[1],x2=mpheno[25])[(length(nt)+1):(2*length(nt))]
  y31<-phasic2(m3$par[-c(1:6)], nt,x1=mpheno[1],x2=mpheno[25])[(length(nt)+1):(2*length(nt))]
  ind.fit<-com.get_mu1(m3$par[-c(1:6)][c(1,2,4,5)],nt,x1=mpheno[1],x2=mpheno[25])[(length(nt)+1):(2*length(nt))]
  #inter.fit<-com.get_mu2(m3$par[-c(1:6)][c(3,6)],nt,x1=mpheno[1],x2=mpheno[25])[(length(nt)+1):(2*length(nt))]
  inter.fit<-y31-ind.fit
  y32 <- phasic1(m3$par[4:6],nt)
  
  mt2<-nt[which(R.lv$RD==max(R.lv$RD))]
  my<-phasic2(m3$par[-c(1:6)], seq(1,mt2,0.1),x1=mpheno[1],x2=mpheno[25])[2*length(seq(1,mt2,0.1))]
  t22<-2*(max(y31)-my)/max(R.lv$RD)
  #t21<-2*(my)/max(R.lv$RD)
  a2<-seq(1,(mt2+t22),0.1)
  a2p<-mt2+t22 
  #b2p<-mt2-t22 
  
  mt2.ind<-nt[which(R.lv.ind$RD==max(R.lv.ind$RD))]
  my.ind<-com.get_mu1(m3$par[-c(1:6)][c(1,2,4,5)], seq(1,mt2.ind,0.1),x1=mpheno[1],x2=mpheno[25])[2*length(seq(1,mt2.ind,0.1))]
  a2p.ind<-mt2.ind+2*(max(ind.fit)-my.ind)/max(R.lv.ind$RD)
  a2.ind<-seq(1,a2p.ind,0.1)
  
  mt1<-m3$par[4]/m3$par[5]
  
  ymax <- apply(dat$pheno_D,2,sort)[66,]
  dat1 <- data.frame(x=1:24,y=c(ymax))
  ymaxl <- loess(y ~ x, dat1)
  ymin <- apply(dat$pheno_D,2,sort)[1,]
  dat2 <- data.frame(x=1:24,y=c(ymin))
  yminl <- loess(y ~ x, dat2)
  
  par(mar=c(4,4,2,1))
  plot(c(0,0), c(2,0),xlim=c(0,25), ylim=c(-6,46), xlab="Age (year)",
       ylab = "Stem Diameter Growth (cm)",type="n", xaxt="n", yaxt="n",
       xaxs="i", yaxs="i",cex.lab=1.2,mgp=c(2.5,1,0));
  polygon(c(seq(1,24,0.1),seq(24,1,-0.1)),c(predict(ymaxl,seq(1,24,0.1)),
                                            predict(yminl,seq(24,1,-0.1))),
          col = rgb(255,165,79,alpha=100,max=255),lwd = 0.01, border = NA)
  lines(nt,d1,col=rgb(0,205,0,alpha=220,max=255),lwd=2)
  lines(nt,d2,col="black",lwd=2)
  lines(nt,d3,col=rgb(0,0,255,alpha=220,max=255),lwd=2)
  lines(nt,y3,col=rgb(255,0,0,alpha=220,max=255),lwd=2)
  lines(nt,y32,col="orchid",lwd=2)
  # lines(a2,y31[1:length(a2)],col="royalblue",lwd=2)
  # lines(a2.ind,ind.fit[1:length(a2.ind)],col="royalblue",lwd=2,lty=2)
  # lines(a2.ind,inter.fit[1:length(a2.ind)],col="royalblue",lwd=2,lty=3)
  lines(nt,y31,col="orangered",lwd=2)
  lines(nt,ind.fit,col="orangered",lwd=2,lty=2)
  lines(nt,inter.fit,col="orangered",lwd=2,lty=3)
  axis(1,seq(1,24,3),seq(1,24,3),cex.axis=1)
  axis(2,seq(-5,46,10),seq(-5,46,10),las = 1,cex.axis=1)
  abline(h=0)
  mtext("B",side=3, cex=1.5,adj=-0.17,padj =-0.2)
  dev.off()
}

Figure2 <- function(ret=res,dat,file1="testcross.csv",file2="intercross.csv"){
  
  require(ggplot2)
  SNP <- dat$snps
  LR <- ret[,1]
  
  LR[which(is.na(LR))] <- 1
  LR[which(LR < 0)] <- 1
  #LR[nnot] <- 1
  snpnames <- 1:dim(SNP)[1]
  testcross <- c()
  intercross <- c()
  for(i in 1:length(LR)){
    index <- which(is.na(SNP[i,]))
    if(length(index)>0){
      nsnp <- SNP[i,-index]
    }else{
      nsnp <- SNP[i,]
    }
    tt <- table(nsnp)
    if(length(tt)==2){
      testcross<- c(testcross,i)
    }
    if(length(tt)==3){
      intercross<- c(intercross,i)
    }
  }
  df <- rep(NA,length(LR))
  df[testcross] <- 12
  df[intercross] <- 24
  p.value <- pchisq(LR,df,lower.tail = F)
  snpinfo <- dat$snpinfo
  map <- snpinfo[,1:5]
  scaffoldsum <- length(table(map[,2]))
  scaffoldc <- as.numeric(names(table(map[,2])))
  
  for(i in 1:scaffoldsum){
    map[which(map[,2]==scaffoldc[i]),2] <- i
  }
  
  
  for (i in 1:scaffoldsum){ 
    ndx <- which(map[, 2]==i)
    lstMrk <- max(map[ndx, 3])
    if (i < scaffoldsum) ndx2 <-which(map[, 2]==i+1)
    if (i < scaffoldsum) map[ndx2, 3] <- map[ndx2, 3] + lstMrk
  }
  
  bpMidVec <- vector(length=scaffoldsum)
  
  for (i in 1:scaffoldsum){
    ndx <- which(map[, 2]==i)
    posSub <- map[ndx, 3]
    bpMidVec[i] <- ((max(posSub) - min(posSub))/2) + min(posSub)
  }
  p.value<-p.adjust(p.value,method="fdr",n=length(p.value))
  tp.value <- p.value
  ip.value <- p.value
  tp.value[intercross] <- NA
  ip.value[testcross] <- NA
  logp1 <- -log(tp.value)
  nmap1 <- cbind(map,logp1)
  colnames(nmap1) <- c('lineNo', 'scaffold', 'loci',"Reference.Base","Alternative.base",'logp')
  logp2 <- -log(ip.value)
  nmap2 <- cbind(map,logp2)
  colnames(nmap2) <- c('lineNo', 'scaffold', 'loci',"Reference.Base","Alternative.base", 'logp')
  
  sig1 <- -log(1e-40/156362)
  tmp1 <- as.matrix(tp.value[which(logp1 > sig1)])
  sigsnp1 <- cbind(snpnames[which(logp1 > sig1)],map[which(logp1 > sig1),],tmp1)
  #write.csv(sigsnp1,file=file1)
  sig2 <- -log(1e-50/156362)
  tmp2 <- as.matrix(ip.value[which(logp2 > sig2)])
  sigsnp2 <- cbind(snpnames[which(logp2 > sig2)],map[which(logp2 > sig2),],tmp2)
  #write.csv(sigsnp2,file=file2)
  
  p1 <- ggplot(nmap1)
  p1 <- p1 + geom_point(aes(x=loci, y=logp,colour=as.factor(scaffold)),size=0.6)  
  p1 <- p1 + scale_color_manual(values=rep(c('#6495ED', '#76EE00'), 668))  
  p1 <- p1 + scale_x_continuous(breaks=c(bpMidVec[1:19]),labels=c(1:19))
  p1 <- p1 + scale_y_continuous(limits=c(0,150),breaks=c(seq(0,150,20)),labels=c(seq(0,150,20)))
  p1 <- p1 + geom_hline(yintercept=sig1,linetype=1, col='red', lwd=0.6)
  p1 <- p1 + annotate("text", x = (bpMidVec[1]), y = 145, label=c("Testcross"),size=3.5)
  p1 <- p1 + xlab('Chromosome') + ylab('-log (p-value)')+theme_zg()
  
  p2 <- ggplot(nmap2)
  p2 <- p2 + geom_point(aes(x=loci, y=logp,colour=as.factor(scaffold)),size=0.6)  
  p2 <- p2 + scale_color_manual(values=rep(c('#6495ED', '#76EE00'), 668))  
  p2 <- p2 + scale_x_continuous(breaks=c(bpMidVec[1:19]),labels=c(1:19))
  p2 <- p2 + scale_y_continuous(limits=c(0,165),breaks=c(seq(0,160,30)),labels=c(seq(0,160,30)))
  p2 <- p2 + geom_hline(yintercept=sig2,linetype=1, col='red', lwd=0.5)
  p2 <- p2 + annotate("text", x = (bpMidVec[1]), y = 160, label=c("Intercross"),size=3.5)
  p2 <- p2 + xlab('Chromosome') + ylab('-log (p-value)')+theme_zg()
  
  tiff("./Figures/figure2.tiff",width=20,height=7, units = "cm",res=300)
  multiplot(p1,p2,cols=2)
  dev.off()
}

H2 <- function(y,SNP){
  IDX<-which(is.na(SNP))
  if(length(IDX) > 0){
    maf1<-mean(SNP[-IDX])/2
    SNP[IDX]<-round(2 * maf1)
  }
  a<-SNP
  asnp <- a-1
  dsnp <- 1-abs(asnp)
  vs <- c()
  va <- c()
  vd <- c()
  for(i in 1:24){
    gsnp <- as.data.frame(cbind(phen=y[,i],asnp,dsnp))
    
    fmla <- as.formula(paste("phen ~ ", paste(colnames(gsnp)[-1], collapse= "+")))
    
    glm.ad <- glm(fmla, family = gaussian(),data=gsnp)
    
    V <- anova(glm.ad)
    snpv <- V$Deviance[-1]
    vs <- c(vs,V$`Resid. Dev`[1])
    va <- c(va,snpv[1])
    vd <- c(vd,snpv[2])
  }
  
  return(list(vs=vs,va=va,vd=vd))
  
}

Figure3<-function(dat,res){
  library(pheatmap)
  test <- read.csv("testcross.csv")[,-1]
  inter <- read.csv("intercross.csv")[,-1]
  
  testloci <- test[,1]
  interloci <- inter[,1]
  
  testloci.chr <- test[,3]
  interloci.chr <- inter[,3]
  
  snpnames <- dat$snpinfo[,3]
  testindex <- snpnames[testloci]
  interindex <- snpnames[interloci]
  
  parin <- res[c(testloci,interloci),8:43]
  parin <- cbind(parin,c(testloci.chr,interloci.chr),c(testloci,interloci))
  rownames(parin) <- c(testindex,interindex)
  
  h2.h<-c()
  h2.d<-c()
  for (i in parin[,38]) {
    SNP <- dat$snps[i,]
    y1 <- dat$pheno_H
    dd1 <- H2(y1,SNP)
    h2.h <- rbind(h2.h,cumsum(dd1$va+dd1$vd)/cumsum(dd1$vs)*100)
    
    y2 <- dat$pheno_D
    dd2 <- H2(y2,SNP)
    h2.d <- rbind(h2.d,cumsum(dd2$va+dd2$vd)/cumsum(dd2$vs)*100)
  }
  
  colnames(h2.h) = paste("yr", 1:24, sep = "")
  rownames(h2.h) = paste("SNP", parin[,38], sep = "")
  
  annotation_col = data.frame(phase = factor(rep(c("phase1", "phase2"), c(6, 18))))
  rownames(annotation_col) = paste("yr", 1:24, sep = "")
  pheatmap(h2.h,treeheight_row=0,treeheight_col=0,
           main = "A", clustering_method = "centroid",
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "canberra",
           annotation_col = annotation_col,
           border_color = "black",angle_col=45,fontsize = 16,
           fontsize_col=13,fontsize_row=10,filename = "./Figures/Figure3-1.tiff",
           width=8,height=10,annotation_legend=FALSE)
  
  
  colnames(h2.d) = paste("yr", 1:24, sep = "")
  rownames(h2.d) = rownames(h2.h) = paste("SNP", parin[,38], sep = "")
  
  annotation_col = data.frame(phase = factor(rep(c("phase1", "phase2"), c(12, 12))))
  rownames(annotation_col) = paste("yr", 1:24, sep = "")
  pheatmap(h2.d,treeheight_row=0,treeheight_col=0,
           main = "B", clustering_method = "centroid",
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "canberra",
           annotation_col = annotation_col,
           border_color = "black",angle_col=45,fontsize = 16,
           fontsize_col=13,fontsize_row=10,filename = "./Figures/Figure3-2.tiff",
           width=8,height=10,annotation_legend=FALSE)
  
}


Figure5<-function(dat,ret=res,index=137076,marktype=2,snp.g=c("AG","AA")){
  ph <- dat$pheno_H
  pd <- dat$pheno_D
  times<-dat$sample_times
  SNP <- dat$snps[index,]
  SNP.index <- as.numeric(names(table(SNP)))
  
  if(marktype==3){
    markerpar <- ret[index,8:43]
  }else{
    markerpar <- ret[index,8:31]
  }
  
  HM <- c();DM <- c();
  for(i in 1:length(SNP.index)){
    HM <- c(HM,mean(ph[which(SNP==SNP.index[i]),1]))
    DM <- c(DM,mean(pd[which(SNP==SNP.index[i]),1]))
  }
  
  nt<-seq(1,40,0.1)
  nt1<-seq(1,24,0.1)
  nt2<-seq(1,-5,-0.1)
  
  j=1
  
  SNP.n <- which(SNP==SNP.index[j])
  npar <- markerpar[(12*(j-1)+1):(j*12)]
  lv.par<-npar[-c(1:6)]
  
  fity <- yy(npar,nt,x1=HM[j],x2=DM[j])[1:length(nt)]
  #fit1y <- phasic1(npar[1:3],nt)
  fit2y <- phasic2(lv.par, nt,x1=HM[j],x2=DM[j])[1:length(nt)]
  ind.fit<-com.get_mu1(lv.par[c(1,2,4,5)],nt,x1=HM[j],x2=DM[j])[1:length(nt)]
  #inter.fit<-com.get_mu2(lv.par[c(1,3,4,6)],nt1,x1=HM[j],x2=DM[j])[1:length(nt1)]
  inter.fit<-fit2y-ind.fit
  
  fit2y_<-phasic2(lv.par, nt2,x1=HM[j],x2=DM[j])[1:length(nt2)]
  fity_ <- yy(npar,nt2,x1=HM[j],x2=DM[j])[1:length(nt2)]
  fit1y_ <- phasic1(npar[1:3],nt2)
  ind.fit_<-com.get_mu1(lv.par[c(1,2,4,5)],nt2,x1=HM[j],x2=DM[j])[1:length(nt2)]
  #inter.fit_<-com.get_mu2(lv.par[c(3,6)],nt2,x1=HM[j],x2=DM[j])[1:length(nt2)]
  inter.fit_<-fit2y_ - ind.fit_
  
  R.lv<-R2(par=lv.par,H1=HM[j],D1=DM[j],nt)
  R.lv_<-R2(par=lv.par,H1=HM[j],D1=DM[j],nt2)
  R.lv.ind<-R2.ind(par=lv.par[c(1,2,4,5)],H1=HM[j],D1=DM[j],nt)
  R.lv.ind_<-R2.ind(par=lv.par[c(1,2,4,5)],H1=HM[j],D1=DM[j],nt2)
  R.Logistic.D<-R1(npar[4:6],nt)
  R.Logistic.D_<-R1(npar[4:6],nt2)
  R.Logistic.H<-R1(npar[1:3],nt)
  R.Logistic.H_<-R1(npar[1:3],nt2)
  
  mt2<-nt[which(R.lv$RH==max(R.lv$RH))]
  my<-phasic2(lv.par, seq(1,mt2,0.1),x1=HM[j],x2=DM[j])[length(seq(1,mt2,0.1))]
  t22<-2*0.5*(max(fit2y)-my)/max(R.lv$RH)
  t21<-2*0.5*(my)/max(R.lv$RH)
  a2<-seq(1,(mt2+t22),0.1)
  a2p<-mt2+t22 
  b2p<-mt2-t22 
  
  mt2.ind<-nt[which(R.lv.ind$RH==max(R.lv.ind$RH))]
  my.ind<-com.get_mu1(lv.par[c(1,2,4,5)], seq(1,mt2.ind,0.1),x1=HM[j],x2=DM[j])[length(seq(1,mt2.ind,0.1))]
  a2p.ind<-mt2.ind+(max(ind.fit)-my.ind)/max(R.lv.ind$RH)
  b2p.ind<-mt2.ind-my.ind/max(R.lv.ind$RH)
  
  mt1<-npar[1]/npar[2]
  t11<-2*0.85/npar[2]
  t12<-2*0.85*(exp(1)-1)/npar[2]
  a1p<-mt1-t11 
  a1<-seq(1,a1p,0.1)
  b1p<-mt1+t12 
  b1<-seq(1,b1p,0.1)
  
  tiff("./Figures/figure5.tiff",width=20,height=20,units = "cm",res=300)
  par(mfcol=c(2,2))
  par(mar=c(0.5,1,2,0.2))
  plot(c(0,0), c(2,0),xlim=c(-9,45), ylim=c(-14,108),xlab="",ylab = "",cex.main=1.5,
       type="n", main=snp.g[j],xaxt="n", yaxt="n", xaxs="i", yaxs="i",frame=F);
  mtext("A",side=3, cex=1.8,adj=0)
  rect(-5,0,40,45);
  rect(-5,55,40,100)
  
  text(20,-11,"Age (year)",cex = 1.3)
  text(-7,77,"Height Growth",cex = 1.3,srt=90)
  text(-7,22,"Height rate",cex = 1.3,srt=90)
  
  for (i in seq(0,35,5)) {
    lines(c(i,i),c(100,101))
    text(i,105,i)
  }
  
  polygon(c(nt,rev(nt)),c(fity,rev(fit2y))+65,
          col = rgb(255,193,193,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt,rev(nt)),c(fit2y,rev(ind.fit))+65,
          col = rgb(152,251,152,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt,rev(nt)),c(rep(0,length(nt)),rev(fit2y))+65,
          col = rgb(191,239,255,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt2,rev(nt2)),c(rep(0,length(nt2)),rev(fit2y_))+65,
          col = rgb(191,239,255,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt2,rev(nt2)),c(rep(0,length(nt2)),rev(inter.fit_))+65,
          col = rgb(255,246,143,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt,rev(nt)),c(rep(0,length(nt)),rev(inter.fit))+65,
          col = rgb(255,246,143,alpha=150,max=255),lwd = 0.01, border = NA)
  
  lines(nt1,fity[1:length(nt1)]+65,type = "l",lwd=2)
  lines(seq(24.1,40,0.1),fity[(length(nt1)+1):length(nt)]+65,type = "l",lwd=2,lty=2)
  lines(nt2,fity_ + 65,type = "l",lwd=2,lty=2)
  lines(nt1,fit2y[1:length(nt1)]+65,type="l",lwd=2,col="black")
  lines(seq(24.1,40,0.1),fit2y[(length(nt1)+1):length(nt)]+65,type="l",lwd=2,col="black",lty=2)
  #lines(nt2,fit2y_ + 65,type = "l",lwd=2,lty=2)
  lines(nt1,ind.fit[1:length(nt1)]+65,type="l",lwd=2,col="green3")
  lines(seq(24.1,40,0.1),ind.fit[(length(nt1)+1):length(nt)]+65,type="l",lwd=2,col="green3",lty=2)
  lines(nt1,inter.fit[1:length(nt1)]+65,type="l",lwd=2,col="darkgoldenrod1")
  lines(seq(24.1,40,0.1),inter.fit[(length(nt1)+1):length(nt)]+65,type="l",lwd=2,col="darkgoldenrod1",lty=2)
  lines(nt2,ind.fit_ + 65,type="l",lwd=2,col="green3",lty=2)
  lines(nt2,inter.fit_ + 65,type="l",lwd=2,col="darkgoldenrod1",lty=2)
  
  lines(c(40,40.5),c(max(fity),max(fity))+65)
  text(43,max(fity)+65,round(max(fity),1))
  lines(c(40,40.5),c(max(fit2y),max(fit2y))+65)
  text(43,max(fit2y)+65,round(max(fit2y),1))
  lines(c(40,40.5),c(max(ind.fit),max(ind.fit))+65)
  text(43,max(ind.fit)+65,round(max(ind.fit),1))
  lines(c(40,40.5),c(inter.fit[length(nt)],inter.fit[length(nt)])+65)
  text(42,inter.fit[length(nt)]+65,round(inter.fit[length(nt)],1))
  
  polygon(c(nt,rev(nt)),c(rep(R.lv.ind$RH[length(nt)],length(nt)),rev(R.lv.ind$RH))*18,
          col = rgb(152,251,152,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt2,rev(nt2)),c(rep(R.lv.ind$RH[length(nt)],length(nt2)),rev(R.lv.ind_$RH))*18,
          col = rgb(152,251,152,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt,rev(nt)),c(rep(R.Logistic.H[length(nt)],length(nt)),rev(R.Logistic.H))*18,
          col = rgb(255,193,193,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt,rev(nt)),c(rep(R.lv$RH[length(nt)],length(nt)),rev(R.lv$RH))*18,
          col = rgb(191,239,255,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt2,rev(nt2)),c(rep(R.lv$RH[length(nt)],length(nt2)),rev(R.lv_$RH))*18,
          col = rgb(191,239,255,alpha=150,max=255),lwd = 0.01, border = NA)
  
  lines(nt1,R.lv.ind$RH[1:length(nt1)]*18,type = "l",col="green3", lwd=2)
  lines(seq(24.1,40,0.1),R.lv.ind$RH[(length(nt1)+1):length(nt)]*18,type = "l",col="green3",lwd=2,lty=2)
  lines(nt2,R.lv.ind_$RH*18,type = "l",col="green3", lwd=2,lty=2)
  lines(nt1,R.lv$RH[1:length(nt1)]*18,type = "l",col="black", lwd=2)
  lines(seq(24.1,40,0.1),R.lv$RH[(length(nt1)+1):length(nt)]*18,type = "l",col="black", lwd=2,lty=2)
  lines(nt2,R.lv_$RH*18,type = "l",col="black", lwd=2,lty=2)
  lines(nt1,R.Logistic.H[1:length(nt1)]*18,lwd=2)
  lines(seq(24.1,40,0.1),R.Logistic.H[(length(nt1)+1):length(nt)]*18,lwd=2,lty=2)
  lines(nt2,R.Logistic.H_*18,lwd=2,lty=2)
  
  lines(c(mt2,40),c(max(R.lv$RH*18),max(R.lv$RH*18)),lty=2)
  lines(c(40,40.5),c(max(R.lv$RH*18),max(R.lv$RH*18)))
  text(42.5,max(R.lv$RH*18),round(max(R.lv$RH),1))
  lines(c(mt1,40),c(max(R.Logistic.H*18),max(R.Logistic.H*18)),lty=2)
  lines(c(40,40.5),c(max(R.Logistic.H*18),max(R.Logistic.H*18)))
  text(42.5,max(R.Logistic.H*18),round(max(R.Logistic.H),1))
  lines(c(mt2.ind,40),c(max(R.lv.ind$RH*18),max(R.lv.ind$RH*18)),lty=2)
  lines(c(40,40.5),c(max(R.lv.ind$RH*18),max(R.lv.ind$RH*18)))
  text(42,max(R.lv.ind$RH*18),round(max(R.lv.ind$RH),1))
  
  text(b2p-1,50,round(b2p,1),cex = 1)
  lines(c(b2p,b2p),c(0,47),lty=2,col="blue")
  lines(c(b2p,b2p),c(52,100),lty=2,col="blue")
  text(a2p,50,round(a2p,1),cex = 1)
  lines(c(a2p,a2p),c(0,47),lty=2,col="blue")
  lines(c(a2p,a2p),c(52,100),lty=2,col="blue")
  lines(c(mt2,mt2),c(0,65+fity[length(seq(1,mt2,0.1))]),col="blue")
  lines(c(mt2,mt2),c(0,-2))
  text(mt2-1,-5,round(mt2,1))
  text(a1p,50,1.0,cex = 1)
  lines(c(a1p,a1p),c(0,47),lty=2,col="lightcoral")
  lines(c(a1p,a1p),c(52,100),lty=2,col="lightcoral")
  text(b1p,50,round(b1p,1),cex = 1)
  lines(c(b1p,b1p),c(0,47),lty=2,col="lightcoral")
  lines(c(b1p,b1p),c(52,100),lty=2,col="lightcoral")
  lines(c(mt1,mt1),c(0,65+fity[length(seq(1,mt1,0.1))]),col="lightcoral")
  lines(c(mt1,mt1),c(0,-2))
  text(mt1,-5,round(mt1,1))
  lines(c(a2p.ind,a2p.ind),c(0,47),lty=2,col="green3")
  lines(c(a2p.ind,a2p.ind),c(52,65+ind.fit[length(seq(1,a2p.ind,0.1))]),lty=2,col="green3")
  text(a2p.ind+1,50,round(a2p.ind,1),cex = 1,col="green3")
  lines(c(b2p.ind,b2p.ind),c(0,65+ind.fit[1]),lty=2,col="green3")
  lines(c(mt2.ind,mt2.ind),c(0,65+ind.fit[length(seq(1,mt2.ind,0.1))]),col="green3")
  lines(c(mt2.ind,mt2.ind),c(0,-2))
  text(mt2.ind+1,-5,round(mt2.ind,1),col = "green3")
  #lines(c(b2p.ind,b2p.ind),c(0,-2))
  text(b2p.ind-1,47,round(b2p.ind,1),col = "green3")
  
  arrows(a1p,95,b1p,95,length = 0.1,code = 3)
  text(15,93,round((b1p-a1p),1))
  arrows(a2p,90,b2p,90,length = 0.1,code = 3)
  text(4,87,round((a2p-b2p),1))
  arrows(a1p,84,a2p,84,length = 0.1,code = 3)
  text((a2p-a1p)/2+a1p,81,round((a2p-a1p),1))
  arrows(a2p.ind,60,b2p.ind,60,length = 0.1,code = 3)
  text(3,63,round((a2p.ind-b2p.ind),1))
  
  
  ##############################################################
  fity <- yy(npar,nt,x1=HM[j],x2=DM[j])[(1+length(nt)):(2*length(nt))]
  fit2y <- phasic2(lv.par, nt,x1=HM[j],x2=DM[j])[(1+length(nt)):(2*length(nt))]
  ind.fit<-com.get_mu1(lv.par[c(1,2,4,5)],nt,x1=HM[j],x2=DM[j])[(1+length(nt)):(2*length(nt))]
  #inter.fit<-com.get_mu2(lv.par[c(3,6)],nt1,x1=HM[j],x2=DM[j])[(1+length(nt1)):(2*length(nt1))]
  inter.fit<-fit2y-ind.fit
  
  fit2y_<-phasic2(lv.par, nt2,x1=HM[j],x2=DM[j])[(1+length(nt2)):(2*length(nt2))]
  fity_ <- yy(npar,nt2,x1=HM[j],x2=DM[j])[(1+length(nt2)):(2*length(nt2))]
  ind.fit_<-com.get_mu1(lv.par[c(1,2,4,5)],nt2,x1=HM[j],x2=DM[j])[(1+length(nt2)):(2*length(nt2))]
  #inter.fit_<-com.get_mu2(lv.par[c(3,6)],nt2,x1=HM[j],x2=DM[j])[(1+length(nt2)):(2*length(nt2))]
  inter.fit_<-fit2y_-ind.fit_
  
  mt2<-nt[which(R.lv$RD==max(R.lv$RD))]
  my<-phasic2(lv.par, seq(1,mt2,0.1),x1=HM[j],x2=DM[j])[2*length(seq(1,mt2,0.1))]
  t22<-(max(fit2y)-my)/max(R.lv$RD)
  t21<-(my)/max(R.lv$RD)
  a2<-seq(1,(mt2+t22),0.1)
  a2p<-mt2+t22 
  b2p<-mt2-t21 
  #b2<-seq(1,(mt2-t21),-0.1)
  
  mt2.ind<-nt[which(R.lv.ind$RD==max(R.lv.ind$RD))]
  my.ind<-com.get_mu1(lv.par[c(1,2,4,5)], seq(1,mt2.ind,0.1),x1=HM[j],x2=DM[j])[2*length(seq(1,mt2.ind,0.1))]
  a2p.ind<-mt2.ind+(max(ind.fit)-my.ind)/max(R.lv.ind$RD)
  b2p.ind<-mt2.ind-my.ind/max(R.lv.ind$RD)
  
  mt1<-npar[4]/npar[5]
  t11<-2*0.85/npar[5]
  t12<-2*0.85*(exp(1)-1)/npar[5]
  a1p<-mt1-t11 
  a1<-seq(1,a1p,0.1)
  b1p<-mt1+t12 
  b1<-seq(1,b1p,0.1)
  
  par(mar=c(0.5,1,2,0.2))
  plot(c(0,0), c(2,0),xlim=c(-9,45), ylim=c(-14,108),xlab="",ylab = "",cex.main=1.5,
       type="n", main=snp.g[j],xaxt="n", yaxt="n", xaxs="i", yaxs="i",frame=F);
  mtext("B",side=3, cex=1.8,adj=0)
  rect(-5,0,40,45);
  rect(-5,55,40,100)
  
  text(20,-11,"Age (year)",cex = 1.3)
  text(-7,77,"Diameter Growth",cex = 1.3,srt=90)
  text(-7,22,"Diameter rate",cex = 1.3,srt=90)
  
  for (i in seq(0,35,5)) {
    lines(c(i,i),c(100,101))
    text(i,105,i)
  }
  
  polygon(c(nt,rev(nt)),c(fity,rev(fit2y))+60,
          col = rgb(255,193,193,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt,rev(nt)),c(rep(0,length(nt)),rev(fit2y))+60,
          col = rgb(191,239,255,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt2,rev(nt2)),c(rep(0,length(nt2)),rev(fit2y_))+60,
          col = rgb(191,239,255,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt2,rev(nt2)),c(rep(0,length(nt2)),rev(inter.fit_))+60,
          col = rgb(255,246,143,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt,rev(nt)),c(rep(0,length(nt)),rev(inter.fit))+60,
          col = rgb(255,246,143,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt,rev(nt)),c(fit2y,rev(ind.fit))+60,
          col = rgb(152,251,152,alpha=150,max=255),lwd = 0.01, border = NA)
  
  lines(nt1,fity[1:length(nt1)]+60,type = "l",lwd=2)
  lines(seq(24.1,40,0.1),fity[(length(nt1)+1):length(nt)]+60,type = "l",lwd=2,lty=2)
  lines(nt2,fity_ + 60,type = "l",lwd=2,lty=2)
  lines(nt1,fit2y[1:length(nt1)]+60,type="l",lwd=2,col="black")
  lines(seq(24.1,40,0.1),fit2y[(length(nt1)+1):length(nt)]+60,type="l",lwd=2,col="black",lty=2)
  #lines(nt2,fit2y_ + 60,type = "l",lwd=2,lty=2)
  lines(nt1,ind.fit[1:length(nt1)]+60,type="l",lwd=2,col="green3")
  lines(seq(24.1,40,0.1),ind.fit[(length(nt1)+1):length(nt)]+60,type="l",lwd=2,col="green3",lty=2)
  lines(nt1,inter.fit[1:length(nt1)]+60,type="l",lwd=2,col="darkgoldenrod1")
  lines(seq(24.1,40,0.1),inter.fit[(length(nt1)+1):length(nt)]+60,type="l",lwd=2,col="darkgoldenrod1",lty=2)
  lines(nt2,ind.fit_ + 60,type="l",lwd=2,col="green3",lty=2)
  lines(nt2,inter.fit_ + 60,type="l",lwd=2,col="darkgoldenrod1",lty=2)
  
  lines(c(40,40.5),c(max(fity),max(fity))+60)
  text(43,max(fity)+60,round(max(fity),1))
  lines(c(40,40.5),c(max(fit2y),max(fit2y))+60)
  text(43,max(fit2y)+62,round(max(fit2y),1))
  lines(c(40,40.5),c(max(ind.fit),max(ind.fit))+60)
  text(43,max(ind.fit)+59,round(max(ind.fit),1))
  lines(c(40,40.5),c(inter.fit[length(nt)],inter.fit[length(nt)])+60)
  text(42,inter.fit[length(nt)]+60,round(inter.fit[length(nt)],1))
  
  polygon(c(nt,rev(nt)),c(rep(0,length(nt)),rev(R.Logistic.D))*16,
          col = rgb(255,193,193,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt,rev(nt)),c(rep(R.lv$RD[length(nt)],length(nt)),rev(R.lv$RD))*16,
          col = rgb(191,239,255,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt2,rev(nt2)),c(rep(R.lv$RD[length(nt)],length(nt2)),rev(R.lv_$RD))*16,
          col = rgb(191,239,255,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt,rev(nt)),c(rep(R.lv.ind$RD[length(nt)],length(nt)),rev(R.lv.ind$RD))*16,
          col = rgb(152,251,152,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt2,rev(nt2)),c(rep(R.lv.ind$RD[length(nt)],length(nt2)),rev(R.lv.ind_$RD))*16,
          col = rgb(152,251,152,alpha=150,max=255),lwd = 0.01, border = NA)
  
  lines(nt1,R.lv.ind$RD[1:length(nt1)]*16,type = "l",col="green3", lwd=2)
  lines(seq(24.1,40,0.1),R.lv.ind$RD[(length(nt1)+1):length(nt)]*16,type = "l",col="green3",lwd=2,lty=2)
  lines(nt2,R.lv.ind_$RD*16,type = "l",col="green3", lwd=2,lty=2)
  lines(nt1,R.lv$RD[1:length(nt1)]*16,type = "l",col="black", lwd=2)
  lines(seq(24.1,40,0.1),R.lv$RD[(length(nt1)+1):length(nt)]*16,type = "l",col="black", lwd=2,lty=2)
  lines(nt2,R.lv_$RD*16,type = "l",col="black", lwd=2,lty=2)
  lines(nt1,R.Logistic.D[1:length(nt1)]*16,lwd=2)
  lines(seq(24.1,40,0.1),R.Logistic.D[(length(nt1)+1):length(nt)]*16,lwd=2,lty=2)
  lines(nt2,R.Logistic.D_*16,lwd=2,lty=2)
  
  lines(c(mt2,40),c(max(R.lv$RD*16),max(R.lv$RD*16)),lty=2)
  lines(c(40,40.5),c(max(R.lv$RD*16),max(R.lv$RD*16)))
  text(42.5,max(R.lv$RD*16),round(max(R.lv$RD),1))
  lines(c(mt1,40),c(max(R.Logistic.D*16),max(R.Logistic.D*16)),lty=2)
  lines(c(40,40.5),c(max(R.Logistic.D*16),max(R.Logistic.D*16)))
  text(42.5,max(R.Logistic.D*16),round(max(R.Logistic.D),1))
  lines(c(mt2.ind,40),c(max(R.lv.ind$RD*16),max(R.lv.ind$RD*16)),lty=2)
  lines(c(40,40.5),c(max(R.lv.ind$RD*16),max(R.lv.ind$RD*16)))
  text(42.5,max(R.lv.ind$RD*16),round(max(R.lv.ind$RD),1))
  
  text(b2p,50,round(b2p,1),cex = 1)
  lines(c(b2p,b2p),c(0,47),lty=2,col="blue")
  lines(c(b2p,b2p),c(52,100),lty=2,col="blue")
  text(a2p+1,50,round(a2p,1),cex = 1)
  lines(c(a2p,a2p),c(0,47),lty=2,col="blue")
  lines(c(a2p,a2p),c(52,100),lty=2,col="blue")
  lines(c(mt2,mt2),c(0,60+fity[length(seq(1,mt2,0.1))]),col="blue")
  lines(c(mt2,mt2),c(0,-2))
  text(mt2,-5,round(mt2,1))
  text(a1p-2,50,round(a1p,1),cex = 1)
  lines(c(a1p-1,a1p-1),c(0,47),lty=2,col="lightcoral")
  lines(c(a1p-1,a1p-1),c(52,100),lty=2,col="lightcoral")
  text(b1p,50,round(b1p,1),cex = 1)
  lines(c(b1p,b1p),c(0,47),lty=2,col="lightcoral")
  lines(c(b1p,b1p),c(52,100),lty=2,col="lightcoral")
  lines(c(mt1,mt1),c(0,60+fity[length(seq(1,mt1,0.1))]),col="lightcoral")
  lines(c(mt1,mt1),c(0,-2))
  text(mt1,-5,round(mt1,1))
  lines(c(a2p.ind,a2p.ind),c(0,60+ind.fit[length(seq(1,a2p.ind,0.1))]),lty=2,col="green3")
  lines(c(a2p.ind,a2p.ind),c(0,-2))
  text(a2p.ind,-5,round(a2p.ind,1),cex = 1,col = "green3")
  lines(c(b2p.ind,b2p.ind),c(0,60+ind.fit[1]),lty=2,col="green3")
  lines(c(mt2.ind,mt2.ind),c(0,60+ind.fit[length(seq(1,mt2.ind,0.1))]),col="green3")
  lines(c(mt2.ind,mt2.ind),c(0,-2))
  text(mt2.ind-1,-7.5,round(mt2.ind,1),col = "green3")
  lines(c(b2p.ind,b2p.ind),c(0,-2))
  text(b2p.ind,-5,round(b2p.ind,1),col = "green3")
  
  arrows(a1p-1,95,b1p,95,length = 0.1,code = 3)
  text(20,93,round((b1p-a1p),1))
  arrows(a2p,90,b2p,90,length = 0.1,code = 3)
  text(5,87,round((a2p-b2p),1))
  arrows(a1p-1,84,a2p,84,length = 0.1,code = 3)
  text((a2p-a1p)/2+a1p,82,round((a2p-a1p),1))
  arrows(a2p.ind,63,b2p.ind,63,length = 0.1,code = 3)
  text(6,65,round((a2p.ind-b2p.ind),1))

  ###################################################
  j=2
  
  SNP.n <- which(SNP==SNP.index[j])
  npar <- markerpar[(12*(j-1)+1):(j*12)]
  lv.par<-npar[-c(1:6)]
  
  fity <- yy(npar,nt,x1=HM[j],x2=DM[j])[1:length(nt)]
  #fit1y <- phasic1(npar[1:3],nt)
  fit2y <- phasic2(lv.par, nt,x1=HM[j],x2=DM[j])[1:length(nt)]
  ind.fit<-com.get_mu1(lv.par[c(1,2,4,5)],nt,x1=HM[j],x2=DM[j])[1:length(nt)]
  #inter.fit<-com.get_mu2(lv.par[c(1,3,4,6)],nt1,x1=HM[j],x2=DM[j])[1:length(nt1)]
  inter.fit<-fit2y-ind.fit
  
  fit2y_<-phasic2(lv.par, nt2,x1=HM[j],x2=DM[j])[1:length(nt2)]
  fity_ <- yy(npar,nt2,x1=HM[j],x2=DM[j])[1:length(nt2)]
  fit1y_ <- phasic1(npar[1:3],nt2)
  ind.fit_<-com.get_mu1(lv.par[c(1,2,4,5)],nt2,x1=HM[j],x2=DM[j])[1:length(nt2)]
  #inter.fit_<-com.get_mu2(lv.par[c(3,6)],nt2,x1=HM[j],x2=DM[j])[1:length(nt2)]
  inter.fit_<-fit2y_ - ind.fit_
  
  R.lv<-R2(par=lv.par,H1=HM[j],D1=DM[j],nt)
  R.lv_<-R2(par=lv.par,H1=HM[j],D1=DM[j],nt2)
  R.lv.ind<-R2.ind(par=lv.par[c(1,2,4,5)],H1=HM[j],D1=DM[j],nt)
  R.lv.ind_<-R2.ind(par=lv.par[c(1,2,4,5)],H1=HM[j],D1=DM[j],nt2)
  R.Logistic.D<-R1(npar[4:6],nt)
  R.Logistic.D_<-R1(npar[4:6],nt2)
  R.Logistic.H<-R1(npar[1:3],nt)
  R.Logistic.H_<-R1(npar[1:3],nt2)
  
  mt2<-nt[which(R.lv$RH==max(R.lv$RH))]
  my<-phasic2(lv.par, seq(1,mt2,0.1),x1=HM[j],x2=DM[j])[length(seq(1,mt2,0.1))]
  t22<-2*0.5*(max(fit2y)-my)/max(R.lv$RH)
  t21<-2*0.5*(my)/max(R.lv$RH)
  a2<-seq(1,(mt2+t22),0.1)
  a2p<-mt2+t22 
  b2p<-mt2-t22 
  
  mt2.ind<-nt[which(R.lv.ind$RH==max(R.lv.ind$RH))]
  my.ind<-com.get_mu1(lv.par[c(1,2,4,5)], seq(1,mt2.ind,0.1),x1=HM[j],x2=DM[j])[length(seq(1,mt2.ind,0.1))]
  a2p.ind<-mt2.ind+(max(ind.fit)-my.ind)/max(R.lv.ind$RH)
  b2p.ind<-mt2.ind-my.ind/max(R.lv.ind$RH)
  
  mt1<-npar[1]/npar[2]
  t11<-2*0.85/npar[2]
  t12<-2*0.85*(exp(1)-1)/npar[2]
  a1p<-mt1-t11 
  #a1<-seq(1,a1p,0.1)
  b1p<-mt1+t12 
  b1<-seq(1,b1p,0.1)
  
  par(mar=c(0.5,1,2,0.2))
  plot(c(0,0), c(2,0),xlim=c(-9,45), ylim=c(-14,108),xlab="",ylab = "",cex.main=1.5,
       type="n", main=snp.g[j],xaxt="n", yaxt="n", xaxs="i", yaxs="i",frame=F);
  rect(-5,0,40,45);
  rect(-5,55,40,100)
  
  text(20,-11,"Age (year)",cex = 1.3)
  text(-7,77,"Height Growth",cex = 1.3,srt=90)
  text(-7,22,"Height rate",cex = 1.3,srt=90)
  
  for (i in seq(0,35,5)) {
    lines(c(i,i),c(100,101))
    text(i,105,i)
  }
  
  polygon(c(nt,rev(nt)),c(fity,rev(fit2y))+65,
          col = rgb(255,193,193,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt,rev(nt)),c(fit2y,rev(ind.fit))+65,
          col = rgb(152,251,152,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt,rev(nt)),c(rep(0,length(nt)),rev(fit2y))+65,
          col = rgb(191,239,255,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt2,rev(nt2)),c(rep(0,length(nt2)),rev(fit2y_))+65,
          col = rgb(191,239,255,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt2,rev(nt2)),c(rep(0,length(nt2)),rev(inter.fit_))+65,
          col = rgb(255,246,143,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt,rev(nt)),c(rep(0,length(nt)),rev(inter.fit))+65,
          col = rgb(255,246,143,alpha=150,max=255),lwd = 0.01, border = NA)
  
  lines(nt1,fity[1:length(nt1)]+65,type = "l",lwd=2)
  lines(seq(24.1,40,0.1),fity[(length(nt1)+1):length(nt)]+65,type = "l",lwd=2,lty=2)
  lines(nt2,fity_ + 65,type = "l",lwd=2,lty=2)
  lines(nt1,fit2y[1:length(nt1)]+65,type="l",lwd=2,col="black")
  lines(seq(24.1,40,0.1),fit2y[(length(nt1)+1):length(nt)]+65,type="l",lwd=2,col="black",lty=2)
  #lines(nt2,fit2y_ + 65,type = "l",lwd=2,lty=2)
  lines(nt1,ind.fit[1:length(nt1)]+65,type="l",lwd=2,col="green3")
  lines(seq(24.1,40,0.1),ind.fit[(length(nt1)+1):length(nt)]+65,type="l",lwd=2,col="green3",lty=2)
  lines(nt1,inter.fit[1:length(nt1)]+65,type="l",lwd=2,col="darkgoldenrod1")
  lines(seq(24.1,40,0.1),inter.fit[(length(nt1)+1):length(nt)]+65,type="l",lwd=2,col="darkgoldenrod1",lty=2)
  lines(nt2,ind.fit_ + 65,type="l",lwd=2,col="green3",lty=2)
  lines(nt2,inter.fit_ + 65,type="l",lwd=2,col="darkgoldenrod1",lty=2)
  
  lines(c(40,40.5),c(max(fity),max(fity))+65)
  text(43,max(fity)+65,round(max(fity),1))
  lines(c(40,40.5),c(max(fit2y),max(fit2y))+65)
  text(43,max(fit2y)+65,round(max(fit2y),1))
  lines(c(40,40.5),c(max(ind.fit),max(ind.fit))+65)
  text(43,max(ind.fit)+65,round(max(ind.fit),1))
  lines(c(40,40.5),c(inter.fit[length(nt)],inter.fit[length(nt)])+65)
  text(42,inter.fit[length(nt)]+65,round(inter.fit[length(nt)],1))
  
  polygon(c(nt,rev(nt)),c(rep(R.lv.ind$RH[length(nt)],length(nt)),rev(R.lv.ind$RH))*18,
          col = rgb(152,251,152,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt2,rev(nt2)),c(rep(R.lv.ind$RH[length(nt)],length(nt2)),rev(R.lv.ind_$RH))*18,
          col = rgb(152,251,152,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt,rev(nt)),c(rep(R.Logistic.H[length(nt)],length(nt)),rev(R.Logistic.H))*18,
          col = rgb(255,193,193,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt,rev(nt)),c(rep(R.lv$RH[length(nt)],length(nt)),rev(R.lv$RH))*18,
          col = rgb(191,239,255,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt2,rev(nt2)),c(rep(R.lv$RH[length(nt)],length(nt2)),rev(R.lv_$RH))*18,
          col = rgb(191,239,255,alpha=150,max=255),lwd = 0.01, border = NA)
  
  lines(nt1,R.lv.ind$RH[1:length(nt1)]*18,type = "l",col="green3", lwd=2)
  lines(seq(24.1,40,0.1),R.lv.ind$RH[(length(nt1)+1):length(nt)]*18,type = "l",col="green3",lwd=2,lty=2)
  lines(nt2,R.lv.ind_$RH*18,type = "l",col="green3", lwd=2,lty=2)
  lines(nt1,R.lv$RH[1:length(nt1)]*18,type = "l",col="black", lwd=2)
  lines(seq(24.1,40,0.1),R.lv$RH[(length(nt1)+1):length(nt)]*18,type = "l",col="black", lwd=2,lty=2)
  lines(nt2,R.lv_$RH*18,type = "l",col="black", lwd=2,lty=2)
  lines(nt1,R.Logistic.H[1:length(nt1)]*18,lwd=2)
  lines(seq(24.1,40,0.1),R.Logistic.H[(length(nt1)+1):length(nt)]*18,lwd=2,lty=2)
  lines(nt2,R.Logistic.H_*18,lwd=2,lty=2)
  
  lines(c(mt2,40),c(max(R.lv$RH*18),max(R.lv$RH*18)),lty=2)
  lines(c(40,40.5),c(max(R.lv$RH*18),max(R.lv$RH*18)))
  text(42.5,max(R.lv$RH*18),round(max(R.lv$RH),1))
  lines(c(mt1,40),c(max(R.Logistic.H*18),max(R.Logistic.H*18)),lty=2)
  lines(c(40,40.5),c(max(R.Logistic.H*18),max(R.Logistic.H*18)))
  text(42.5,max(R.Logistic.H*18),round(max(R.Logistic.H),1))
  lines(c(mt2.ind,40),c(max(R.lv.ind$RH*18),max(R.lv.ind$RH*18)),lty=2)
  lines(c(40,40.5),c(max(R.lv.ind$RH*18),max(R.lv.ind$RH*18)))
  text(42.5,max(R.lv.ind$RH*18),round(max(R.lv.ind$RH),1))
  
  text(b2p-1,50,round(b2p,1),cex = 1)
  lines(c(b2p,b2p),c(0,47),lty=2,col="blue")
  lines(c(b2p,b2p),c(52,100),lty=2,col="blue")
  text(a2p,50,round(a2p,1),cex = 1)
  lines(c(a2p,a2p),c(0,47),lty=2,col="blue")
  lines(c(a2p,a2p),c(52,100),lty=2,col="blue")
  lines(c(mt2,mt2),c(0,65+fity[length(seq(1,mt2,0.1))]),col="blue")
  lines(c(mt2,mt2),c(0,-2))
  text(mt2-1,-5,round(mt2,1))
  text(a1p+1,50,round(a1p,1),cex = 1)
  lines(c(a1p+0.5,a1p+0.5),c(0,47),lty=2,col="lightcoral")
  lines(c(a1p+0.5,a1p+0.5),c(52,100),lty=2,col="lightcoral")
  text(b1p,50,round(b1p,1),cex = 1)
  lines(c(b1p,b1p),c(0,47),lty=2,col="lightcoral")
  lines(c(b1p,b1p),c(52,100),lty=2,col="lightcoral")
  lines(c(mt1,mt1),c(0,65+fity[length(seq(1,mt1,0.1))]),col="lightcoral")
  lines(c(mt1,mt1),c(0,-2))
  text(mt1,-5,round(mt1,1))
  lines(c(a2p.ind,a2p.ind),c(0,65+ind.fit[length(seq(1,a2p.ind,0.1))]),lty=2,col="green3")
  text(a2p.ind,-7.5,round(a2p.ind,1),cex = 1,col = "green3")
  lines(c(b2p.ind,b2p.ind),c(0,65+ind.fit[1]),lty=2,col="green3")
  lines(c(mt2.ind,mt2.ind),c(0,65+ind.fit[length(seq(1,mt2.ind,0.1))]),col="green3")
  lines(c(mt2.ind,mt2.ind),c(0,-2))
  text(mt2.ind+1,-5,round(mt2.ind,1),col = "green3")
  lines(c(b2p.ind,b2p.ind),c(0,-2))
  text(b2p.ind-0.5,-7.5,round(b2p.ind,1),col = "green3")
  
  arrows(a1p+0.5,95,b1p,95,length = 0.1,code = 3)
  text(15,93,round((b1p-a1p),1))
  arrows(a2p,90,b2p,90,length = 0.1,code = 3)
  text(3,87,round((a2p-b2p),1))
  arrows(a1p+0.5,84,a2p,84,length = 0.1,code = 3)
  text((a2p-a1p)/2+a1p+0.5,81,round((a2p-a1p),1))
  arrows(a2p.ind,60,b2p.ind,60,length = 0.1,code = 3)
  text(3,63,round((a2p.ind-b2p.ind),1))

  ##############################################################
  fity <- yy(npar,nt,x1=HM[j],x2=DM[j])[(1+length(nt)):(2*length(nt))]
  fit2y <- phasic2(lv.par, nt,x1=HM[j],x2=DM[j])[(1+length(nt)):(2*length(nt))]
  ind.fit<-com.get_mu1(lv.par[c(1,2,4,5)],nt,x1=HM[j],x2=DM[j])[(1+length(nt)):(2*length(nt))]
  #inter.fit<-com.get_mu2(lv.par[c(3,6)],nt1,x1=HM[j],x2=DM[j])[(1+length(nt1)):(2*length(nt1))]
  inter.fit<-fit2y-ind.fit
  
  fit2y_<-phasic2(lv.par, nt2,x1=HM[j],x2=DM[j])[(1+length(nt2)):(2*length(nt2))]
  fity_ <- yy(npar,nt2,x1=HM[j],x2=DM[j])[(1+length(nt2)):(2*length(nt2))]
  ind.fit_<-com.get_mu1(lv.par[c(1,2,4,5)],nt2,x1=HM[j],x2=DM[j])[(1+length(nt2)):(2*length(nt2))]
  #inter.fit_<-com.get_mu2(lv.par[c(3,6)],nt2,x1=HM[j],x2=DM[j])[(1+length(nt2)):(2*length(nt2))]
  inter.fit_<-fit2y_-ind.fit_
  
  mt2<-nt[which(R.lv$RD==max(R.lv$RD))]
  my<-phasic2(lv.par, seq(1,mt2,0.1),x1=HM[j],x2=DM[j])[2*length(seq(1,mt2,0.1))]
  t22<-(max(fit2y)-my)/max(R.lv$RD)
  t21<-(my)/max(R.lv$RD)
  a2<-seq(1,(mt2+t22),0.1)
  a2p<-mt2+t22 
  b2p<-mt2-t21 
  #b2<-seq(1,(mt2-t21),-0.1)
  
  mt2.ind<-nt[which(R.lv.ind$RD==max(R.lv.ind$RD))]
  my.ind<-com.get_mu1(lv.par[c(1,2,4,5)], seq(1,mt2.ind,0.1),x1=HM[j],x2=DM[j])[2*length(seq(1,mt2.ind,0.1))]
  a2p.ind<-mt2.ind+(max(ind.fit)-my.ind)/max(R.lv.ind$RD)
  b2p.ind<-mt2.ind-my.ind/max(R.lv.ind$RD)
  
  mt1<-npar[4]/npar[5]
  t11<-2*0.85/npar[5]
  t12<-2*0.85*(exp(1)-1)/npar[5]
  a1p<-mt1-t11 
  a1<-seq(1,a1p,0.1)
  b1p<-mt1+t12 
  b1<-seq(1,b1p,0.1)
  
  par(mar=c(0.5,1,2,0.2))
  plot(c(0,0), c(2,0),xlim=c(-9,45), ylim=c(-14,108),xlab="",ylab = "",cex.main=1.5,
       type="n", main=snp.g[j],xaxt="n", yaxt="n", xaxs="i", yaxs="i",frame=F);
  rect(-5,0,40,45);
  rect(-5,55,40,100)
  
  text(20,-11,"Age (year)",cex = 1.3)
  text(-7,77,"Diameter Growth",cex = 1.3,srt=90)
  text(-7,22,"Diameter rate",cex = 1.3,srt=90)
  
  for (i in seq(0,35,5)) {
    lines(c(i,i),c(100,101))
    text(i,105,i)
  }
  
  polygon(c(nt,rev(nt)),c(fity,rev(fit2y))+60,
          col = rgb(255,193,193,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt,rev(nt)),c(rep(0,length(nt)),rev(fit2y))+60,
          col = rgb(191,239,255,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt2,rev(nt2)),c(rep(0,length(nt2)),rev(fit2y_))+60,
          col = rgb(191,239,255,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt2,rev(nt2)),c(rep(0,length(nt2)),rev(inter.fit_))+60,
          col = rgb(255,246,143,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt,rev(nt)),c(rep(0,length(nt)),rev(inter.fit))+60,
          col = rgb(255,246,143,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt,rev(nt)),c(fit2y,rev(ind.fit))+60,
          col = rgb(152,251,152,alpha=150,max=255),lwd = 0.01, border = NA)
  
  lines(nt1,fity[1:length(nt1)]+60,type = "l",lwd=2)
  lines(seq(24.1,40,0.1),fity[(length(nt1)+1):length(nt)]+60,type = "l",lwd=2,lty=2)
  lines(nt2,fity_ + 60,type = "l",lwd=2,lty=2)
  lines(nt1,fit2y[1:length(nt1)]+60,type="l",lwd=2,col="black")
  lines(seq(24.1,40,0.1),fit2y[(length(nt1)+1):length(nt)]+60,type="l",lwd=2,col="black",lty=2)
  #lines(nt2,fit2y_ + 60,type = "l",lwd=2,lty=2)
  lines(nt1,ind.fit[1:length(nt1)]+60,type="l",lwd=2,col="green3")
  lines(seq(24.1,40,0.1),ind.fit[(length(nt1)+1):length(nt)]+60,type="l",lwd=2,col="green3",lty=2)
  lines(nt1,inter.fit[1:length(nt1)]+60,type="l",lwd=2,col="darkgoldenrod1")
  lines(seq(24.1,40,0.1),inter.fit[(length(nt1)+1):length(nt)]+60,type="l",lwd=2,col="darkgoldenrod1",lty=2)
  lines(nt2,ind.fit_ + 60,type="l",lwd=2,col="green3",lty=2)
  lines(nt2,inter.fit_ + 60,type="l",lwd=2,col="darkgoldenrod1",lty=2)
  
  lines(c(40,40.5),c(max(fity),max(fity))+60)
  text(43,max(fity)+60,round(max(fity),1))
  lines(c(40,40.5),c(max(fit2y),max(fit2y))+60)
  text(43,max(fit2y)+61,round(max(fit2y),1))
  lines(c(40,40.5),c(max(ind.fit),max(ind.fit))+60)
  text(43,max(ind.fit)+60,round(max(ind.fit),1))
  lines(c(40,40.5),c(inter.fit[length(nt)],inter.fit[length(nt)])+60)
  text(43,inter.fit[length(nt)]+60,round(inter.fit[length(nt)],1))
  
  polygon(c(nt,rev(nt)),c(rep(0,length(nt)),rev(R.Logistic.D))*12,
          col = rgb(255,193,193,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt,rev(nt)),c(rep(R.lv$RD[length(nt)],length(nt)),rev(R.lv$RD))*12,
          col = rgb(191,239,255,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt2,rev(nt2)),c(rep(R.lv$RD[length(nt)],length(nt2)),rev(R.lv_$RD))*12,
          col = rgb(191,239,255,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt,rev(nt)),c(rep(R.lv.ind$RD[length(nt)],length(nt)),rev(R.lv.ind$RD))*12,
          col = rgb(152,251,152,alpha=150,max=255),lwd = 0.01, border = NA)
  polygon(c(nt2,rev(nt2)),c(rep(R.lv.ind$RD[length(nt)],length(nt2)),rev(R.lv.ind_$RD))*12,
          col = rgb(152,251,152,alpha=150,max=255),lwd = 0.01, border = NA)
  
  lines(nt1,R.lv.ind$RD[1:length(nt1)]*12,type = "l",col="green3", lwd=2)
  lines(seq(24.1,40,0.1),R.lv.ind$RD[(length(nt1)+1):length(nt)]*12,type = "l",col="green3",lwd=2,lty=2)
  lines(nt2,R.lv.ind_$RD*12,type = "l",col="green3", lwd=2,lty=2)
  lines(nt1,R.lv$RD[1:length(nt1)]*12,type = "l",col="black", lwd=2)
  lines(seq(24.1,40,0.1),R.lv$RD[(length(nt1)+1):length(nt)]*12,type = "l",col="black", lwd=2,lty=2)
  lines(nt2,R.lv_$RD*12,type = "l",col="black", lwd=2,lty=2)
  lines(nt1,R.Logistic.D[1:length(nt1)]*12,lwd=2)
  lines(seq(24.1,40,0.1),R.Logistic.D[(length(nt1)+1):length(nt)]*12,lwd=2,lty=2)
  lines(nt2,R.Logistic.D_*12,lwd=2,lty=2)
  
  lines(c(mt2,40),c(max(R.lv$RD*12),max(R.lv$RD*12)),lty=2)
  lines(c(40,40.5),c(max(R.lv$RD*12),max(R.lv$RD*12)))
  text(42.5,max(R.lv$RD*12),round(max(R.lv$RD),1))
  lines(c(mt1,40),c(max(R.Logistic.D*12),max(R.Logistic.D*12)),lty=2)
  lines(c(40,40.5),c(max(R.Logistic.D*12),max(R.Logistic.D*12)))
  text(42.5,max(R.Logistic.D*12),round(max(R.Logistic.D),1))
  lines(c(mt2.ind,40),c(max(R.lv.ind$RD*12),max(R.lv.ind$RD*12)),lty=2)
  lines(c(40,40.5),c(max(R.lv.ind$RD*12),max(R.lv.ind$RD*12)))
  text(42.5,max(R.lv.ind$RD*12),round(max(R.lv.ind$RD),1))
  
  text(b2p,50,round(b2p,1),cex = 1)
  lines(c(b2p,b2p),c(0,47),lty=2,col="blue")
  lines(c(b2p,b2p),c(52,100),lty=2,col="blue")
  text(a2p+1,50,round(a2p,1),cex = 1)
  lines(c(a2p,a2p),c(0,47),lty=2,col="blue")
  lines(c(a2p,a2p),c(52,100),lty=2,col="blue")
  lines(c(mt2,mt2),c(0,60+fity[length(seq(1,mt2,0.1))]),col="blue")
  lines(c(mt2,mt2),c(0,-2))
  text(mt2,-5,round(mt2,1))
  text(a1p-1,50,round(a1p,1),cex = 1)
  lines(c(a1p-0.5,a1p-0.5),c(0,47),lty=2,col="lightcoral")
  lines(c(a1p-0.5,a1p-0.5),c(52,100),lty=2,col="lightcoral")
  text(b1p,50,round(b1p,1),cex = 1)
  lines(c(b1p,b1p),c(0,47),lty=2,col="lightcoral")
  lines(c(b1p,b1p),c(52,100),lty=2,col="lightcoral")
  lines(c(mt1,mt1),c(0,60+fity[length(seq(1,mt1,0.1))]),col="lightcoral")
  lines(c(mt1,mt1),c(0,-2))
  text(mt1,-5,round(mt1,1))
  lines(c(a2p.ind,a2p.ind),c(0,60+ind.fit[length(seq(1,a2p.ind,0.1))]),lty=2,col="green3")
  lines(c(a2p.ind,a2p.ind),c(0,-2))
  text(a2p.ind,-5,round(a2p.ind,1),cex = 1,col="green3")
  lines(c(b2p.ind,b2p.ind),c(0,60+ind.fit[1]),lty=2,col="green3")
  lines(c(mt2.ind,mt2.ind),c(0,60+ind.fit[length(seq(1,mt2.ind,0.1))]),col="green3")
  lines(c(mt2.ind,mt2.ind),c(0,-2))
  text(mt2.ind-1,-7.5,round(mt2.ind,1),col = "green3")
  lines(c(b2p.ind,b2p.ind),c(0,-2))
  text(b2p.ind,-5,round(b2p.ind,1),col = "green3")
  
  arrows(a1p-0.5,95,b1p,95,length = 0.1,code = 3)
  text(18,93,round((b1p-a1p),1))
  arrows(a2p,90,b2p,90,length = 0.1,code = 3)
  text(4.5,93,round((a2p-b2p),1))
  arrows(a1p,84,a2p,84,length = 0.1,code = 3)
  text(5.5,84,round((a2p-a1p),1))
  arrows(a2p.ind,60,b2p.ind,60,length = 0.1,code = 3)
  text(4.5,58,round((a2p.ind-b2p.ind),1))
  dev.off()
  
  
  
}

H2 <- function(y,SNP){
  IDX<-which(is.na(SNP))
  if(length(IDX) > 0){
    maf1<-mean(SNP[-IDX])/2
    SNP[IDX]<-round(2 * maf1)
  }
  a<-SNP
  asnp <- a-1
  dsnp <- 1-abs(asnp)
  vs <- c()
  va <- c()
  vd <- c()
  for(i in 1:24){
    gsnp <- as.data.frame(cbind(phen=y[,i],asnp,dsnp))
    
    fmla <- as.formula(paste("phen ~ ", paste(colnames(gsnp)[-1], collapse= "+")))
    
    glm.ad <- glm(fmla, family = gaussian(),data=gsnp)
    
    V <- anova(glm.ad)
    snpv <- V$Deviance[-1]
    vs <- c(vs,V$`Resid. Dev`[1])
    va <- c(va,snpv[1])
    vd <- c(vd,snpv[2])
  }
  
  return(list(vs=vs,va=va,vd=vd))
  
}


Figure6 <- function(dat,ret,index=137076){
  
  testpar<-ret[index,8:31]
  nt <- seq(1,24,0.1)
  
  ph <- dat$pheno_H
  pd <- dat$pheno_D
  
  SNP <- dat$snps[index,]
  SNP.index <- as.numeric(names(table(SNP)))
  HM <- c();DM <- c();
  for(i in 1:length(SNP.index)){
    HM <- c(HM,mean(ph[which(SNP==SNP.index[i]),1]))
    DM <- c(DM,mean(pd[which(SNP==SNP.index[i]),1]))
  }
  npar1 <- testpar[1:12]
  npar2 <- testpar[13:24]
  
  fit1 <- yy(npar1,nt,x1=HM[1],x2=DM[1])
  fit11<-phasic2(npar1[-c(1:6)], nt,x1=HM[1],x2=DM[1])
  fit.ind1 <- com.get_mu1(npar1[-c(1:6)][c(1,2,4,5)],nt,x1=HM[1],x2=DM[1])
  fit.inter1 <- fit1 - fit.ind1
  fit21.h<-phasic1(npar1[1:3],nt)
  fit21.d<-phasic1(npar1[4:6],nt)
  
  fit2 <- yy(npar2,nt,x1=HM[2],x2=DM[2])
  fit12<-phasic2(npar2[-c(1:6)], nt,x1=HM[2],x2=DM[2])
  fit.ind2 <- com.get_mu1(npar2[-c(1:6)][c(1,2,4,5)],nt,x1=HM[2],x2=DM[2])
  fit.inter2 <- fit2 - fit.ind2
  fit22.h<-phasic1(npar2[1:3],nt)
  fit22.d<-phasic1(npar2[4:6],nt)
  
  effect <- apply(cbind(fit1,fit2),1,sd)
  effect1 <- apply(cbind(fit11,fit12),1,sd)
  effect.ind <- apply(cbind(fit.ind1,fit.ind2),1,sd)
  effect.inter <- apply(cbind(fit.inter1,fit.inter2),1,sd)
  
  tiff("./Figures/figure6.tiff",width=20,height=7,units = "cm",res=300)
  par(mfrow=c(1,3))
  par(mar=c(4,4,2,1))
  plot(c(0,0), c(2,0),xlim=c(0,25), ylim=c(-0.2,1.7), xlab="Age (year)",
       ylab = "Height (m)",type="n", xaxt="n", yaxt="n",
       xaxs="i", yaxs="i",cex.lab=1.3,mgp=c(2.5,1,0))
  axis(1,seq(1,24,3),seq(1,24,3),cex.axis=1)
  axis(2,seq(-0.5,2.5,0.5),seq(-0.5,2.5,0.5),las = 1,cex.axis=1)
  lines(nt,effect[1:length(nt)],lwd=2,col="steelblue")
  lines(nt,effect1[1:length(nt)],col="tomato",lwd=2)
  lines(nt,effect.ind[1:length(nt)],col="tomato",lwd=2,lty=2)
  lines(nt,effect.inter[1:length(nt)],col="tomato",lwd=2,lty=3)
  lines(nt,apply(cbind(fit22.h,fit21.h),1,sd),col="orchid3",lwd=2)
  mtext("A",side=3, cex=1.2,adj=-0.17)

  par(mar=c(4,4,2,1))
  plot(c(0,0), c(2,0),xlim=c(0,25), ylim=c(-0.5,4.7), xlab="Age (year)",
       ylab = "Diameter (cm)",type="n", xaxt="n", yaxt="n",
       xaxs="i", yaxs="i",cex.lab=1.3,mgp=c(2.5,1,0))
  axis(1,seq(1,24,3),seq(1,24,3),cex.axis=1)
  axis(2,seq(0,5,1),seq(0,5,1),las = 1,cex.axis=1)
  lines(nt,effect[(1+length(nt)):(2*length(nt))],lwd=2,col="steelblue")
  lines(nt,effect1[(1+length(nt)):(2*length(nt))],col="tomato",lwd=2)
  lines(nt,effect.ind[(1+length(nt)):(2*length(nt))],col="tomato",lwd=2,lty=2)
  lines(nt,effect.inter[(1+length(nt)):(2*length(nt))],col="tomato",lwd=2,lty=3)
  lines(nt,apply(cbind(fit22.d,fit21.d),1,sd),col="orchid3",lwd=2)
  mtext("B",side=3, cex=1.2,adj=-0.17)

  par(mar=c(4,4,2,1))
  plot(c(0,0), c(2,0),xlim=c(-0.2,1.7), ylim=c(-0.3,4.2), xlab="Height (m)",
       ylab = "Diameter (cm)",type="n", xaxt="n", yaxt="n",
       xaxs="i", yaxs="i",cex.lab=1.3,mgp=c(2.5,1,0))
  axis(1,seq(-0.5,2.5,0.5),seq(-0.5,2.5,0.5),cex.axis=1)
  axis(2,seq(0,5,1),seq(0,5,1),las = 1,cex.axis=1)
  lines(effect[1:length(nt)],effect[(1+length(nt)):(2*length(nt))],lwd=2,col="steelblue")
  lines(effect1[1:length(nt)],effect1[(1+length(nt)):(2*length(nt))],col="tomato",lwd=2)
  lines(effect.ind[1:length(nt)],effect.ind[(1+length(nt)):(2*length(nt))],col="tomato",lwd=2,lty=2)
  lines(effect.inter[1:length(nt)],effect.inter[(1+length(nt)):(2*length(nt))],col="tomato",lwd=2,lty=3)
  lines(apply(cbind(fit22.h,fit21.h),1,sd),apply(cbind(fit22.d,fit21.d),1,sd),col="orchid3",lwd=2)
  mtext("C",side=3, cex=1.2,adj=-0.17)
  dev.off()
  
}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...,base_family="serif") +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(c(0.5,0.5,0.5,1.0), 'lines'),
          panel.background=element_rect(fill='transparent', color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),
          plot.title=element_text(size=15,vjust=2),
          axis.title = element_text(color='black',size=12),
          axis.title.x = element_text(vjust=0),
          axis.title.y = element_text(vjust=1.5),
          #axis.ticks.length = unit(0.1,"lines"),
          axis.text.x=element_text(vjust=0.5,size=8),
          axis.text.y=element_text(hjust=0.5,size=8),
          axis.ticks = element_line(color='black',size=0.3),
          #axis.ticks.margin = unit(0.8,"lines"),
          axis.line.x=element_line(linetype=1,color="black",size=0.2),
          axis.line.y=element_line(linetype=1,color="black",size=0.3),
          legend.position='none',
          legend.title=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
          strip.background=element_rect(fill='transparent', color='transparent'),
          strip.text=element_text(color='black',vjust=-3,hjust=0.06,size=16))
  
}


Figure7<-function(auc005,auc01){
  
  tiff("./Figures/Figure7.tiff",width=20,height=10,units = "cm",res=300)
  par(mfrow=c(1,2))
  par(mar=c(4,4.4,2,1),las=1)
  plot(c(0,auc005$TF.66[,2],1),c(0,auc005$TF.66[,1],1),type = "l",
       xlim=c(0,1),ylim=c(0,1),xlab="FPR",ylab="TPR",cex.lab=1.1,
       main = "Heritability=0.05",cex.main=1.2)
  lines(c(0,auc005$TF.100[,2],1),c(0,auc005$TF.100[,1],1),type = "l")
  lines(c(0,auc005$TF.200[,2],1),c(0,auc005$TF.200[,1],1),type = "l")
  polygon(c(0,1,rev(c(auc005$TF.200[,2],1))),c(0,0,rev(c(auc005$TF.200[,1],1))),
          col = rgb(254, 67, 101, 255, maxColorValue=255),lwd = 0.01, border = NA)
  polygon(c(0,1,rev(c(auc005$TF.100[,2],1))),c(0,0,rev(c(auc005$TF.100[,1],1))),
          col = rgb(252, 157, 154, 255, maxColorValue=255),lwd = 0.01, border = NA)
  polygon(c(0,1,rev(c(auc005$TF.66[,2],1))),c(0,0,rev(c(auc005$TF.66[,1],1))),
          col = rgb(249, 205, 173, 255, maxColorValue=255),lwd = 0.01, border = NA)
  lines(c(0,auc005$TF.100[,2],1),c(0,auc005$TF.100[,1],1),type = "l",lwd=2)
  lines(c(0,auc005$TF.200[,2],1),c(0,auc005$TF.200[,1],1),type = "l",lwd=2)
  lines(c(0,auc005$TF.66[,2],1),c(0,auc005$TF.66[,1],1),type = "l",lwd=2)
  lines(c(0,1),c(0,1),type="l",lwd=1,lty=2)
  text(0.25,0.93,"AUC1=0.3732",cex = 1)
  text(0.25,0.86,"AUC2=0.7076",cex = 1)
  text(0.25,0.79,"AUC3=0.9935",cex = 1)
  
  
  par(mar=c(4,4.4,2,1),las=1)
  plot(c(0,auc01$TF.66[,2],1),c(0,auc01$TF.66[,1],1),type = "l",
       xlim=c(0,1),ylim=c(0,1),xlab="FPR",ylab="TPR",cex.lab=1.1,
       main = "Heritability=0.10",cex.main=1.2)
  lines(c(0,auc01$TF.100[,2],1),c(0,auc01$TF.100[,1],1),type = "l")
  lines(c(0,auc01$TF.200[,2],1),c(0,auc01$TF.200[,1],1),type = "l")
  polygon(c(0,1,rev(c(auc01$TF.200[,2],1))),c(0,0,rev(c(auc01$TF.200[,1],1))),
          col = rgb(254, 67, 101, 255, maxColorValue=255),lwd = 0.01, border = NA)
  polygon(c(0,1,rev(c(auc01$TF.100[,2],1))),c(0,0,rev(c(auc01$TF.100[,1],1))),
          col = rgb(252, 157, 154, 255, maxColorValue=255),lwd = 0.01, border = NA)
  polygon(c(0,1,rev(c(auc01$TF.66[,2],1))),c(0,0,rev(c(auc01$TF.66[,1],1))),
          col = rgb(249, 205, 173, 255, maxColorValue=255),lwd = 0.01, border = NA)
  lines(c(0,auc01$TF.100[,2],1),c(0,auc01$TF.100[,1],1),type = "l",lwd=2)
  lines(c(0,auc01$TF.200[,2],1),c(0,auc01$TF.200[,1],1),type = "l",lwd=2)
  lines(c(0,auc01$TF.66[,2],1),c(0,auc01$TF.66[,1],1),type = "l",lwd=2)
  lines(c(0,1),c(0,1),type="l",lwd=1,lty=2)
  text(0.25,0.83,"AUC1=0.9659",cex = 1)
  text(0.25,0.76,"AUC2=0.9878",cex = 1)
  text(0.25,0.69,"AUC3=1.0000",cex = 1)
  
  dev.off()
}