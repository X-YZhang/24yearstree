setwd("E:/G/林木24年新")

library(mvtnorm)
require("deSolve")
source("double.load.R")
source("ode2.gompertz.R")
source("double.covar.R")
source("double.curve.R")
source("double.optim.R")
source("double.plot.R")


dat <- com.DH.load(datt="./data/poplarMF65-201204-geno-v2.dat",
                   pheno="./data/phe-1.csv",
                   pheno1="./data/phe-4.csv",
                   snpinfo="./data/poplarMF65-201204-geno-stat.r2.csv",
                   nstart=1,nlen=156363, times=c(1:24))
load("dat.RData")

rets<-com.DH.est1(dat,interval=c(1,156362))
save(res,file = "res.RData")
load("res.RData")

load("E:/G/林木24年新/new-simulation/0.05/auc.RData")
auc005<-auc
load("E:/G/林木24年新/new-simulation/24yr-simulation/h01_data/auc.RData")
auc01<-auc
