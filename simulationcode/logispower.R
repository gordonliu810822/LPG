rm(list=ls())
library("MASS")
library("Matrix")
library("glmnet")
library("broman")
library("LPG")
library("ROCR")

setwd("/home/gmsv1024/gmsv1024/logissimu")
source("/home/gmsv1024/gmsv1024/functionset.R")
source("/home/gmsv1024/gmsv1024/liability.R")


nrep = 100;          ## NO. of iterations
gamma = 0.3;         ## parameter to control pleiotropy
pi1 = 0.02;          ## proportion of nonzero random effects
rho=0.2;             ## parameter to control LD
L=100;               ## parameter to control LD
p = 20000;           ## NO. of SNPs
nGWAS=2;             ## NO. of GWAS
K=0.1;               ## disease prevalence
h=0.5;               ## heritability
n = c(3500, 3500)    ## NO. of samples include prediction samples


pi0 = 1-pi1;
alpha0 = c(pi0*(pi0+gamma*pi1),pi1*pi0*(1-gamma),pi1*pi0*(1-gamma),pi1*(pi1+gamma*pi0));
id_tra = list();
id_pre = list();
cutoff = c(3000, 3000);
id_tra[[1]] = c(1:1500,2001:3500);
id_tra[[2]] = c(1:1500,2001:3500);
id_pre[[1]] = 1501:2000;
id_pre[[2]] = 1501:2000;

pvalue = matrix(-1,nrep,1)
irep=1;
#set.seed(1);   ##work for all the following  random funciton
maf = runif(p,0.05,0.5);
MU=matrix(0,1,L);
SIGMA = matrix(0,L,L);
for (i in 1:L){
  for (j in 1:L){
    SIGMA[i,j] = rho^(abs(i-j));
  }
}
R = SIGMA;

## computation
for (irep in 1:nrep)
{ 
  print(irep);
  ## coefficients
  gamma0 = genegamma0(p,pi1,gamma);
  X00 = geneX(gamma0,p,pi1,h,K,MU,R,n,maf,L)
  X01 = X00[[2]];
  X02 = X00[[3]];
  rm(list = "X00");
  
  X01 = scale(X01)/sqrt(p*pi1);
  X02 = scale(X02)/sqrt(p*pi1);
  
  y1 = rep(c(1,0),each = n[1]/2);
  y2 = rep(c(1,0),each = n[2]/2);
  
  ly1 = y1;
  ly2 = y2;
  y1[y1==0]=-1;
  y2[y2==0]=-1;
  
  X1_tra = X01[id_tra[[1]], ];
  X2_tra = X02[id_tra[[2]], ];
  
  X1_pre = X01[id_pre[[1]],];
  X2_pre = X02[id_pre[[2]],];
  
  rm(list = c("X01","X02"))
  
  Y1_tra = y1[id_tra[[1]]];
  Y2_tra = y2[id_tra[[2]]]; 
  
  lY1_tra = ly1[id_tra[[1]]];
  lY2_tra = ly2[id_tra[[2]]];
  
  ## variational bayesian method
  print("variational bayesian joint analysis")
  obj = Lpg(X1_tra,Y1_tra, x2=X2_tra,y2=Y2_tra,family = "binomial");
  
  opts1 = list(max_iter = 1000,dispF = 1,display_gap = 10,epsStopLogLik = 1e-5, constraintalpha = 1)
  obj1 = Lpg(X1_tra,Y1_tra, x2=X2_tra,y2=Y2_tra,family = "binomial", opt = opts1);
  
  down = max(obj$Lq);
  up = max(obj1$Lq);
  
  stat = -2*(up-down)
  pvalue[irep] = pchisq(stat, 1, lower.tail = F);
  
  save(pvalue, file = "logisrho2power01.RData");
  
}


