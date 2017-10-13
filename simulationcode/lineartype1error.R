rm(list=ls())
library("MASS")
library("Matrix")
library("glmnet")
library("broman")
library("LPG")
library("ROCR")

setwd("/home/gmsv1024/gmsv1024/linsimu")
source("/home/gmsv1024/gmsv1024/functionset.R")


nrep = 500;                       ## NO. of iterations
gamma = 0;                        ## parameter to control pleiotropy
pi1 = 0.005;                      ## proportion of nonzero random effects
rho=0.2;                          ## parameter to control LD
L=100;                            ## parameter to control LD
p = 20000;                        ## NO. of SNPs
nGWAS=2;                          ## NO. of GWAS
h=0.5;                            ## heritability
n = c(3500, 3500)                 ## NO. of samples include prediction samples


pi0 = 1-pi1;
alpha0 = c(pi0*(pi0+gamma*pi1),pi1*pi0*(1-gamma),pi1*pi0*(1-gamma),pi1*(pi1+gamma*pi0));
id_tra = list();
id_pre = list();
cutoff = c(3000, 3000);
id_tra[[1]] = 1:cutoff[1];
id_tra[[2]] = 1:cutoff[2];
id_pre[[1]] = (cutoff[1] + 1):n[1];
id_pre[[2]] = (cutoff[2] + 1):n[2];

pvalue = matrix(-1,nrep,1)
irep=1;
#set.seed(13);   ##work for all the following  random funciton
maf = runif(p,0.05,0.5);
sp=500;
MU=matrix(0,1,sp);
R = geneSIGMA(L,rho,sp);

## computation
for (irep in 1:nrep)
{ 
  print(irep)
  X01=mvrnorm(n[1],MU,R);
  X02=mvrnorm(n[2],MU,R);
  
  for (k in 1:(p/sp-1))
  {
    sX1=mvrnorm(n[1],MU,R);
    sX2=mvrnorm(n[2],MU,R); 
    X01 = cbind(X01,sX1);
    X02 = cbind(X02,sX2);
  }
  
  X00 = generatediscreteData(X01,X02,p,n[1],maf);
  X1 = scale(X00[[1]])/sqrt(pi1*p);
  X2 = scale(X00[[2]])/sqrt(pi1*p);
  
  ## coefficients
  gamma0 = genegamma0(p,pi1,gamma);
  
  # coefficients
  b = matrix(0,p,2);
  Beta = matrix(0,sum(gamma0[,1]==1),2);
  Beta[,1] =rnorm(sum(gamma0[,1]==1));
  Beta[,2] =rnorm(sum(gamma0[,2]==1));
  
  b[gamma0[,1]==1,1] = Beta[,1];  
  b[gamma0[,2]==1,2] = Beta[,2];
  NULLset1 = (gamma0[,1]==0);    ## attention this technique
  NULLset2 = (gamma0[,2]==0);    
  
  Y1_0 = X1%*%b[,1];
  Err1 = rnorm(n[1]);
  Y1 = X1%*%b[,1] + Err1*sd(Y1_0)*sqrt((1-h)/h);
  
  Y2_0 = X2%*%b[,2];
  Err2 = rnorm(n[2]);
  Y2 = X2%*%b[,2] + Err2*sd(Y2_0)*sqrt((1-h)/h);
  
  X1_tra = X1[id_tra[[1]], ];
  X2_tra = X2[id_tra[[2]], ];

  Y1_tra = Y1[id_tra[[1]],];
  Y2_tra = Y2[id_tra[[2]],]; 
  
  
  ## variational bayesian method
  print("variational bayesian joint analysis")
  obj = Lpg(X1_tra,Y1_tra,x2 = X2_tra,y2 = Y2_tra);
  
  opts1 = list(max_iter = 1000,dispF = 1,display_gap = 10,epsStopLogLik = 1e-5, constraintalpha = 1)
  obj1 = Lpg(X1_tra,Y1_tra, x2 = X2_tra,y2 = Y2_tra, opt = opts1);
  
  down = max(obj$Lq);
  up = max(obj1$Lq);
  
  stat = -2*(up-down)
  pvalue[irep] = pchisq(stat, 1, lower.tail = F);
  
  save(pvalue, file = "linrho2type1error_1.RData");

}

