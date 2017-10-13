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


nrep = 50;           ## NO. of iterations
gamma = 0;           ## parameter to control pleiotropy
pi1 = 0.0025;        ## proportion of nonzero random effects
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

AUC_rep = matrix(0,nrep,nGWAS); 
Fdrc0 = matrix(0,nrep,nGWAS);
power0 = matrix(0,nrep,nGWAS);
AUC_rep1 = matrix(0,nrep,1);
power01  = matrix(0,nrep,1);
Fdrc01  = matrix(0,nrep,1);
AUC_rep2 = matrix(0,nrep,1);
power02  = matrix(0,nrep,1);
Fdrc02  = matrix(0,nrep,1);
lasso_AUC_rep= matrix(0,nrep,nGWAS);
Cor_rep4 = matrix(0,nrep,2);
Cor_rep2 = matrix(0,nrep,2);
Cor_rep_lasso = matrix(0,nrep,2);

irep=1;
set.seed(1);   ##work for all the following  random funciton
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
  obj = Lpg(X1_tra,Y1_tra,x2=X2_tra,y2=as.matrix(Y2_tra), family = "binomial")
  
  print("variational bayesian separate analysis")
  obj1 = Lpg(X1_tra, Y1_tra, family = "binomial");
  obj2 = Lpg(X2_tra, Y2_tra, family = "binomial");
  
  ## lasso
  print("lasso analysis")
  Lasso_obj1 = cv.glmnet(X1_tra,lY1_tra,family="binomial");
  Lasso_obj2 = cv.glmnet(X2_tra,lY2_tra,family="binomial");
  
  #computate AUC
  #study 1
  beta = Lasso_obj1$glmnet.fit$beta[,Lasso_obj1$glmnet.fit$lambda==Lasso_obj1$lambda.min];
  Perf = laperfcurve(gamma0[,1],beta);
  lasso_AUC_rep[irep,1] = Perf$AUC;
  
  #study 2
  beta = Lasso_obj2$glmnet.fit$beta[,Lasso_obj2$glmnet.fit$lambda==Lasso_obj2$lambda.min];
  Perf = laperfcurve(gamma0[,2],beta);   
  lasso_AUC_rep[irep,2] = Perf$AUC;

  ## VB four group
  # AUC
  # study 1
  phi1 = obj$vardist_gamma[,1];
  pred <- prediction(phi1,gamma0[,1])  
  AUC = performance(pred,'auc')@y.values[[1]]; #AUCֵ  
  perf <- performance(pred,'tpr','fpr')  
  AUC_rep[irep,1] = AUC;
  
  phi1 = obj$vardist_gamma[,2];
  pred <- prediction(phi1,gamma0[,2])  
  AUC = performance(pred,'auc')@y.values[[1]]; #AUCֵ  
  perf <- performance(pred,'tpr','fpr')  
  AUC_rep[irep,2] = AUC;
  
  # True Global False Discovery Rate
  fdr_1=1-obj$vardist_gamma[,1];
  fdr_2=1-obj$vardist_gamma[,2];
  
  indx = 1:p;
  out = FDR(fdr_1,0.2,indx[gamma0[,1]==1]);
  power0[irep,1] = out$power;
  Fdrc0[irep,1] = out$fdr;
  
  indx = 1:p;
  out = FDR(fdr_2,0.2,indx[gamma0[,2]==1]);
  power0[irep,2] = out$power;
  Fdrc0[irep,2] = out$fdr;
  
  ## VB two
  # study 1
  pred <- prediction(obj1$vardist_gamma,gamma0[,1])  
  AUC = performance(pred,'auc')@y.values[[1]]; #AUCֵ  
  perf <- performance(pred,'tpr','fpr')  
  AUC_rep1[irep] = AUC;
  
  indx = 1:p;
  out = FDR(1-obj1$vardist_gamma,0.2,indx[gamma0[,1]==1]);
  power01[irep,1] = out$power;
  Fdrc01[irep,1] = out$fdr;
  
  # study 2
  pred <- prediction(obj2$vardist_gamma,gamma0[,2])  
  AUC = performance(pred,'auc')@y.values[[1]]; #AUCֵ  
  perf <- performance(pred,'tpr','fpr')  
  AUC_rep2[irep] = AUC;
  
  indx = 1:p;
  out = FDR(1-obj2$vardist_gamma,0.2,indx[gamma0[,2]==1]);
  power02[irep,1] = out$power;
  Fdrc02[irep,1] = out$fdr;
  
  ## prediction
  Y1_pre = y1[id_pre[[1]]];
  Y2_pre = y2[id_pre[[2]]];
  
  c = list();
  c[[1]] = apply(X1_pre,2,mean);  ## attention
  c[[2]] = apply(X2_pre,2,mean);
  X00c = list();
  X00c[[1]] = matrix(0,500,p);
  X00c[[2]] = matrix(0,500,p);
  for (k in 1:p)
  {
    X00c[[1]][,k] = X1_pre[,k] - c[[1]][k];
    X00c[[2]][,k] = X2_pre[,k] - c[[2]][k];  
  }
  
  # four groups
  # study 1
  Y1_hat = obj$u[1,1] + X00c[[1]]%*%(obj$vardist_mu[,1]*obj$vardist_gamma[,1]);
  Y1_p = 1/(1 + exp(-Y1_hat));
  pred <- prediction(Y1_p,Y1_pre);
  AUC = performance(pred,'auc')@y.values[[1]]; #AUCֵ 
  Cor_rep4[irep,1] = AUC;
  
  # study 2
  Y2_hat = obj$u[1,2] + X00c[[2]]%*%(obj$vardist_mu[,2]*(obj$vardist_gamma[,2]));
  Y2_p = 1/(1 + exp(-Y2_hat));
  pred <- prediction(Y2_p,Y2_pre);
  AUC = performance(pred,'auc')@y.values[[1]]; #AUCֵ 
  Cor_rep4[irep,2] = AUC;
  
  # two groups
  # study 1
  Y1_hat = obj1$u + X00c[[1]]%*%(obj1$vardist_mu*obj1$vardist_gamma);
  Y1_p = 1/(1 + exp(-Y1_hat));
  pred <- prediction(Y1_p,Y1_pre);
  AUC = performance(pred,'auc')@y.values[[1]]; #AUCֵ 
  Cor_rep2[irep,1] = AUC;
  
  # study 2
  Y2_hat = obj2$u +  X00c[[2]]%*%(obj2$vardist_mu*obj2$vardist_gamma);
  Y2_p = 1/(1 + exp(-Y2_hat));
  pred <- prediction(Y2_p,Y2_pre);
  AUC = performance(pred,'auc')@y.values[[1]]; #AUCֵ
  Cor_rep2[irep,2] = AUC;
  
  # lasso
  # study 1
  beta = Lasso_obj1$glmnet.fit$beta[,Lasso_obj1$glmnet.fit$lambda==Lasso_obj1$lambda.min];
  Y1_hat = X00c[[1]]%*%beta + Lasso_obj1$glmnet.fit$a0[Lasso_obj1$glmnet.fit$lambda==Lasso_obj1$lambda.min];
  Y1_p = 1/(1 + exp(-Y1_hat));
  pred <- prediction(Y1_p,Y1_pre);
  AUC = performance(pred,'auc')@y.values[[1]]; #AUCֵ 
  Cor_rep_lasso[irep,1] = AUC;
  
  # study 2
  beta = Lasso_obj2$glmnet.fit$beta[,Lasso_obj2$glmnet.fit$lambda==Lasso_obj2$lambda.min];
  Y2_hat = X00c[[2]]%*%beta + Lasso_obj2$glmnet.fit$a0[Lasso_obj2$glmnet.fit$lambda==Lasso_obj2$lambda.min];
  Y2_p = 1/(1 + exp(-Y2_hat));
  pred <- prediction(Y2_p,Y2_pre);
  AUC = performance(pred,'auc')@y.values[[1]]; #AUCֵ
  Cor_rep_lasso[irep,2] = AUC;
  
}

save(  AUC_rep , power0 , Fdrc0 ,
       AUC_rep1 , power01 , Fdrc01 ,
       AUC_rep2 , power02 , Fdrc02 ,
       lasso_AUC_rep, Cor_rep4 , Cor_rep2, Cor_rep_lasso,
       file = "logisrho2gamma00.RData");
