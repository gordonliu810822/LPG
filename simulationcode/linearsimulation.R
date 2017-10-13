rm(list=ls())
library("MASS")
library("Matrix")
library("glmnet")
library("broman")
library("LPG")
library("ROCR")

setwd("/home/gmsv1024/gmsv1024/linsimu")
source("/home/gmsv1024/gmsv1024/functionset.R")


nrep = 50;                        ## NO. of iterations
gamma = 1;                        ## parameter to control pleiotropy
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
lasso_power0 = matrix(0,nrep,nGWAS);
lasso_Fdrc0 = matrix(0,nrep,nGWAS);
Cor_rep4 = matrix(0,nrep,2);
Cor_rep2 = matrix(0,nrep,2);
Cor_rep_lasso = matrix(0,nrep,2);
irep=1;
#set.seed(1);   ##work for all the following  random funciton
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
  he_true_rep[irep,1] = var(Y1_0)/var(Y1);
  
  Y2_0 = X2%*%b[,2];
  Err2 = rnorm(n[2]);
  Y2 = X2%*%b[,2] + Err2*sd(Y2_0)*sqrt((1-h)/h);
  he_true_rep[irep,2] = var(Y2_0)/var(Y2);
  
  X1_tra = X1[id_tra[[1]], ];
  X2_tra = X2[id_tra[[2]], ];
  
  Y1_tra = Y1[id_tra[[1]],];
  Y2_tra = Y2[id_tra[[2]],]; 
  
  ## variational bayesian method
  print("variational bayesian joint analysis")
  obj = Lpg(X1_tra,Y1_tra,x2=X2_tra,y2=Y2_tra);
  print("variational bayesian separate analysis")
  obj1 = Lpg(X1_tra, Y1_tra);
  obj2 = Lpg(X2_tra, Y2_tra);
  
  ## lasso
  print("lasso analysis")
  Lasso_obj1 = cv.glmnet(X1_tra,Y1_tra);
  Lasso_obj2 = cv.glmnet(X2_tra,Y2_tra);
  
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
  phi1 = obj$vardist_gamma[,4]+obj$vardist_gamma[,3];
  pred <- prediction(phi1,gamma0[,1])  
  AUC = performance(pred,'auc')@y.values[[1]]; #AUCֵ  
  perf <- performance(pred,'tpr','fpr')  
  AUC_rep[irep,1] = AUC;
  
  phi1 = obj$vardist_gamma[,4]+obj$vardist_gamma[,2];
  pred <- prediction(phi1,gamma0[,2])  
  AUC = performance(pred,'auc')@y.values[[1]]; #AUCֵ  
  perf <- performance(pred,'tpr','fpr')  
  AUC_rep[irep,2] = AUC;
  
  # True Global False Discovery Rate
  fdr_1=1-obj$vardist_gamma[,4]-obj$vardist_gamma[,3];
  fdr_2=1-obj$vardist_gamma[,4]-obj$vardist_gamma[,2];
  
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
  X1_pre = X1[id_pre[[1]],];
  X2_pre = X2[id_pre[[2]],];
  Y1_pre = Y1[id_pre[[1]]];
  Y2_pre = Y2[id_pre[[2]]];
  
  c0 = c(mean(Y1_pre),mean(Y2_pre));
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
  Y1_hat =c0[1] + X00c[[1]]%*%(obj$vardist_mu[,1]*obj$vardist_gamma[,1]);
  Cor_rep4[irep,1] = cor(Y1_pre,Y1_hat);
  
  # study 2
  Y2_hat =c0[2] + X00c[[2]]%*%(obj$vardist_mu[,2]*obj$vardist_gamma[,2]);
  Cor_rep4[irep,2] = cor(Y2_pre,Y2_hat);
  
  # two groups
  # study 1
  Y1_hat = c0[1] + X00c[[1]]%*%(obj1$vardist_mu*obj1$vardist_gamma);
  Cor_rep2[irep,1] = cor(Y1_pre,Y1_hat);
  # study 2
  Y2_hat = c0[2] + X00c[[2]]%*%(obj2$vardist_mu*obj2$vardist_gamma);
  Cor_rep2[irep,2] = cor(Y2_pre,Y2_hat);
  
  # lasso
  # study 1
  beta = Lasso_obj1$glmnet.fit$beta[,Lasso_obj1$glmnet.fit$lambda==Lasso_obj1$lambda.min];
  Y1_hat = c0[1] + X00c[[1]]%*%beta + Lasso_obj1$glmnet.fit$a0[Lasso_obj1$glmnet.fit$lambda==Lasso_obj1$lambda.min];
  Cor_rep_lasso[irep,1] = cor(Y1_pre,Y1_hat);
  
  # study 2
  beta = Lasso_obj2$glmnet.fit$beta[,Lasso_obj2$glmnet.fit$lambda==Lasso_obj2$lambda.min];
  Y2_hat = c0[2] + X00c[[2]]%*%beta + Lasso_obj2$glmnet.fit$a0[Lasso_obj1$glmnet.fit$lambda==Lasso_obj1$lambda.min];
  Cor_rep_lasso[irep,2] = cor(Y2_pre,Y2_hat);
  
}

save(  AUC_rep , power0 , Fdrc0 ,
       AUC_rep1 , power01 , Fdrc01 ,
       AUC_rep2 , power02 , Fdrc02 ,
       lasso_AUC_rep, Cor_rep4 , Cor_rep2 , Cor_rep_lasso,
       file = "linrho2gamma04.RData");

