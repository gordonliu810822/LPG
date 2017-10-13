

liability<-function(n,nsnp,m,h,K,geno)
{
# n = 3500;
# nsnp = 20000;
# m = 30;
# h = 0.5;
# K = 0.1;
# geno = X00[[1]][,loc1]
geno_ = geno;
for (j in 1:m)
{
geno_[, j] = (geno[, j] - mean(geno[, j]))/sd(geno[, j])/sqrt(m);
}
beta = rnorm(m);
L0 = geno_%*%beta;
L = L0 + rnorm(n/K+100)*sd(L0)*sqrt((1-h)/h);
he = var(L0)/var(L)
T = quantile(L, 1-K);
case_idx = sample(which(L > T), n/2);
ctrl_idx = sample(which(L < T), n/2);
idx = list();
idx$case = case_idx;
idx$ctrl = ctrl_idx;
idx$he = he
idx
}


gediscreteData = function(X,p,n,maf)
{
  L = 1;M = p;nGWAS=1;
  X00 = list(); Y  = list();
  X00[[1]] = X;
  for (k in 1:nGWAS){
    for (l in 1:L){
      index = (M*(l-1)+1): (M*l);
      AAprob = maf[index]^2;
      Aaprob = 2*maf[index]*(1-maf[index]);
      quanti = matrix(c(1-Aaprob-AAprob, 1- AAprob),p,2);  ## attention
      Xrefc = matrix(0,n,M);
      for (j in 1:M){
        cutoff = qnorm(quanti[j,]);
        Xrefc[X00[[k]][,j] < cutoff[1],j] = 0;
        Xrefc[X00[[k]][,j] >= cutoff[1] & X00[[k]][,j] < cutoff[2],j] = 1;  ## attention
        Xrefc[X00[[k]][,j] >= cutoff[2],j] = 2;
      }
      X00[[k]] = Xrefc;
    }
  }
  X00;
}


interpart = function(loc)
{
  nn = length(loc)
  out1 = loc
  for (i in 1:nn)
  {
    if (loc[i]%%100==0)  ## loc[i] is divisible by 100
    {
      out1[i] = loc[i]%/%100
    }else{
      out1[i] = (loc[i]+100)%/%100
    }
  }
  out1
}




geneX = function(gamma0,p,pi1,h,K,MU,R,n,maf,L)
{
  source("/home/gmsv1024/gmsv1024/functionset.R")
  source("/home/gmsv1024/gmsv1024/liability.R")
  ############################################
  ###   list non_block zero_block store block
  ###   matrix non_para zero_para
  ###   non_para$ Lnb Freq SFreq EFreq Smaf Emaf
  ###   zero_para$ Lzb Smaf Emaf
  ###   nnonblock denotes the number of nonblock
  ###   nzeroblock denotes the number of zeroblock  
  
  ###   loc is location in p
  ###   position is location in each block
loc1 = which(gamma0[,1]==1)
loc2 = which(gamma0[,2]==1)

position1 = (loc1+100)%%100
position1[position1==0]=100

position2 = (loc2+100)%%100
position2[position2==0]=100


non_para1 = as.data.frame(table(interpart(loc1)))
colnames(non_para1)[1] = "Lnb"
non_para1[,1] = as.numeric(as.character(non_para1[,1]))
non_para1$EFreq = cumsum(non_para1[,2])
non_para1$SFreq = non_para1$EFreq - non_para1$Freq + 1
nnonblock1 = dim(non_para1)[1]
non_para1$Smaf = (non_para1$Lnb-1)*100+1
non_para1$Emaf = (non_para1$Lnb)*100

zero_para1 = as.data.frame(setdiff(1:(p/L),non_para1$Lnb))
colnames(zero_para1)[1]="Lzb"
nzeroblock1 = dim(zero_para1)[1]
zero_para1$Smaf = (zero_para1$Lzb-1)*100+1
zero_para1$Emaf = (zero_para1$Lzb)*100

non_block1 = list();
zero_block1 = list();


non_para2 = as.data.frame(table(interpart(loc2)))
colnames(non_para2)[1] = "Lnb"
non_para2[,1] = as.numeric(as.character(non_para2[,1]))
non_para2$EFreq = cumsum(non_para2[,2])
non_para2$SFreq = non_para2$EFreq - non_para2$Freq + 1
nnonblock2 = dim(non_para2)[1]
non_para2$Smaf = (non_para2$Lnb-1)*100+1
non_para2$Emaf = (non_para2$Lnb)*100

zero_para2 = as.data.frame(setdiff(1:(p/L),non_para2$Lnb))
colnames(zero_para2)[1]="Lzb"
nzeroblock2 = dim(zero_para2)[1]
zero_para2$Smaf = (zero_para2$Lzb-1)*100+1
zero_para2$Emaf = (zero_para2$Lzb)*100

non_block2 = list();
zero_block2 = list(); 


for (kk in 1:nnonblock1)
{
  sX1=mvrnorm(n[1]/K+100,MU,R);
  sX00 = gediscreteData(sX1,L,n[1]/K+100,maf[non_para1$Smaf[kk]:non_para1$Emaf[kk]]);
  non_block1[[kk]] = sX00[[1]]
  if (kk==1)
  {
    gene = sX00[[1]][,position1[non_para1$SFreq[kk]:non_para1$EFreq[kk]]]
  }else{
    gene = cbind(gene,sX00[[1]][,position1[non_para1$SFreq[kk]:non_para1$EFreq[kk]]])
  }
}
idx1 = liability(n[1],p,p*pi1,h,K,gene)


for (kk in 1:(p/L))
{ 
  if (kk%in%non_para1$Lnb)
  {
    if (kk==1)
    {
      X01 = non_block1[[which(non_para1$Lnb==kk)]][c(idx1$case,idx1$ctrl),]
      non_block1[[which(non_para1$Lnb==kk)]]=0
    }else{
      X01 = cbind(X01,non_block1[[which(non_para1$Lnb==kk)]][c(idx1$case,idx1$ctrl),])
      non_block1[[which(non_para1$Lnb==kk)]]=0
    }
  }else{
    if (kk==1)
    {
      sX1=mvrnorm(n[1],MU,R);
      sX01 = gediscreteData(sX1,L,n[1],maf[zero_para1$Smaf[which(zero_para1$Lzb==kk)]:zero_para1$Emaf[which(zero_para1$Lzb==kk)]]);
      X01 = sX01[[1]]
    }else{
      sX1=mvrnorm(n[1],MU,R);
      sX01 = gediscreteData(sX1,L,n[1],maf[zero_para1$Smaf[which(zero_para1$Lzb==kk)]:zero_para1$Emaf[which(zero_para1$Lzb==kk)]])
      X01 = cbind(X01,sX01[[1]])
    }
  }
}
rm(list = c("sX1","sX01"))


for (kk in 1:nnonblock2)
{
  sX1=mvrnorm(n[2]/K+100,MU,R);
  sX00 = gediscreteData(sX1,L,n[2]/K+100,maf[non_para2$Smaf[kk]:non_para2$Emaf[kk]]);
  non_block2[[kk]] = sX00[[1]]
  if (kk==1)
  {
    gene = sX00[[1]][,position2[non_para2$SFreq[kk]:non_para2$EFreq[kk]]]
  }else{
    gene = cbind(gene,sX00[[1]][,position2[non_para2$SFreq[kk]:non_para2$EFreq[kk]]])
  }
}
idx2 = liability(n[2],p,p*pi1,h,K,gene)


for (kk in 1:(p/L))
{ 
  if (kk%in%non_para2$Lnb)
  {
    if (kk==1)
    {
      X02 = non_block2[[which(non_para2$Lnb==kk)]][c(idx2$case,idx2$ctrl),]
      non_block2[[which(non_para2$Lnb==kk)]]=0
    }else{
      X02 = cbind(X02,non_block2[[which(non_para2$Lnb==kk)]][c(idx2$case,idx2$ctrl),])
      non_block2[[which(non_para2$Lnb==kk)]]=0
    }
  }else{
    if (kk==1)
    {
      sX1=mvrnorm(n[2],MU,R);
      sX01 = gediscreteData(sX1,L,n[2],maf[zero_para2$Smaf[which(zero_para2$Lzb==kk)]:zero_para2$Emaf[which(zero_para2$Lzb==kk)]]);
      X02 = sX01[[1]]
    }else{
      sX1=mvrnorm(n[2],MU,R);
      sX01 = gediscreteData(sX1,L,n[2],maf[zero_para2$Smaf[which(zero_para2$Lzb==kk)]:zero_para2$Emaf[which(zero_para2$Lzb==kk)]])
      X02 = cbind(X02,sX01[[1]])
    }
  }
}
rm(list = c("sX1","sX01"))
################################################################################ 
out = list()
out$he = c(idx1$he,idx2$he)
out$X1 = X01
out$X2 = X02
out
}








