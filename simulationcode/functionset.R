genegamma0 = function(p,pi1,gamma){
gamma0 = matrix(0,p,2);    # true gamma for each run
pi0 = 1-pi1;
b1=sample(p,round(p*pi1));
gamma0[b1,1]=1;
matched = round(p*pi1*(pi1+ gamma*pi0));
nmatched= length(b1) - matched;
indx=1:p;
gamma0[sample(indx[gamma0[,1]==1],matched),2] = 1;
gamma0[sample(indx[gamma0[,1]==0],nmatched),2] = 1;
gamma0;
}



geneSIGMA = function(L,rho,p)
{
groupinfo = matrix(0,L,2);
blocksize = p/L;
for (l in 1:L)
{
    groupinfo[l,1] = (l-1)*blocksize + 1;
    groupinfo[l,2] = groupinfo[l,1]+ blocksize -1;
}

SIGMA = matrix(0,L,L);

for (i in 1:L){
    for (j in 1:L){
        SIGMA[i,j] = rho^(abs(i-j));
    }
}
R = SIGMA;
for (l in 1:(blocksize-1)){
    R = bdiag(R,SIGMA);
}
R;
}




generatediscreteData = function (X01,X02,p,n,maf)
{
L = 1;M = p;nGWAS = 2;
X00 = list(); Y  = list(); # Xo = cell(nGWAS,1);
# RRaw = cell(L,2); 
# stderror = zeros(nGWAS,1);

X00[[1]] = X01;
X00[[2]] = X02;
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



laperfcurve = function(gamma0,beta)
{
# compute lasso's ROC and AUC
nsnp = length(gamma0);
IX = sort(abs(beta),decreasing = TRUE,index.return = T);
nonzero = (which(gamma0==1));

TPR = matrix(0,nsnp,1);
FPR = matrix(0,nsnp,1);

for (k in 1:nsnp){
  
kTP = intersect(IX$ix[1:k],nonzero);

TPR[k] = length(kTP)/length(nonzero);

# significance  bisect some belong to nonzero, others belongs to zero
FPR[k] = (length(IX$ix[1:k]) - length(kTP))/(nsnp - length(nonzero));
}

# plot(FPR,TPR)
#AUC = trapz(FPR,TPR);###attention have no corresponding function
area = 0;
area = area + TPR[1]*FPR[1]/2;
for (i in 2:nsnp)
  {
    area = area + (TPR[i] + TPR[i-1])*(FPR[i] - FPR[i-1])/2;
  }
AUC = area;
result = list();
result$FPR = FPR;
result$TPR = TPR;
result$AUC = AUC;
result;
}



FDR = function(Zmarg, fdr0,nonzero) 
{
#Zmarg = fdr_1;fdr0=0.2;nonzero = indx[gamma0[,1]==1];
Sort = sort(Zmarg,index.return = T);
Zsort = Sort$x;
IX = Sort$ix;
nsnp = length(Zsort);
# fDR = cumsum(Zsort.*ones(nsnp,1))./(1:nsnp)';
fDR = cumsum(Zsort*matrix(1,nsnp,1))/(1:nsnp);## globe fdr the range of globe is limited
tmp = fDR - fdr0;
markerSel = IX[tmp<0];
tp = intersect(markerSel,nonzero);
power = length(tp) /length(nonzero); ## positive discovery rate
fdr = (length(markerSel) - length(tp))/length(markerSel);
if (length(markerSel) == 0)
fdr = 0;
result = list();
result$power = power;
result$fdr = fdr;
result;
## comment
# tp is true positive
# markersel means the selected genes, also is significant
# according to the formula in the paper
}

