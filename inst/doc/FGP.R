## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  install.packages("devtools")
#  library(devtools)
#  install_github("Shufeyangyi2015310117/LPG")

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
library(LPG)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
data(simulation)  

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  opts = list(max_iter = 1e5, dispF = 1, display_gap = 20, epsStopLogLik = 1e-5)
#  fit<-Lpg(X, y, opts = opts)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  opts = list(max_iter = 1e5, dispF = 1, display_gap = 20, epsStopLogLik = 1e-5)
#  fit<-Lpg(X, y, family = "binomial", opts = opts)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  opts = list(max_iter = 1e5, dispF = 1, display_gap = 20, epsStopLogLik = 1e-5, constraintalpha = 0)
#  fit<-Lpg(X, y, x2 = X2, y2 = y2, opts = opts)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  opts = list(max_iter = 1e5, dispF = 1, display_gap = 20, epsStopLogLik = 1e-5, constraintalpha = 0)
#  fit<-Lpg(X, y, x2 = X2, y2 = y2, family = "binomial", opts = opts)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  opts = list(max_iter = 1e5, dispF = 1, display_gap = 20, epsStopLogLik = 1e-5)
#  fit<-Lpg(X, y, z=z, family = "binomial", opts = opts)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  opts = list(max_iter = 1e5, dispF = 1, display_gap = 20, epsStopLogLik = 1e-5, constraintalpha = 0)
#  fit<-Lpg(X, y, z = z, x2 = X2, y2 = y2, z2 = z2, family = "binomial", opts = opts)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  opts = list(max_iter = 1e5, dispF = 1, display_gap = 20, epsStopLogLik = 1e-5, constraintalpha = 1)
#  fit<-Lpg(X, y, z = z, x2 = X2, y2 = y2, z2 = z2, family = "binomial", opts = opts)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
str(fit)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  opts = list(max_iter = 1e5, dispF = 1, display_gap = 20, epsStopLogLik = 1e-5)
#  fit <- Lpg.Plink("D:/realdata/WTCCC_all/BDqc37", opts = opts)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  opts = list(max_iter = 1e5, dispF = 1, display_gap = 20, epsStopLogLik = 1e-5)
#  fit <- Lpg.Plink("D:/realdata/WTCCC_all/BDqc37", family = "binomial", opts = opts)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  opts = list(max_iter = 1e5, dispF = 1, display_gap = 20, epsStopLogLik = 1e-5, constraintalpha = 0)
#  fit <- Lpg.Plink("D:/realdata/WTCCC_all/BDqc37",file2 = "D:/realdata/WTCCC_all/BDqc36", opts = opts)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  opts = list(max_iter = 1e5, dispF = 1, display_gap = 20, epsStopLogLik = 1e-5, constraintalpha = 0)
#  fit <- Lpg.Plink("D:/realdata/WTCCC_all/BDqc37", file2 =
#                     "D:/realdata/WTCCC_all/BDqc36", family = "binomial", opts = opts)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  opts = list(max_iter = 1e5, dispF = 1, display_gap = 20, epsStopLogLik = 1e-5)
#  fit <- Lpg.Plink("D:/realdata/WTCCC_all/BDqc37", z, family = "binomial", opts = opts)

## ---- eval=FALSE, echo=TRUE, tidy = FALSE--------------------------------
#  opts = list(max_iter = 1e5, dispF = 1, display_gap = 20, epsStopLogLik = 1e-5, constraintalpha = 0)
#  fit <- Lpg.Plink("D:/realdata/WTCCC_all/BDqc37", z = z, file2 =
#                     "D:/realdata/WTCCC_all/BDqc36", z2 = z2, family = "binomial", opts = opts)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  data <- Read.Plink("D:/realdata/WTCCC_all/BDqc37")

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
str(data) 

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  opts = list(max_iter = 1e5, dispF = 1, display_gap = 20, epsStopLogLik = 1e-5, constraintalpha = 0)
#  fit<-Lpg(X, y, x2 = X2, y2 = y2, family = "binomial", opts = opts)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
riskmat<-assoc(fit, FDR = 0.2, fdrControl="global")

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
str(riskmat)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  opts = list(max_iter = 1e5, dispF = 1, display_gap = 20, epsStopLogLik = 1e-5, constraintalpha = 1)
#  fit0<-Lpg(X, y, x2 = X2, y2 = y2, family = "binomial", opts = opts)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
out<-Pleiotropy.test(fit0,fit)

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
str(out)

