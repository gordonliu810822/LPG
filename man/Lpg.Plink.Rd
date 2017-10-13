\name{Lpg.Plink}
\alias{Lpg.Plink}
\title{
Fit Lpg model
}
\description{
Fit Lpg model
}
\usage{
Lpg.Plink(file,z=NULL,file2=NULL,z2=NULL,family="nomral",opts = NULL)
}
\arguments{
    \item{file}{is a character denotes the plink file name of first study
}  
    \item{z}{is a n-by-q matrix denotes the covariates of first study
}  
    \item{file2}{is a character denotes the plink file name of second study
}  
    \item{z2}{is a n-by-q matrix denotes the covariates of second study
}  
   \item{family}{is a character either "normal" or "binomial",
        where "normal" and "binomial" correspond to quantative phenotype and binary phenotype, respectively.
        The default value is "normal"} 
   \item{opts}{is a list consists of "max_iter", "dispF", "display_gap", "epsStopLogLik", "constraintalpha". 
        max_iter is an 
        positive integer which denotes the maximal iterations, default value 1e5. dispF is a logic value either 
        1 or 0, controlling whether or not 
        showing interations information, default value 1. display_gap is a positive integer which denotes the gap 
        showing interations information, default value 10. epsStopLogLik is a positive value which control 
        interation digits, default value 1e-5. constraintalpha is a logic value either 1 or 0, controlling whether
        alpha is constrained, default value 0. Noted that for separate analsis, constraintalpha does't exist in opts. 
        }
}

\value{
  The returned object is a list,
  \item{vardist_gamma}{For separate analysis, it is a p-by-1 matrix denoting the true posterior probability
  of a SNP associated with phenotype, For joint analysis, it is a p-by-2 matrix whose first and second column denote
  the true posterior probability of SNP in study 1 associated with phenotype 1 and SNP in study 2 associated 
  with phenotype 2, respectively.}
  \item{vardist_mu}{For separate analysis, it is a p-by-1 matrix denoting the posterior mean of random effects.
  For joint analysis, it is a p-by-2 matrix whose first and second column denots the posterior mean of random effects
  in first and second study, respectively.}
  \item{vardist_sigma2beta}{For separate analysis, it is a p-by-1 matrix denoting the posterior variance of random 
  effects. For joint analysis, it is a p-by-2 matrix whose first and second column denote the posterior variance
   of random effects in first and second study, respectively.}
  \item{sigma2beta}{is the variance of the random effects for corresponding study.}
  \item{sigma2e}{exists only for quatitative phenotype model, is the variance of error for corresponding study.}
  \item{alpha}{It denotes the parameter of distribution of indicator variable, for separate analysis, it is a scalar,
  for joint analysis, it is a vector.}
  \item{Lq}{is a row vector denoting the lower bound value of every iteration.}
  \item{u}{is the coefficients of the covriates, noted that the first row is the constant term.}
}
\author{
Yang Yi Jin Liu
}

\examples{
Separate analysis for quantitative phenotype
opts = list(max_iter = 1e5, dispF = 1, display_gap = 20, epsStopLogLik = 1e-5)
fit <- Lpg.Plink("D:/realdata/WTCCC_all/BDqc37", opts = opts)

Separate analysis for binary phenotype
opts = list(max_iter = 1e5, dispF = 1, display_gap = 20, epsStopLogLik = 1e-5)
fit <- Lpg.Plink("D:/realdata/WTCCC_all/BDqc37",family = "binomial", opts = opts)

Joint analysis for quantitative phenotype
opts = list(max_iter = 1e5, dispF = 1, display_gap = 20, epsStopLogLik = 1e-5, constraintalpha = 0)
fit <- Lpg.Plink("D:/realdata/WTCCC_all/BDqc37",file2 = "D:/realdata/WTCCC_all/BDqc36", opts = opts)

Joint analysis for binary phenotype
opts = list(max_iter = 1e5, dispF = 1, display_gap = 20, epsStopLogLik = 1e-5, constraintalpha = 0)
fit <- Lpg.Plink("D:/realdata/WTCCC_all/BDqc37", file2 = "D:/realdata/WTCCC_all/BDqc36", family = "binomial", opts = opts)

Separate analysis for binary phenotype with covariates
opts = list(max_iter = 1e5, dispF = 1, display_gap = 20, epsStopLogLik = 1e-5)
fit <- Lpg.Plink("D:/realdata/WTCCC_all/BDqc37",z,family = "binomial", opts = opts)

Joint analysis for binary phenotype with covariates
opts = list(max_iter = 1e5, dispF = 1, display_gap = 20, epsStopLogLik = 1e-5, constraintalpha = 0)
fit <- Lpg.Plink("D:/realdata/WTCCC_all/BDqc37", z = z, file2 = "D:/realdata/WTCCC_all/BDqc36", z2 = z2, family = "binomial", opts = opts)
}


