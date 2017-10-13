\name{assoc}
\alias{Lpg}
\title{
Fit Lpg model
}
\description{
Fit Lpg model
}
\usage{
Lpg <-function(x,y,z=NULL,x2=NULL,y2=NULL,z2=NULL,family="normal",opts = NULL)
}
\arguments{
    \item{x}{The first n-by-p genotype matrix from GWAS data,
		where row and column correspond to individual and SNP, respectively. }
    \item{y}{The first n-by-1 phenotype matrix from GWAS data,
		where row corresponds to phenotype, the column is only equal to one, 
        for binary scenario, the case is 1, the control is -1}
    \item{z}{The first n-by-q covariates matrix,
        where row and column correspond to individual and the number of covariates, respectively } 
    \item{x2}{The Second n-by-p genotype matrix from GWAS data,
		where row and column correspond to individual and SNP, respectively. }
    \item{y2}{The Second n-by-1 phenotype matrix from GWAS data,
		where row corresponds to phenotype, the column is only equal to one,
        for binary scenario, the case is 1, the control is -1}
    \item{z2}{The Second n-by-q covariates matrix,
        where row and column correspond to individual and the number of covariates, respectively }            
    \item{family}{is a character either "normal" or "binomial",
        where normal and binomial correspond to quantative phenotype and binary phenotype, respectively
        The default value is normal}
    \item{opts}{refer to Lpg.Plink}        
}

\value{refer to value in Lpg.Plink}

\author{
Yang Yi Jin Liu
}

\examples{
# Ftp for separate analysis with quantitative phenotype
opts = list(max_iter = 1e5, dispF = 1, display_gap = 20, epsStopLogLik = 1e-5)
fit<-Lpg(X, y, opts = opts)

# Lpg for separate analysis with binary phenotype
opts = list(max_iter = 1e5, dispF = 1, display_gap = 20, epsStopLogLik = 1e-5)
fit<-Lpg(X, y, family = "binomial", opts = opts)

# Lpg for joint analysis with quantitative phenotype
opts = list(max_iter = 1e5, dispF = 1, display_gap = 20, epsStopLogLik = 1e-5, epsStopLogLik = 0)
fit<-Lpg(X,y,x2 = X2,y2 = y2, opts = opts)

# Lpg for joint analysis with binary phenotype
opts = list(max_iter = 1e5, dispF = 1, display_gap = 20, epsStopLogLik = 1e-5, epsStopLogLik = 0)
fit<-Lpg(X,y,x2 = X2,y2 = y2,family = "binomial", opts = opts)

# Lpg for separate analysis with binary phenotype
opts = list(max_iter = 1e5, dispF = 1, display_gap = 20, epsStopLogLik = 1e-5)
fit<-Lpg(X,y,z=z,family = "binomial", opts = opts)

# Lpg for joint analysis with binary phenotype
opts = list(max_iter = 1e5, dispF = 1, display_gap = 20, epsStopLogLik = 1e-5, epsStopLogLik = 0)
fit<-Lpg(X,y,z = z,x2 = X2,y2 = y2,z2 = z2,family = "binomial", opts = opts)
}

