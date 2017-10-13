\name{Read.Plink}
\alias{Read.Plink}
\title{
Read Plink data
}
\description{
Read Plink data
}
\usage{
Read.Plink(file)
}
\arguments{
  \item{file}{is a character of Plink file}
}

\value{ 
The outcome is a list contains a n-by-p genotype matrix "X" and n-by-1 phenotype matrix "Y"
}

\author{
Yang Yi Jin Liu
}

\examples{
# read GWAS data
data<-Read.Plink("D:/realdata/WTCCC_all/BDqc37")
}
