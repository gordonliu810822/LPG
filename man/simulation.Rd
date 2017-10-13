\name{simulation}
\alias{simulation}
\docType{data}
\title{
simulation data
}
\description{
simulation data
}
\usage{data("simulation")}
\format{
Two n-by-p genotpye matrix X, X2, two n-by-1 binary phenotype matrix y, y2, two n-by-q covariates z, z2.
}

\examples{
data(simulation)
fit<-Lpg(X1,y1, family = "binomial")
fit<-Lpg(X1,y1,x2 = X2,y2 = y2, family = "binomial")
}
\keyword{datasets}
