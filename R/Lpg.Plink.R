## function Lpg.Plink
Lpg.Plink<-function(file,z=NULL,file2=NULL,z2=NULL,family="normal", opts=NULL)
{
  if (is.null(z)&is.null(file2)&is.null(z2))
  {
    if (family!="normal" && family!="binomial"){
      stop("Inappropriate value for 'family' argument")
    }

    if (family=="normal"){
      print("Separate analysis for linear model")
      .Call('LPG_ILinearMVS2GVB', PACKAGE = 'LPG', file, opts)
    }else {
      print("Separate analysis for logistic model")
      .Call('LPG_ILogisMVS2GVB', PACKAGE = 'LPG', file, opts)
    }

  }else if (is.null(z)&is.null(z2)){
    if (family!="normal" && family!="binomial"){
      stop("Inappropriate value for 'family' argument")
    }

    if (family=="normal"){
      print("Joint analysis for linear model")
      .Call('LPG_ILinearMVS4GVB', PACKAGE = 'LPG', file, file2, opts)
    }else{
      print("Joint analysis for logistic model")
      .Call('LPG_ILogisMVS4GVB', PACKAGE = 'LPG', file, file2, opts)
    }
  }else if (is.null(file2)&is.null(z2)){
    if (family!="normal" && family!="binomial"){
      stop("Inappropriate value for 'family' argument")
    }

    if (family=="normal"){
      stop("The function hasn't developed!!")
    }

    if (family=="binomial"){
      print("Separate analysis for logistic model with adjusted covariates")
      .Call('LPG_ILogisWFMVS2GVB', PACKAGE = 'LPG', file,z, opts)
    }

  }else if ((!is.null(file))&&(!is.null(z))&&(!is.null(file2))&&(!is.null(z2))){
    if (family!="normal" && family!="binomial"){
      stop("Inappropriate value for 'family' argument")
    }

    if (family=="normal"){
      stop("The function hasn't developed!!")
    }

    if (family=="binomial")
    {
      print("Joint analysis for logistic model with adjusted covariants")
      .Call('LPG_ILogisWFMVS4GVB', PACKAGE = 'LPG', file,z,file2,z2, opts)
    }

  }else{
    stop("Inappropriate argument")
  }
}
