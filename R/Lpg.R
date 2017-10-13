## fucntion Lpg
Lpg <-function(x,y,z=NULL,x2=NULL,y2=NULL,z2=NULL,family="normal", opts=NULL)
{
  ####################linear and logistic two group #############################
  if (is.null(z)&is.null(x2)&is.null(y2)&is.null(z2))
  {
    if (dim(x)[1]!=dim(as.matrix(y))[1]){
      stop("The rows of x is not equal to the rows of y")
    }

    if (dim(as.matrix(y))[2]!=1) {
      stop("the column of y is not equal to one")
    }

    if (family!="normal" && family!="binomial"){
      stop("Inappropriate value for 'family' argument")
    }

    if (family=="normal"){
      print("Separate analysis for linear model")
      .Call('LPG_LinearMVS2GVB', PACKAGE = 'LPG', x,y, opts)
    }else{
      print("Separate analysis for logistic model")
      .Call('LPG_LogisMVS2GVB', PACKAGE = 'LPG', x,y, opts)
    }

    ####################linear and logistic four group #########################
  }else if (is.null(z)&&is.null(z2)){
    if (dim(x)[1]!=dim(as.matrix(y))[1]){
      stop("The rows of x is not equal to the rows of y")
    }

    if (dim(as.matrix(as.matrix(y)))[2]!=1) {
      stop("the column of y is not equal to one")
    }

    if (dim(x2)[1]!=dim(as.matrix(y2))[1]){
      stop("The rows of x2 is not equal to the rows of y2")
    }

    if (dim(as.matrix(y2))[2]!=1) {
      stop("the column of y2 is not equal to one")
    }

    if (dim(x)[2]!=dim(x2)[2]){
      stop("The columns of x is not equal to The columns of x2")
    }

    if (family!="normal" && family!="binomial"){
      stop("Inappropriate value for 'family' argument")
    }

    if (family=="normal"){
      print("Joint analysis for linear model")
      .Call('LPG_LinearMVS4GVB', PACKAGE = 'LPG', x,x2,y,y2, opts)
    }else{
      print("Joint analysis for logistic model")
      .Call('LPG_LogisMVS4GVB', PACKAGE = 'LPG', x,x2,y,y2, opts)
    }

    ####################logistic two group with fixed term #########################
  }else if (is.null(x2)&&is.null(y2)&&is.null(z2)){

    if (dim(x)[1]!=dim(as.matrix(y))[1]){
      stop("The rows of x is not equal to the rows of y")
    }

    if (dim(as.matrix(y))[2]!=1) {
      stop("the column of y is not equal to one")
    }

    if (dim(z)[1]!=dim(x)[1]){
      stop("The rows of x is not equal to the rows of z")
    }

    if (family!="normal" && family!="binomial"){
      stop("Inappropriate value for 'family' argument")
    }

    if (family=="normal"){
      stop("The function hasn't developed!!")
    }

    if (family=="binomial"){
      print("Separate analysis for logistic model with adjusted covariates")
      .Call('LPG_LogisWFMVS2GVB', PACKAGE = 'LPG', x,z,y, opts)
    }

    ####################logistic four group with fixed term #########################
  }else if ((!is.null(x))&&(!is.null(y))&&(!is.null(z))&&
            (!is.null(x2))&&(!is.null(y2))&&(!is.null(z2))){

    if (family!="normal" && family!="binomial"){
      stop("Inappropriate value for 'family' argument")
    }

    if (family=="normal"){
      stop("The function hasn't developed!!")
    }

    if (family=="binomial")
    {
      print("Joint analysis for logistic model with adjusted covariates")
      .Call('LPG_LogisWFMVS4GVB', PACKAGE = 'LPG',x,x2,z,z2,y,y2, opts)
    }

  }else{
    stop("Inappropriate argument")
  }

}
