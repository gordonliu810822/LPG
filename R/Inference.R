assoc<-function(object, FDR=0.2, fdrControl="global") {
    
    # check arguments
    
    if ( fdrControl != "global" & fdrControl != "local" ) {
      stop( "Invalid value for 'fdrControl' argument! It should be either 'global' or 'local'." )
    }
    
    if ( FDR < 0 | FDR > 1 ) {
      stop( "Invalid value for 'FDR' argument! It should be between zero and one." )
    }
    
    # association mapping
    
    fdrmat <- 1 - object$vardist_gamma
    

      # based on marginal FDR
      
      amat <- matrix( 0, nrow(fdrmat), ncol(fdrmat) )
      
      if ( fdrControl == "local" ) {
        # local FDR control
        
        message( "Info: Association mapping based on the local FDR control at level ", FDR, "." )
        
        amat[ fdrmat <= FDR ] <- 1
      } else if ( fdrControl == "global" ) {
        # global FDR control
        
        message( "Info: Association mapping based on the global FDR control at level ", FDR, "." )
        
        # direct approach for FDR control
        
        for ( j_pheno in 1:ncol(amat) ) {
          pp <- fdrmat[,j_pheno]
          pp.ordered <- sort(pp)
          pp.cum <- cumsum( pp.ordered ) / c(1:length(pp))
          cutoff <- max( pp.ordered[ pp.cum <= FDR ] )
          amat[ pp <= cutoff, j_pheno ] <- 1
        }  			
      }

    return(amat)
  }



##### pleiotropy test #####
Pleiotropy.test<-function(H0_fit, H1_fit){
  LRT = -2*(max(H0_fit$Lq) - max(H1_fit$Lq))
  p.value = pchisq(LRT,1,lower.tail = F)
  out = list(LRT = LRT, pvalue = p.value)
  return(out)
}
