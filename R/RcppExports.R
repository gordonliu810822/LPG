LogisWFMVS2GVB <- function(x,z,y, opts=NULL) {
  .Call('LPG_LogisWFMVS2GVB', PACKAGE = 'LPG', x,z,y, opts)
}

LogisWFMVS4GVB <- function(x,x2,z,z2,y,y2, opts=NULL) {
  .Call('LPG_LogisWFMVS4GVB', PACKAGE = 'LPG',x,x2,z,z2,y,y2,opts)
}

ILogisMVS4GVB <- function(stringname1,stringname2, opts=NULL) {
  .Call('LPG_ILogisMVS4GVB', PACKAGE = 'LPG', stringname1,stringname2, opts)
}

ILogisMVS2GVB <- function(stringname1, opts=NULL) {
  .Call('LPG_ILogisMVS2GVB', PACKAGE = 'LPG', stringname1, opts)
}

LogisMVS4GVB <- function(x,x2,y,y2,opts=NULL) {
  .Call('LPG_LogisMVS4GVB', PACKAGE = 'LPG', x,x2,y,y2,opts)
}

LogisMVS2GVB <- function(x,y, opts=NULL) {
  .Call('LPG_LogisMVS2GVB', PACKAGE = 'LPG', x,y, opts)
}

ReadPlinkRcpp <- function(stringname) {
  .Call('LPG_ReadPlinkRcpp', PACKAGE = 'LPG', stringname)
}

ILinearMVS2GVB <- function(stringname, opts=NULL) {
  .Call('LPG_ILinearMVS2GVB', PACKAGE = 'LPG', stringname, opts)
}

ILinearMVS4GVB <- function(stringname1,stringname2, opts=NULL) {
  .Call('LPG_ILinearMVS4GVB', PACKAGE = 'LPG', stringname1, stringname2, opts)
}

VBsl4 <- function(x1,x2,y1,y2,opts=NULL) {
  .Call('LPG_LinearMVS4GVB', PACKAGE = 'LPG', x1,x2,y1,y2,opts)
}

VBsl2 <- function(x,y, opts=NULL) {
  .Call('LPG_LinearMVS2GVB', PACKAGE = 'LPG', x,y, opts)
}

