## function Read.Plink
Read.Plink <- function(stringname) {
  .Call('LPG_ReadPlinkRcpp', PACKAGE = 'LPG', stringname)
}
