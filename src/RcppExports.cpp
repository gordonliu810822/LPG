#include <RcppArmadillo.h>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <bitset>
#include <string>


using namespace Rcpp;
using namespace std;
using namespace arma;


/////////////////////////////////////////////////////////////////////////////////////
/////                                                                           /////
/////                            Linear model                                   /////
/////                                                                           /////
/////////////////////////////////////////////////////////////////////////////////////


// LinearMVS2GVB
Rcpp::List LinearMVS2GVB(arma::mat& Xr, arma::colvec& yr, SEXP opts);
RcppExport SEXP LPG_LinearMVS2GVB(SEXP XrSEXP, SEXP yrSEXP, SEXP optsSEXP) {
	BEGIN_RCPP
	Rcpp::RObject rcpp_result_gen;
	Rcpp::RNGScope rcpp_rngScope_gen;
	Rcpp::traits::input_parameter< arma::mat& >::type Xr(XrSEXP);
	Rcpp::traits::input_parameter< arma::colvec& >::type yr(yrSEXP);
	Rcpp::traits::input_parameter< SEXP >::type opts(optsSEXP);
	rcpp_result_gen = Rcpp::wrap(LinearMVS2GVB(Xr, yr, opts));
	return rcpp_result_gen;
	END_RCPP
}

// LinearMVS4GVB
Rcpp::List LinearMVS4GVB(arma::mat&  Xr, arma::mat&  Xr2, arma::colvec&  yr, arma::colvec&  yr2, SEXP opts);
RcppExport SEXP LPG_LinearMVS4GVB(SEXP XrSEXP, SEXP Xr2SEXP, SEXP yrSEXP, SEXP yr2SEXP, SEXP optsSEXP) {
	BEGIN_RCPP
	Rcpp::RObject rcpp_result_gen;
	Rcpp::RNGScope rcpp_rngScope_gen;
	Rcpp::traits::input_parameter< arma::mat& >::type Xr(XrSEXP);
	Rcpp::traits::input_parameter< arma::mat& >::type Xr2(Xr2SEXP);
	Rcpp::traits::input_parameter< arma::colvec& >::type yr(yrSEXP);
	Rcpp::traits::input_parameter< arma::colvec& >::type yr2(yr2SEXP);
	Rcpp::traits::input_parameter< SEXP >::type opts(optsSEXP);
	rcpp_result_gen = Rcpp::wrap(LinearMVS4GVB(Xr, Xr2, yr, yr2, opts));
	return rcpp_result_gen;
	END_RCPP
}

// ILinearMVS2GVB
Rcpp::List ILinearMVS2GVB(string stringname, SEXP opts);
RcppExport SEXP LPG_ILinearMVS2GVB(SEXP stringnameSEXP, SEXP optsSEXP) {
	BEGIN_RCPP
	Rcpp::RObject rcpp_result_gen;
	Rcpp::RNGScope rcpp_rngScope_gen;
	Rcpp::traits::input_parameter< string >::type stringname(stringnameSEXP);
	Rcpp::traits::input_parameter< SEXP >::type opts(optsSEXP);
	rcpp_result_gen = Rcpp::wrap(ILinearMVS2GVB(stringname, opts));
	return rcpp_result_gen;
	END_RCPP
}

// ILinearMVS4GVB
Rcpp::List ILinearMVS4GVB(string stringname1, string stringname2, SEXP opts);
RcppExport SEXP LPG_ILinearMVS4GVB(SEXP stringname1SEXP, SEXP stringname2SEXP, SEXP optsSEXP) {
	BEGIN_RCPP
	Rcpp::RObject rcpp_result_gen;
	Rcpp::RNGScope rcpp_rngScope_gen;
	Rcpp::traits::input_parameter< string >::type stringname1(stringname1SEXP);
	Rcpp::traits::input_parameter< string >::type stringname2(stringname2SEXP);
	Rcpp::traits::input_parameter< SEXP >::type opts(optsSEXP);
	rcpp_result_gen = Rcpp::wrap(ILinearMVS4GVB(stringname1, stringname2, opts));
	return rcpp_result_gen;
	END_RCPP
}


/////////////////////////////////////////////////////////////////////////////////////
/////                                                                           /////
/////                            Logistic model                                 /////
/////                                                                           /////
/////////////////////////////////////////////////////////////////////////////////////


// LogisMVS2GVB
Rcpp::List LogisMVS2GVB(arma::mat& Xr, arma::colvec& yr, SEXP opts);
RcppExport SEXP LPG_LogisMVS2GVB(SEXP XrSEXP, SEXP yrSEXP, SEXP optsSEXP) {
	BEGIN_RCPP
		Rcpp::RObject rcpp_result_gen;
	Rcpp::RNGScope rcpp_rngScope_gen;
	Rcpp::traits::input_parameter< arma::mat& >::type Xr(XrSEXP);
	Rcpp::traits::input_parameter< arma::colvec& >::type yr(yrSEXP);
	Rcpp::traits::input_parameter< SEXP >::type opts(optsSEXP);
	rcpp_result_gen = Rcpp::wrap(LogisMVS2GVB(Xr, yr, opts));
	return rcpp_result_gen;
	END_RCPP
}

// LogisMVS4GVB
Rcpp::List LogisMVS4GVB(arma::mat & Xr, arma::mat & Xr2, arma::colvec & yr, arma::colvec & yr2, SEXP opts);
RcppExport SEXP LPG_LogisMVS4GVB(SEXP XrSEXP, SEXP Xr2SEXP, SEXP yrSEXP, SEXP yr2SEXP, SEXP optsSEXP) {
	BEGIN_RCPP
		Rcpp::RObject rcpp_result_gen;
	Rcpp::RNGScope rcpp_rngScope_gen;
	Rcpp::traits::input_parameter< arma::mat& >::type Xr(XrSEXP);
	Rcpp::traits::input_parameter< arma::mat& >::type Xr2(Xr2SEXP);
	Rcpp::traits::input_parameter< arma::colvec& >::type yr(yrSEXP);
	Rcpp::traits::input_parameter< arma::colvec& >::type yr2(yr2SEXP);
	Rcpp::traits::input_parameter< SEXP >::type opts(optsSEXP);
	rcpp_result_gen = Rcpp::wrap(LogisMVS4GVB(Xr, Xr2, yr, yr2, opts));
	return rcpp_result_gen;
	END_RCPP
}

// LogisWFMVS2GVB
Rcpp::List LogisWFMVS2GVB(arma::mat& Xr, arma::mat& Zr, arma::colvec& yr, SEXP opts);
RcppExport SEXP LPG_LogisWFMVS2GVB(SEXP XrSEXP, SEXP ZrSEXP, SEXP yrSEXP, SEXP optsSEXP) {
	BEGIN_RCPP
		Rcpp::RObject rcpp_result_gen;
	Rcpp::RNGScope rcpp_rngScope_gen;
	Rcpp::traits::input_parameter< arma::mat& >::type Xr(XrSEXP);
	Rcpp::traits::input_parameter< arma::mat& >::type Zr(ZrSEXP);
	Rcpp::traits::input_parameter< arma::colvec& >::type yr(yrSEXP);
	Rcpp::traits::input_parameter< SEXP >::type opts(optsSEXP);
	rcpp_result_gen = Rcpp::wrap(LogisWFMVS2GVB(Xr, Zr, yr, opts));
	return rcpp_result_gen;
	END_RCPP
}

// LogisWFMVS4GVB
Rcpp::List LogisWFMVS4GVB(arma::mat& Xr, arma::mat& Xr2, arma::mat& Zr, arma::mat& Zr2, arma::colvec& yr, arma::colvec& yr2, SEXP opts);
RcppExport SEXP LPG_LogisWFMVS4GVB(SEXP XrSEXP, SEXP Xr2SEXP, SEXP ZrSEXP, SEXP Zr2SEXP, SEXP yrSEXP, SEXP yr2SEXP, SEXP optsSEXP) {
	BEGIN_RCPP
		Rcpp::RObject rcpp_result_gen;
	Rcpp::RNGScope rcpp_rngScope_gen;
	Rcpp::traits::input_parameter< arma::mat& >::type Xr(XrSEXP);
	Rcpp::traits::input_parameter< arma::mat& >::type Zr(ZrSEXP);
	Rcpp::traits::input_parameter< arma::colvec& >::type yr(yrSEXP);
	Rcpp::traits::input_parameter< arma::mat& >::type Xr2(Xr2SEXP);
	Rcpp::traits::input_parameter< arma::mat& >::type Zr2(Zr2SEXP);
	Rcpp::traits::input_parameter< arma::colvec& >::type yr2(yr2SEXP);
	Rcpp::traits::input_parameter< SEXP >::type opts(optsSEXP);
	rcpp_result_gen = Rcpp::wrap(LogisWFMVS4GVB(Xr, Xr2, Zr, Zr2, yr, yr2, opts));
	return rcpp_result_gen;
	END_RCPP
}


// ILogisMVS2GVB
Rcpp::List ILogisMVS2GVB(string stringname, SEXP opts);
RcppExport SEXP LPG_ILogisMVS2GVB(SEXP stringnameSEXP, SEXP optsSEXP) {
	BEGIN_RCPP
		Rcpp::RObject rcpp_result_gen;
	Rcpp::RNGScope rcpp_rngScope_gen;
	Rcpp::traits::input_parameter< string >::type stringname(stringnameSEXP);
	Rcpp::traits::input_parameter< SEXP >::type opts(optsSEXP);
	rcpp_result_gen = Rcpp::wrap(ILogisMVS2GVB(stringname, opts));
	return rcpp_result_gen;
	END_RCPP
}

// ILogisMVS4GVB
Rcpp::List ILogisMVS4GVB(string stringname, string stringname2, SEXP opts);
RcppExport SEXP LPG_ILogisMVS4GVB(SEXP stringnameSEXP, SEXP stringname2SEXP, SEXP optsSEXP) {
	BEGIN_RCPP
		Rcpp::RObject rcpp_result_gen;
	Rcpp::RNGScope rcpp_rngScope_gen;
	Rcpp::traits::input_parameter< string >::type stringname(stringnameSEXP);
	Rcpp::traits::input_parameter< string >::type stringname2(stringname2SEXP);
	Rcpp::traits::input_parameter< SEXP >::type opts(optsSEXP);
	rcpp_result_gen = Rcpp::wrap(ILogisMVS4GVB(stringname, stringname2, opts));
	return rcpp_result_gen;
	END_RCPP
}


// ILogisWFMVS2GVB
Rcpp::List ILogisWFMVS2GVB(string stringname, arma::mat & Zr, SEXP opts);
RcppExport SEXP LPG_ILogisWFMVS2GVB(SEXP stringnameSEXP, SEXP ZrSEXP, SEXP optsSEXP) {
	BEGIN_RCPP
	Rcpp::RObject rcpp_result_gen;
	Rcpp::RNGScope rcpp_rngScope_gen;
	Rcpp::traits::input_parameter< string >::type stringname(stringnameSEXP);
	Rcpp::traits::input_parameter< arma::mat& >::type Zr(ZrSEXP);
	Rcpp::traits::input_parameter< SEXP >::type opts(optsSEXP);
	rcpp_result_gen = Rcpp::wrap(ILogisWFMVS2GVB(stringname, Zr, opts));
	return rcpp_result_gen;
	END_RCPP
}


// ILogisWFMVS4GVB
Rcpp::List ILogisWFMVS4GVB(string stringname, arma::mat & Zr, string stringname2, arma::mat & Zr2, SEXP opts);
RcppExport SEXP LPG_ILogisWFMVS4GVB(SEXP stringnameSEXP, SEXP ZrSEXP, SEXP stringname2SEXP, SEXP Zr2SEXP, SEXP optsSEXP) {
	BEGIN_RCPP
	Rcpp::RObject rcpp_result_gen;
	Rcpp::RNGScope rcpp_rngScope_gen;
	Rcpp::traits::input_parameter< string >::type stringname(stringnameSEXP);
	Rcpp::traits::input_parameter< string >::type stringname2(stringname2SEXP);
	Rcpp::traits::input_parameter< arma::mat& >::type Zr(ZrSEXP);
	Rcpp::traits::input_parameter< arma::mat& >::type Zr2(Zr2SEXP);
	Rcpp::traits::input_parameter< SEXP >::type opts(optsSEXP);
	rcpp_result_gen = Rcpp::wrap(ILogisWFMVS4GVB(stringname, Zr, stringname2, Zr2, opts));
	return rcpp_result_gen;
	END_RCPP
}

/////////////////////////////////////////////////////////////////////////////////////
/////                                                                           /////
/////                         read plink function                               /////
/////                                                                           /////
/////////////////////////////////////////////////////////////////////////////////////

// ReadPlinkRcpp
Rcpp::List ReadPlinkRcpp(string stringname);
RcppExport SEXP LPG_ReadPlinkRcpp(SEXP stringnameSEXP) {
	BEGIN_RCPP
		Rcpp::RObject rcpp_result_gen;
	Rcpp::RNGScope rcpp_rngScope_gen;
	Rcpp::traits::input_parameter< string >::type stringname(stringnameSEXP);
	rcpp_result_gen = Rcpp::wrap(ReadPlinkRcpp(stringname));
	return rcpp_result_gen;
	END_RCPP
}