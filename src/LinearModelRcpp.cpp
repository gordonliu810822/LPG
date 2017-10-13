#include "RcppArmadillo.h"
#include <Rcpp.h>
#include "LinearModel.hpp"
#include "readPlink.hpp"
#include "readY.hpp"
#include "StandardizeData.hpp"
#include "LogisModel.hpp"

using namespace Rcpp;
using namespace std;
using namespace arma;


// [[Rcpp::export]]
Rcpp::List LinearMVS4GVB(mat&  Xr, mat&  Xr2, colvec&  yr, colvec&  yr2, SEXP opts = R_NilValue) {
	
	int n1 = Xr.n_rows;
	int n2 = Xr2.n_rows;
	int p = Xr.n_cols;

	mat X1(Xr.begin(), n1, p, false);// pointer
	colvec y1(yr.begin(), yr.size(), false);
	mat X2(Xr2.begin(), n2, p, false);// pointer
	colvec y2(yr2.begin(), yr2.size(), false);

	Options* lp_opt = NULL;
	if (!Rf_isNull(opts)){
		Rcpp::List opt(opts);
		lp_opt = new Options(opt["max_iter"], opt["dispF"], opt["display_gap"], opt["epsStopLogLik"], opt["constraintalpha"]);
	}

	if (Rf_isNull(opts)){
		lp_opt = new Options();
	}

	ObjLinearMVS4GVB obj = rcpparma_LinearMVS4GVB(X1, X2, y1, y2, lp_opt); //obj is c++ language

	List ret;
	ret["vardist_gamma"] = Rcpp::wrap(obj.vardist_gamma);
	ret["vardist_mu"] = Rcpp::wrap(obj.vardist_mu);
	ret["vardist_sigma2beta"] = Rcpp::wrap(obj.vardist_sigma2beta);
	ret["sigma2beta"] = Rcpp::wrap(obj.sigma2beta);
	ret["alpha"] = Rcpp::wrap(obj.alpha);
	ret["sigma2e"] = Rcpp::wrap(obj.sigma2e);
	ret["Lq"] = Rcpp::wrap(obj.Lq);

	return ret;
}


// [[Rcpp::export]]
Rcpp::List LinearMVS2GVB(mat& Xr, colvec& yr, SEXP opts = R_NilValue) {

	int n = Xr.n_rows;
	int p = Xr.n_cols;
	

	mat X(Xr.begin(), n, p, false);// pointer
	colvec y(yr.begin(), yr.size(), false);

	Options* lp_opt = NULL;
	if (!Rf_isNull(opts)){
		Rcpp::List opt(opts);
		lp_opt = new Options(opt["max_iter"], opt["dispF"], opt["display_gap"], opt["epsStopLogLik"]);
	}

	if (Rf_isNull(opts)){
		lp_opt = new Options();
	}

	ObjLinearMVS2GVB obj = rcpparma_LinearMVS2GVB(X, y, lp_opt);

	List ret;
	ret["vardist_gamma"] = Rcpp::wrap(obj.vardist_gamma);
	ret["vardist_mu"] = Rcpp::wrap(obj.vardist_mu);
	ret["vardist_sigma2beta"] = Rcpp::wrap(obj.vardist_sigma2beta);
	ret["sigma2beta"] = Rcpp::wrap(obj.sigma2beta);
	ret["alpha"] = Rcpp::wrap(obj.alpha);
	ret["sigma2e"] = Rcpp::wrap(obj.sigma2e);
	ret["Lq"] = Rcpp::wrap(obj.Lq);

	return ret;

	
}
