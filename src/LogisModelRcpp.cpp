#include "RcppArmadillo.h"
#include <Rcpp.h>
#include "LogisModel.hpp"
#include "readPlink.hpp"
#include "readY.hpp"
#include "StandardizeData.hpp"

using namespace Rcpp;
using namespace std;
using namespace arma;


// [[Rcpp::export]]
Rcpp::List LogisMVS2GVB(arma::mat& Xr, arma::colvec& yr, SEXP opts = R_NilValue) {

	int n = Xr.n_rows;
	int p = Xr.n_cols;
	
	arma::mat X(Xr.begin(), n, p, false);// pointer
	arma::colvec y(yr.begin(), yr.size(), false);

	Options* lp_opt = NULL;
	if (!Rf_isNull(opts)){
		Rcpp::List opt(opts);
		lp_opt = new Options(opt["max_iter"], opt["dispF"], opt["display_gap"], opt["epsStopLogLik"]);
	}

	if (Rf_isNull(opts)){
		lp_opt = new Options();
	}

	ObjLogisMVS2GVB obj = rcpparma_LogisMVS2GVB(X, y, lp_opt);

	List ret;
	ret["vardist_gamma"] = Rcpp::wrap(obj.vardist_gamma);
	ret["vardist_mu"] = Rcpp::wrap(obj.vardist_mu);
	ret["vardist_sigma2beta"] = Rcpp::wrap(obj.vardist_sigma2beta);
	ret["sigma2beta"] = Rcpp::wrap(obj.sigma2beta);
	ret["alpha"] = Rcpp::wrap(obj.alpha);
	ret["Lq"] = Rcpp::wrap(obj.Lq);
	ret["u"] = Rcpp::wrap(obj.beta0);

	return ret;	
}



Rcpp::List LogisMVS4GVB(arma::mat& Xr, arma::mat& Xr2, arma::colvec& yr, arma::colvec& yr2, SEXP opts = R_NilValue) {

	int n1 = Xr.n_rows;
	int p = Xr.n_cols;
	int n2 = Xr2.n_rows;

	arma::mat X(Xr.begin(), n1, p, false);// pointer
	arma::colvec y(yr.begin(), yr.size(), false);
	arma::mat X2(Xr2.begin(), n2, p, false);// pointer
	arma::colvec y2(yr2.begin(), yr2.size(), false);

	Options* lp_opt = NULL;
	if (!Rf_isNull(opts)){
		Rcpp::List opt(opts);
		lp_opt = new Options(opt["max_iter"], opt["dispF"], opt["display_gap"], opt["epsStopLogLik"], opt["constraintalpha"]);
	}

	if (Rf_isNull(opts)){
		lp_opt = new Options();
	}

	ObjLogisMVS4GVB obj = rcpparma_LogisMVS4GVB(X, X2, y, y2, lp_opt);

	List ret;
	ret["vardist_gamma"] = Rcpp::wrap(obj.vardist_gamma);
	ret["vardist_mu"] = Rcpp::wrap(obj.vardist_mu);
	ret["vardist_sigma2beta"] = Rcpp::wrap(obj.vardist_sigma2beta);
	ret["sigma2beta"] = Rcpp::wrap(obj.sigma2beta);
	ret["alpha"] = Rcpp::wrap(obj.alpha);
	ret["Lq"] = Rcpp::wrap(obj.Lq);
	ret["u"] = Rcpp::wrap(obj.beta0);

	return ret;
}

Rcpp::List LogisWFMVS2GVB(arma::mat& Xr, arma::mat& Zr, arma::colvec& yr, SEXP opts = R_NilValue) {

	int n = Xr.n_rows;
	int p = Xr.n_cols;
	int q = Zr.n_cols;

	arma::mat X(Xr.begin(), n, p, false);// pointer
	arma::mat Z(Zr.begin(), n, q, false);// pointer
	arma::colvec y(yr.begin(), yr.size(), false);

	Options* lp_opt = NULL;
	if (!Rf_isNull(opts)){
		Rcpp::List opt(opts);
		lp_opt = new Options(opt["max_iter"], opt["dispF"], opt["display_gap"], opt["epsStopLogLik"]);
	}

	if (Rf_isNull(opts)){
		lp_opt = new Options();
	}

	ObjLogisWFMVS2GVB obj = rcpparma_LogisWFMVS2GVB(X, Z, y, lp_opt);

	List ret;
	ret["vardist_gamma"] = Rcpp::wrap(obj.vardist_gamma);
	ret["vardist_mu"] = Rcpp::wrap(obj.vardist_mu);
	ret["vardist_sigma2beta"] = Rcpp::wrap(obj.vardist_sigma2beta);
	ret["sigma2beta"] = Rcpp::wrap(obj.sigma2beta);
	ret["alpha"] = Rcpp::wrap(obj.alpha);
	ret["Lq"] = Rcpp::wrap(obj.Lq);
	ret["u"] = Rcpp::wrap(obj.u);

	return ret;
}

Rcpp::List LogisWFMVS4GVB(arma::mat& Xr, arma::mat& Xr2, arma::mat& Zr, arma::mat& Zr2, arma::colvec& yr, arma::colvec& yr2, SEXP opts = R_NilValue) {



	int n1 = Xr.n_rows;
	int n2 = Xr2.n_rows;
	int p = Xr.n_cols;
	int q1 = Zr.n_cols;
	int q2 = Zr2.n_cols;

	arma::mat X(Xr.begin(), n1, p, false);// pointer
	arma::mat X2(Xr2.begin(), n2, p, false);// pointer
	arma::mat Z(Zr.begin(), n1, q1, false);// pointer
	arma::mat Z2(Zr2.begin(), n2, q2, false);// pointer
	arma::colvec y(yr.begin(), yr.size(), false);
	arma::colvec y2(yr2.begin(), yr2.size(), false);

	Options* lp_opt = NULL;
	if (!Rf_isNull(opts)){
		Rcpp::List opt(opts);
		lp_opt = new Options(opt["max_iter"], opt["dispF"], opt["display_gap"], opt["epsStopLogLik"], opt["constraintalpha"]);
	}

	if (Rf_isNull(opts)){
		lp_opt = new Options();
	}

	ObjLogisWFMVS4GVB obj = rcpparma_LogisWFMVS4GVB(X, X2, Z, Z2, y, y2, lp_opt);

	List ret;
	ret["vardist_gamma"] = Rcpp::wrap(obj.vardist_gamma);
	ret["vardist_mu"] = Rcpp::wrap(obj.vardist_mu);
	ret["vardist_sigma2beta"] = Rcpp::wrap(obj.vardist_sigma2beta);
	ret["sigma2beta"] = Rcpp::wrap(obj.sigma2beta);
	ret["alpha"] = Rcpp::wrap(obj.alpha);
	ret["Lq"] = Rcpp::wrap(obj.Lq);
	ret["u"] = Rcpp::wrap(obj.u);

	return ret;
}