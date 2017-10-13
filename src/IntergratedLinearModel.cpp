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


// LinearMVS2GVB
Rcpp::List ILinearMVS2GVB(string stringname, SEXP opts = R_NilValue){

	Rcout << "begin to read genetype and phenotype" << endl;
	ObjXY obj_XY = ReadDataFromFile(stringname);

	arma::Mat<unsigned>* X = &obj_XY.X;
	arma::mat standardX = StandardX(X);

	Options* lp_opt = NULL;
	if (!Rf_isNull(opts)){
		Rcpp::List opt(opts);
		lp_opt = new Options(opt["max_iter"], opt["dispF"], opt["display_gap"], opt["epsStopLogLik"]);
	}

	if (Rf_isNull(opts)){
		lp_opt = new Options();
	}

	Rcout << "begin to fit LPG model" << endl;
	ObjLinearMVS2GVB obj = rcpparma_LinearMVS2GVB(standardX, obj_XY.Y, lp_opt);

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



// LinearMVS4GVB
Rcpp::List ILinearMVS4GVB(string stringname1, string stringname2, SEXP opts = R_NilValue){

	Rcout << "begin to read genetype and phenotype" << endl;
	ObjXY obj_XY1 = ReadDataFromFile(stringname1);
	ObjXY obj_XY2 = ReadDataFromFile(stringname2);

	arma::Mat<unsigned>* X1 = &obj_XY1.X;
	arma::mat standardX1 = StandardX(X1);
	arma::Mat<unsigned>* X2 = &obj_XY2.X;
	arma::mat standardX2 = StandardX(X2);

	Options* lp_opt = NULL;
	if (!Rf_isNull(opts)){
		Rcpp::List opt(opts);
		lp_opt = new Options(opt["max_iter"], opt["dispF"], opt["display_gap"], opt["epsStopLogLik"], opt["constraintalpha"]);
	}

	if (Rf_isNull(opts)){
		lp_opt = new Options();
	}

	Rcout << "begin to fit LPG model" << endl;
	ObjLinearMVS4GVB obj = rcpparma_LinearMVS4GVB(standardX1, standardX2, obj_XY1.Y, obj_XY2.Y, lp_opt); //obj is c++ language

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