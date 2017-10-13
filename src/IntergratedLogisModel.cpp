#include "RcppArmadillo.h"
#include <Rcpp.h>
#include "LogisModel.hpp"
#include "readPlink.hpp"
#include "readY.hpp"
#include "StandardizeData.hpp"

using namespace Rcpp;
using namespace std;
using namespace arma;

// ILogisMVS2GVB
Rcpp::List ILogisMVS2GVB(string stringname, SEXP opts = R_NilValue){

	Rcout << "begin to read genetype and phenotype" << endl;
	ObjXY obj1 = ReadDataFromFile(stringname);

	Rcout << "begin to acquire the number of sample" << endl;
	string famfileqc3 = stringname;
	famfileqc3 += ".fam";
	int Nqc3 = getLineNum(famfileqc3);
	Rcout << "number of sample is " << Nqc3 << endl;


	//cout << "begin to scale the genetype" << endl;
	arma::Mat<unsigned>* X = &obj1.X;
	arma::mat standardX = StandardX(X);

	arma::colvec Y = arma::zeros(Nqc3);
	for (int i = 0; i < Nqc3; i++)
	{
		if (obj1.Y(i) >1.5)
			Y(i) = 1;
		else
			Y(i) = -1;
	}

	Options* lp_opt = NULL;
	if (!Rf_isNull(opts)){
		Rcpp::List opt(opts);
		lp_opt = new Options(opt["max_iter"], opt["dispF"], opt["display_gap"], opt["epsStopLogLik"]);
	}

	if (Rf_isNull(opts)){
		lp_opt = new Options();
	}

	Rcout << "begin to fit LPG model" << endl;
	ObjLogisMVS2GVB obj = rcpparma_LogisMVS2GVB(standardX, Y, lp_opt);


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

// ILogisMVS4GVB
Rcpp::List ILogisMVS4GVB(string stringname, string stringname2, SEXP opts = R_NilValue){
	
	Rcout << "begin to read genetype and phenotype" << endl;
	ObjXY obj1 = ReadDataFromFile(stringname);
	ObjXY obj2 = ReadDataFromFile(stringname2);

	Rcout << "begin to acquire the number of sample" << endl;
	string famfileqc3 = stringname;
	famfileqc3 += ".fam";
	int Nqc3 = getLineNum(famfileqc3);
	Rcout <<  "number of sample in study 1 is "<< Nqc3 << endl;

	string X2famfileqc3 = stringname2;
	X2famfileqc3 += ".fam";
	int X2Nqc3 = getLineNum(X2famfileqc3);
	Rcout << "number of sample in study 2 is " << X2Nqc3 << endl;


	arma::colvec Y = arma::zeros(Nqc3);
	for (int i = 0; i < Nqc3; i++)
	{
		if (obj1.Y(i) >1.5)
			Y(i) = 1;
		else
			Y(i) = -1;
	}

	arma::colvec Y2 = arma::zeros(X2Nqc3);
	for (int i = 0; i < X2Nqc3; i++)
	{
		if (obj2.Y(i) >1.5)
			Y2(i) = 1;
		else
			Y2(i) = -1;
	}

	//cout << "begin to scale the genetype" << endl;
	arma::Mat<unsigned>* X = &obj1.X;
	arma::mat standardX = StandardX(X);

	arma::Mat<unsigned>* X2 = &obj2.X;
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
	ObjLogisMVS4GVB obj = rcpparma_LogisMVS4GVB(standardX, standardX2, Y, Y2, lp_opt);


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

// ILogisWFMVS2GVB
Rcpp::List ILogisWFMVS2GVB(string stringname, arma::mat& Zr, SEXP opts = R_NilValue){

	Rcout << "begin to read genetype and phenotype" << endl;
	ObjXY obj1 = ReadDataFromFile(stringname);

	Rcout << "begin to acquire the number of sample" << endl;
	string famfileqc3 = stringname;
	famfileqc3 += ".fam";
	int Nqc3 = getLineNum(famfileqc3);
	Rcout << "number of sample is " << Nqc3 << endl;


	//cout << "begin to scale the genetype" << endl;
	arma::Mat<unsigned>* X = &obj1.X;
	arma::mat standardX = StandardX(X);

	arma::colvec Y = arma::zeros(Nqc3);
	for (int i = 0; i < Nqc3; i++)
	{
		if (obj1.Y(i) >1.5)
			Y(i) = 1;
		else
			Y(i) = -1;
	}

	Rcout << "begin to LPG analyze" << endl;
	int n = Zr.n_rows;
	int q = Zr.n_cols;
	arma::mat Z(Zr.begin(), n, q, false);// pointer

	Options* lp_opt = NULL;
	if (!Rf_isNull(opts)){
		Rcpp::List opt(opts);
		lp_opt = new Options(opt["max_iter"], opt["dispF"], opt["display_gap"], opt["epsStopLogLik"]);
	}

	if (Rf_isNull(opts)){
		lp_opt = new Options();
	}

	ObjLogisWFMVS2GVB obj = rcpparma_LogisWFMVS2GVB(standardX, Z, Y, lp_opt);


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


// ILogisWFMVS4GVB
Rcpp::List ILogisWFMVS4GVB(string stringname, arma::mat& Zr, string stringname2, arma::mat& Zr2, SEXP opts = R_NilValue){

	Rcout << "begin to read genetype and phenotype" << endl;
	ObjXY obj1 = ReadDataFromFile(stringname);
	ObjXY obj2 = ReadDataFromFile(stringname2);

	Rcout << "begin to acquire the number of sample" << endl;
	string famfileqc3 = stringname;
	famfileqc3 += ".fam";
	int Nqc3 = getLineNum(famfileqc3);
	Rcout << "number of sample in study 1 is " << Nqc3 << endl;

	string X2famfileqc3 = stringname2;
	X2famfileqc3 += ".fam";
	int X2Nqc3 = getLineNum(X2famfileqc3);
	Rcout << "number of sample in study 2 is " << X2Nqc3 << endl;


	arma::colvec Y = arma::zeros(Nqc3);
	for (int i = 0; i < Nqc3; i++)
	{
		if (obj1.Y(i) >1.5)
			Y(i) = 1;
		else
			Y(i) = -1;
	}

	arma::colvec Y2 = arma::zeros(X2Nqc3);
	for (int i = 0; i < X2Nqc3; i++)
	{
		if (obj2.Y(i) >1.5)
			Y2(i) = 1;
		else
			Y2(i) = -1;
	}

	//cout << "begin to scale the genetype" << endl;
	arma::Mat<unsigned>* X = &obj1.X;
	arma::mat standardX = StandardX(X);

	arma::Mat<unsigned>* X2 = &obj2.X;
	arma::mat standardX2 = StandardX(X2);

	Rcout << "begin to fit LPG model" << endl;
	int n1 = Zr.n_rows;
	int q1 = Zr.n_cols;
	int n2 = Zr2.n_rows;
	int q2 = Zr2.n_cols;
	arma::mat Z(Zr.begin(), n1, q1, false);// pointer
	arma::mat Z2(Zr2.begin(), n2, q2, false);// pointer

	Options* lp_opt = NULL;
	if (!Rf_isNull(opts)){
		Rcpp::List opt(opts);
		lp_opt = new Options(opt["max_iter"], opt["dispF"], opt["display_gap"], opt["epsStopLogLik"], opt["constraintalpha"]);
	}

	if (Rf_isNull(opts)){
		lp_opt = new Options();
	}

	ObjLogisWFMVS4GVB obj = rcpparma_LogisWFMVS4GVB(standardX, standardX2, Z, Z2, Y, Y2, lp_opt);


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