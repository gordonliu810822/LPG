#ifndef LogisModel_hpp
#define LogisModel_hpp

#include <RcppArmadillo.h>
#include<iostream>
#include<fstream>
#include <stdio.h>
#include "plinkfun.hpp"



using namespace std;
using namespace arma;


class Options{
public:
	Options(){
		this->max_iter = 1e6;
		this->dispF = 1;
		this->display_gap = 10;
		this->epsStopLogLik = 1e-5;
		this->constraintalpha = 0;
	}

	Options(int max_iter, int dispF, int display_gap, double epsStopLogLik){
		this->max_iter = max_iter;
		this->dispF = dispF;
		this->display_gap = display_gap;
		this->epsStopLogLik = epsStopLogLik;
	}

	Options(int max_iter, int dispF, int display_gap, double epsStopLogLik, int constraintalpha){
		this->max_iter = max_iter;
		this->dispF = dispF;
		this->display_gap = display_gap;
		this->epsStopLogLik = epsStopLogLik;
		this->constraintalpha = constraintalpha;
	}
	int max_iter;
	int dispF;
	int display_gap;
	double epsStopLogLik;
	int constraintalpha;

};

struct ObjLogisMVS2GVB{
	arma::colvec vardist_gamma;
	arma::colvec vardist_mu;
	arma::colvec vardist_sigma2beta;
	arma::colvec sigma2beta;
	arma::colvec alpha;
	arma::rowvec Lq;
	double beta0;
};

struct ObjLogisMVS4GVB{
	arma::mat vardist_gamma;
	arma::mat vardist_mu;
	arma::mat vardist_sigma2beta;
	arma::colvec sigma2beta;
	arma::rowvec alpha;
	arma::rowvec Lq;
	arma::mat beta0;
};

struct ObjLogisWFMVS2GVB{
	arma::mat vardist_gamma;
	arma::mat vardist_mu;
	arma::mat vardist_sigma2beta;
	arma::colvec sigma2beta;
	arma::rowvec alpha;
	arma::rowvec Lq;
	arma::colvec u;
};

struct ObjLogisWFMVS4GVB{
	arma::mat vardist_gamma;
	arma::mat vardist_mu;
	arma::mat vardist_sigma2beta;
	arma::colvec sigma2beta;
	arma::rowvec alpha;
	arma::rowvec Lq;
	arma::mat u;
};

ObjLogisMVS2GVB rcpparma_LogisMVS2GVB(arma::mat & Xr, arma::colvec & yr, Options* opts);

ObjLogisMVS4GVB rcpparma_LogisMVS4GVB(arma::mat& Xr, arma::mat& Xr2, arma::colvec& yr, arma::colvec& yr2, Options* opts);

ObjLogisWFMVS2GVB rcpparma_LogisWFMVS2GVB(arma::mat & Xr, arma::mat & Zr, arma::colvec & yr, Options* opts);

ObjLogisWFMVS4GVB rcpparma_LogisWFMVS4GVB(arma::mat & Xr, arma::mat& Xr2, arma::mat & Zr, arma::mat & Zr2, arma::colvec & yr, arma::colvec & yr2, Options* opts);
#endif
