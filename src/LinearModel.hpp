#ifndef LinearModel_hpp
#define LinearModel_hpp

#include <RcppArmadillo.h>
#include<iostream>
#include<fstream>
#include <stdio.h>
#include "plinkfun.hpp"
#include "LogisModel.hpp"


using namespace std;
using namespace arma;


struct ObjLinearMVS2GVB{
	arma::colvec vardist_gamma;
	arma::colvec vardist_mu;
	arma::colvec vardist_sigma2beta;
	arma::colvec sigma2beta;
	arma::colvec alpha;
	arma::colvec sigma2e;
	arma::rowvec Lq;
};


struct ObjLinearMVS4GVB{
	arma::mat vardist_gamma;
	arma::mat vardist_mu;
	arma::mat vardist_sigma2beta;
	arma::colvec sigma2beta;
	arma::rowvec alpha;
	arma::colvec sigma2e;
	arma::rowvec Lq;
};


ObjLinearMVS2GVB rcpparma_LinearMVS2GVB(arma::mat & Xr, arma::colvec & yr, Options* opts);
ObjLinearMVS4GVB rcpparma_LinearMVS4GVB(arma::mat & Xr, arma::mat & Xr2, arma::colvec & yr, arma::colvec & yr2, Options* opts);

#endif
