#ifndef StandardizeData_hpp
#define StandardizeData_hpp

#include <RcppArmadillo.h>
#include<iostream>
#include<fstream>
#include <stdio.h>
#include "plinkfun.hpp"


using namespace std;
using namespace arma;

arma::mat StandardX(arma::Mat<unsigned>* X);

#endif
