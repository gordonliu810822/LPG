//
//  readPlink.cpp
//  ReadPlinkGit


#ifndef readY_hpp
#define readY_hpp

#include <RcppArmadillo.h>
#include<iostream>
#include<fstream>
#include <stdio.h>
#include "plinkfun.hpp"



using namespace std;
using namespace arma;

struct ObjY{
    arma::vec Y;  
};


ObjY ReadYFromFile(string stringname);


#endif
