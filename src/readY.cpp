#include <RcppArmadillo.h>
#include "readY.hpp"
#include <iostream>
#include <fstream>
#include <Rcpp.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

ObjY ReadYFromFile(string stringname) {
	string famfile = stringname;
	famfile += ".fam";
	int N = getLineNum(famfile);
	arma::colvec fammat = arma::zeros(N, 1);
	fstream file;
	char arr[100];
	strcpy(arr, famfile.c_str());
	file.open(arr, ios::in);
	string s;
	unsigned i = 1;
	while (file >> s)
	{
		if (i % 6 == 0)
		{
			//string->double
			double a;
			a = atof(s.c_str());
			fammat[i / 6 - 1] = a;
		}
		i++;
	}
	arma::colvec Y = fammat;
	ObjY obj;
	obj.Y = Y;
	return obj;
}