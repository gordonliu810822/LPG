#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "readPlink.hpp"

using namespace Rcpp;
using namespace arma;
using namespace std;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


// [[Rcpp::export]]
Rcpp::List ReadPlinkRcpp(string stringname) {
	//index of R vector substract 1 for the C indices
	//uvec	 = 	ucolvec	 = 	Col<uword>
	//uword is a typedef for an unsigned integer type
	ObjXY obj = ReadDataFromFile(stringname);
	List ret;
	ret["X"] = Rcpp::wrap(obj.X);
	ret["Y"] = Rcpp::wrap(obj.Y);
	return ret;
}



