#include <RcppArmadillo.h>
#include "readPlink.hpp"
#include "readY.hpp"


// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//



ObjXY ReadDataFromFile(string stringname) {
  string famfile = stringname;
  famfile += ".fam";
  string bimfile = stringname;
  bimfile += ".bim";
  int N =  getLineNum(famfile);
  int P =  getLineNum(bimfile);
  arma::mat bimmat;
  bimmat.load(bimfile,arma::csv_ascii);
  arma::vec c1 = bimmat.col(0);
  cout << bimmat.size() << endl;
  arma::vec index(22);
  double sum2 = 0;
  for(int i=1; i <= 22; i++ ){
    index[i-1] = sum(c1 == i);
    sum2 += index[i-1];
    cout <<"number of snps in chromsome "<<i <<" =" << index[i-1] << endl;
  }
  unsigned* X = new unsigned[ N * P];
  readPlink(stringname,N, P, X);
  ObjY Obj_Y = ReadYFromFile(stringname);
  arma::Mat<unsigned>* Xdata = new arma::Mat<unsigned>(X, N, P, false,false);
  ObjXY obj;
  obj.X = *Xdata;
  obj.Y = Obj_Y.Y;
  return obj;
}
