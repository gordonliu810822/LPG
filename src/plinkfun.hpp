#ifndef plinkfun_hpp
#define plinkfun_hpp

#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <string.h>
using namespace std;


void readPlink(string stringname,int N, int P, unsigned* X);
int getLineNum(string filename);

#endif /* plinkfun_hpp */
