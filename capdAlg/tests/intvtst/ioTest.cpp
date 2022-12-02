/*
 * ioTest.cpp
 *
 *  Created on: Jul 28, 2012
 *      Author: kapela
 */

#include <iostream>
#include <fstream>
#include "capd/intervals/lib.h"
typedef capd::DInterval DInterval;
using namespace std;

bool checkEqual(const DInterval & a, const DInterval & b){
	if(a == b)
		return true;
	std::cout << "\n Intervals not equal  a = " << a << " b = " << b << std::endl;
	return false;
}
bool checkContained(const DInterval & a, const DInterval & b){
	if(a == b)
		return true;
        if(b.contains(a)){
        	std::cout << "\n Warning: we wrote " << a << " and read " << b << std::endl;
                return true;
        }
        std::cout << "\nError : we wrote " << a << " and read " << b << std::endl;
	return false;
}

bool testIO(const DInterval & x){
  std::cout <<  "iotest for : " << x << std::endl;
	std::ofstream out("tmp.txt");
	out.precision(17);
	out << x << endl;
	hexWrite(out, x) << endl;
	bitWrite(out, x) <<endl;
	out.close();
	DInterval dec, hex, bit, bin;
	std::ifstream in("tmp.txt");
	in >> dec;
	hexRead(in, hex);
	bitRead(in, bit);
  in.close();
	out.open("tmpb.txt",std::ios::binary|std::ios::out);
	binWrite(out,x);
	out.close();
	in.open("tmpb.txt",std::ios::binary|std::ios::in);
  binRead(in, bin);
  in.close();
    
  return  checkContained(x, dec) && checkEqual(x, hex) &&
    	   checkEqual(x, bit) && checkEqual(x, bin);
}

int main() {

	DInterval x(1.0), y(0.5), z(3.0), a(-3.0, 8.0);
  x /= 3;
	z /= 7;
	a /= 100;
  double xx= 1./3.;
  DInterval b(xx);
  std::cout.precision(17);
  bool result  = testIO(x) && testIO(y)
    		    && testIO(z) && testIO(a) && testIO(b);
	return result? 0 : -1;
}



