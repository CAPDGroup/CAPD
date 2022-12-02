/**
 *  @file multiprecExample.cpp
 *
 *  you can compile with:
 *    g++ -o multiprecExample -D__HAVE_MPFR__ -s -O2 -Wall -frounding-math -I. -ICAPD_DIR/include  multiprecExample.cpp CAPD_DIR/lib/libcapd.a  -lmpfr -lgmp
 *  where
 *    CAPD_DIR is a root directory of the CAPD library
 *  Created on: 2009-12-01
 *      Author: kapela
 */

#include <iostream>
#include "capd/map/Function.hpp"
#include "capd/capdlib.h"
#include "capd/mpcapdlib.h"
using namespace std;
using namespace capd;

int main(){
  // We set precision to 100 mantisa bits (about 30 decimal digits)
  MpFloat::setDefaultPrecision(100);
  // We set MpFloat to round numbers to the nearest representable
  MpFloat::setDefaultRndMode(MpFloat::RoundNearest);
  MpFloat a = 1,
      b("1.2345678901234567890123456789"),       // string allows to input any number of digits
      c(1.12112, MpFloat::RoundUp, 10),  // use only 10 mantisa bits and round up the initial value 1.2112
      d(1.234567890123456789),                   // only 15 digits is used to initialize (1.234... is converted to double)
      e(1.234567890123456789L);                  // here we say (L suffix) that parameter is of a long double type
  cout << "Precision used (in mantisa bits): " << MpFloat::getDefaultPrecision() <<"\n";
  cout << "\n a = " << a << "\n b = " << b << "\n c = " << c << "\n d = " << d << "\n e = " << e ;

  MpVector v(3);
  v[0] = a;  v[1] = b; v[2] = c;
  cout << "\n\n v = " << v;

  MpFunction f("var:a,b,c;fun:a+b+c;");
  cout << "\n f(a,b,c) = " << f(v);

  MpInterval  ia(a),
              ib(b,c),
              ic("-1.2345678901234567890","2020.202020202020202002");
  cout << "\n\n ia = " << ia << "\n ib = " << ib << "\n ic " << ic << "\n";
  return 0;
}
