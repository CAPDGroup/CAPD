
/*!
 * @file mpExample1.cpp
 * Example how to use MpReal
 *
 * @author kapela  @date 2010-01-07
 *
 */

// Copyright (C) CAPD group
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.


#include <iostream>
#include "capd/capdlib.h"
#include "capd/mpcapdlib.h"
using namespace std;
using namespace capd;

int main(){
  /// We set precision to 100 mantisa bits (about 30 decimal digits)
  MpFloat::setDefaultPrecision(100);
  /// We set MpFloat to round numbers to the nearest representable
  MpFloat::setDefaultRndMode(MpFloat::RoundNearest);
  MpFloat a = 1,
      b("1.2345678901234567890123456789"),       // string allows to input any number of digits
      c(1.12112, MpFloat::RoundUp, 10),  // use only 10 mantisa bits and round up the initial value 1.2112
      d(1.234567890123456789),                   // only 15 digits is used to initialize (1.234... is converted to double)
      e(1.234567890123456789L);                  // here we say (L suffix) that parameter is of a long double type
  cout << "\n Precision used (in mantisa bits): " << capd::multiPrec::MpReal::getDefaultPrecision() <<"\n";
  cout.precision(30);
  cout << "\n a = " << a << "\n b = " << b << "\n c = " << c
       << "\n d = " << d << "\n e = " << e << endl;
  return 0;
}
