
/////////////////////////////////////////////////////////////////////////////
/// @file mpExample2.cpp
/// Example how to use intervals with multiple precision endpoints
///
/// @author kapela  @date 2010-01-07
///
// ///////////////////////////////////////////////////////////////////////////

// Copyright (C) CAPD group
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include <iostream>
using namespace std;

// symbol __HAVE_MPFR__ is defined if you compile with multiple-precision support (e.g. make target=X11-gmp)
#ifdef __HAVE_MPFR__

// definition of an interval with multiple-precision end points
#include "capd/intervals/MpInterval.h"
#include "capd/multiPrec/MpReal.h"
//typedef capd::multiPrec::MpReal          MpFloat;        // multiple-precision floating point number
//typedef capd::intervals::Interval<MpFloat>  MpInterval;     // interval with MpFloat endpoints
//typedef capd::multiPrec::MpReal          MpReal;
using namespace capd::intervals;
int main(){
  
  
  MpReal::setDefaultPrecision(200);
  MpInterval sum = capd::TypeTraits<MpInterval>::zero();
  
  MpInterval  x("-1e-50","1e50");
  cout.precision(60);
  for(int i=0; i<100; i++){
    sum += x;
  }
  cout << "x = " << x << "\nsum = " <<sum<<endl;
  MpInterval expected("-1e-48","1e-48");
  if(sum.contains(expected))
    return 0;
  else
    return 1;
}

#else
int main(){
  cout << "This program was compiled without multiple precision support."
       << "\nThe __HAVE_MPFR__ symbol was not defined during compilation "
       << "\n(most probably during CAPD compilation mpfr library was not found)\n    ";
  return 0;
}
#endif

