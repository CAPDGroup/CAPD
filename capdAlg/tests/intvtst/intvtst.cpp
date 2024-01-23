/// @addtogroup intvtst
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file intvtst.cpp
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details. 
#include <fstream>
#include <cmath>
#include <string>
#include <cstdlib>

// To enable compilation with KRAK uncomment the following line
//#define KRAK
// Note that to  compile with KRAK you will also have to add the krak library to the makefile

#define RESTORE_ROUNDING

#include "capd/rounding/DoubleRounding.h"
#include "capd/intervals/lib.h"

using namespace std;
using capd::rounding::DoubleRounding;
using capd::interval;

void checkResult(interval & iresult, double & dresult){
  if(! iresult.contains(dresult)){
    std::cout << "ERROR: a rigorous interval result does not contain a non rigorous result"
      << "\n  interval result = " << iresult
      << "\n  double result   = " << dresult <<"\n";
    exit(1);
  }
}
void testExp(interval z)
{
 std::cout << "z=" << z << "\n";
 interval expintv=exp(z);
 #ifdef RESTORE_ROUNDING
  DoubleRounding::roundNearest();
 #endif
 double dexp=exp(z.rightBound());
 std::cout << "diff in exp: our intervalexp - standardexp " << (expintv - dexp) << "\n";
 checkResult(expintv, dexp);
}

void testExp()
{
 std::cout << "Testing exp(interval)\n\n";

 testExp(1.01);
 testExp(0.99);
 testExp(0.5);
// ponizej objawia sie blad, ale tylko przy wylaczonej obsludze niedozwolonej operacji przez
//   procesor INTEL (IM = 0 w bitach maskowania bledow numerycznych

 testExp(12.0);
}


void test_sin(interval z)
{
 std::cout << "z=" << z << "\n";
 interval sinintv=sin(z);
 #ifdef RESTORE_ROUNDING
  DoubleRounding::roundNearest();
 #endif 
 double dsin=sin(z.rightBound());
 interval diff = sinintv - dsin;
 std::cout << "diff in sin: our interval sin  - standard sin " << diff << "\n";
 checkResult(sinintv, dsin);
}

void test_sin()
{
 std::cout << "Testing sin(interval)\n\n";
 test_sin(interval::pi());
 test_sin(6.0);
 test_sin(1.0);
}

void show_rounding_test(void)
{
 string result;
 std::cout << "\nTesting rounding:\n\n";

 DoubleRounding::roundNearest();
 result=(DoubleRounding::test() == 0 ? "OK" : "BAD");
 std::cout << "Rounding nearest: " << result << "\n";

 DoubleRounding::roundDown();
 result=(DoubleRounding::test() == 1 ? "OK" : "BAD");
 std::cout << "Rounding down:    " << result << "\n";

 DoubleRounding::roundUp();
 result=(DoubleRounding::test() == 2 ? "OK" : "BAD");
 std::cout << "Rounding up:      " << result << "\n";

 DoubleRounding::roundCut();
 result=(DoubleRounding::test() == 3 ? "OK" : "BAD");
 std::cout << "Rounding cut:     " << result << "\n";

 std::cout << "isWorking:        " << DoubleRounding::isWorking() << "\n";
}

int main(int, char*[])
{
 //init_fpunit();
 
// test verifying if rounding mechanism works properly
 show_rounding_test();

#ifdef UP_ROUNDING
 DoubleRounding::roundUp();
#endif

// test for the exp(interval) function - it has problems!!!
 testExp();
    
// simple test of sin function
 test_sin();

 return 0;
}

/// @}
