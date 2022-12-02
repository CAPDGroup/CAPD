//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file divtst.cpp
///
///
/// @author kapela  @date 2010-02-23
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) CAPD group 
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details. 

#include <iostream>
#include "capd/intervals/lib.h"
using namespace capd;
int main ()
{
    interval x3 = interval(1.0, 1.0) / 3.0;
    std::cout << "[1/3] = ";
    hexWrite(std::cout, x3) << std::endl;
    if (x3. leftBound () != x3. rightBound ()){
        std::cout << "Result seems to be OK\n";
    } else {
        std::cout << "Incorrect result in division (1.0,1.0)/3.0 for native CAPD intervals.\n";
#ifndef __USE_FILIB__
        return 1;
#else
        std::cout << "Therefore CAPD library is using filib intervals.\n\n";
#endif
    }
    return 0;
}


