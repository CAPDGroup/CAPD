/// @file tightSetupTest.cpp
///
/// TIGHT intervals require the FPU in the round up mode all the time.
/// The mode is set automatically at program start by an inline variable defined in capd/tight/Interval.h.
/// This test checks (without any explicit call of Interval::setup()) that
///  - the FPU is in the round up mode when main() starts,
///  - the mode is kept by TIGHT interval operations,
///  - decimal strings are parsed to rigorous enclosures (which requires round up).

#include <cfenv>
#include <iostream>
#include <sstream>
#include "capd/intervals/lib.h"

// remembered before main() is entered
static const int roundingAtStatic = std::fegetround();

int main(){
  int errors = 0;
  if(std::fegetround() != FE_UPWARD){
    std::cerr << "FPU is not in the round up mode at the start of main()\n";
    ++errors;
  }
  // a static object of another translation unit could be initialized earlier; this one is after the header
  if(roundingAtStatic != FE_UPWARD){
    std::cerr << "FPU is not in the round up mode during static initialization\n";
    ++errors;
  }

  capd::DInterval x(1.0), y(3.0);
  capd::DInterval z = x / y;
  if(std::fegetround() != FE_UPWARD){
    std::cerr << "FPU left the round up mode after an interval operation\n";
    ++errors;
  }
  if(!(z.leftBound() < z.rightBound())){
    std::cerr << "1/3 is not enclosed in a nondegenerate interval: " << z << "\n";
    ++errors;
  }

  capd::DInterval a;
  std::istringstream s("[3.21312312, 4.324324324]");
  s >> a;
  if(!(a.rightBound() > 4.324324324)){
    std::cerr << "decimal string was not parsed to a rigorous enclosure\n";
    ++errors;
  }
  return errors;
}
