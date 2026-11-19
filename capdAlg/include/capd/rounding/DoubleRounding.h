
/////////////////////////////////////////////////////////////////////////////
/// @file DoubleRounding.h
///
/// @author CAPD group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2006 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details. 

#ifndef CAPD_ROUNDING_DOUBLEROUNDING_H
#define CAPD_ROUNDING_DOUBLEROUNDING_H

// This allows automatic detection of operating system and compiler
#include <cfenv>
#include "capd/settings/archSetting.h"
//#include "capd/operatingSystemSetting.h"
//#include "capd/compilerSetting.h"

namespace capd{
namespace rounding{
/// @addtogroup rounding
/// @{

///  Rounding modes
enum RoundingMode { RoundUnknown= -1, RoundNearest, RoundDown, RoundUp,  RoundCut };
 
//////////////////////////////////////////////////////////////////////////////
//   DoubleRounding
///
///  Definition of class that switches rounding modes of double numbers
///
///   @author Tomasz Kapela   @date 11-01-2006
//////////////////////////////////////////////////////////////////////////////
class DoubleRounding{

    void initFpUnit();           ///<  Initialize Floating Point Unit

   public:
    DoubleRounding();            ///<  Call initialization of Floating Point Unit

    static void roundNearest();  ///<  Sets rounding to nearest mode
    
    static void roundUp();       ///<  Sets rounding up mode
    
    static void roundDown();     ///<  Sets rounding down mode
    
    static void roundCut();      ///<  Sets rounding towards zero mode

    /// Sets the rounding mode in which the selected interval library expects the FPU between its operations:
    /// round up for TIGHT (__USE_TIGHT__, it relies on the FPU being always in this mode), round to nearest otherwise.
    /// Code that needs only an approximate double arithmetic should call this function instead of roundNearest().
    static void roundDefault();

    static RoundingMode test();  ///<  Returns the actual rounding mode
    
    static bool isWorking();  ///<  Tests that all rounding modes works correctly. Returns false is they don't
};

//////////////////////////////////////////////////////////////////////////////
//   ScopedRounding
///
///  RAII guard for the rounding mode of double numbers. The constructor remembers the current
///  rounding mode of the FPU and switches to the requested one (round to nearest by default),
///  the destructor restores the remembered mode - also when the scope is left by 'continue',
///  'return' or an exception. Use it around plain double computations (std::exp, std::log, file I/O, ...)
///  that must be done with round to nearest, in code that may be called in any rounding mode
///  (e.g. round up, which is required by TIGHT intervals between their operations):
///
///      { capd::rounding::ScopedRounding nearest;   // round to nearest here
///        double step = std::exp(std::log(x)/n); }  // the previous mode is back here
//////////////////////////////////////////////////////////////////////////////
class ScopedRounding{
  public:
    explicit ScopedRounding(RoundingMode mode = RoundNearest) : m_saved(current()) {
      set(mode);
    }
    ~ScopedRounding(){
      set(m_saved);
    }
    ScopedRounding(const ScopedRounding&) = delete;
    ScopedRounding& operator=(const ScopedRounding&) = delete;

    /// the rounding mode currently set in the FPU (RoundUnknown if it is none of the four modes)
    static RoundingMode current(){
      switch(std::fegetround()){
        case FE_TONEAREST:  return RoundNearest;
        case FE_UPWARD:     return RoundUp;
        case FE_DOWNWARD:   return RoundDown;
        case FE_TOWARDZERO: return RoundCut;
        default:            return RoundUnknown;
      }
    }
    /// sets the given rounding mode (RoundUnknown leaves the FPU unchanged)
    static void set(RoundingMode mode){
      switch(mode){
        case RoundNearest: DoubleRounding::roundNearest(); break;
        case RoundUp:      DoubleRounding::roundUp();      break;
        case RoundDown:    DoubleRounding::roundDown();    break;
        case RoundCut:     DoubleRounding::roundCut();     break;
        default: ;
      }
    }
  private:
    RoundingMode m_saved;
};

#if CAPD_CPU_ARCH==CAPD_CPU_ARCH_SPARC
/* on SPARC rounding_* are inline wrappers around standard functions
 * on other architectures the functions are implemented in asm and shouldn't
 * be inline for the compiler not to reorder the instructions */

#include <ieeefp.h>

inline void DoubleRounding::roundNearest(){
  fpsetround(FP_RN);
}
inline void DoubleRounding::roundUp(){
  fpsetround(FP_RP);
}
inline void DoubleRounding::roundDown(){
  fpsetround(FP_RM);
}
inline void DoubleRounding::roundCut(){
  fpsetround(FP_RZ);
}

#endif

/// @}
}} // namespace capd::rounding

extern capd::rounding::DoubleRounding initRounding;

#endif // CAPD_ROUNDING_DOUBLEROUNDING_H

