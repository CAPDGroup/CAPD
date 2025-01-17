//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file intervals/lib.h
///
/// @author Tomasz Kapela   @date 2010-01-23
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) Tomasz Kapela 2010
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef CAPD_INTERVAL_LIB_H
#define CAPD_INTERVAL_LIB_H

#include "minmax_interval.h"

#include <cmath>
#include "capd/basicalg/doubleFun.h"

#if defined(__USE_CXSC__)

#include "capd/cxsc/Interval.h"
namespace capd
{

typedef ::capd::cxsc::Interval  DInterval;

}

#elif defined(__USE_NATIVE__)

using std::log;
#include "capd/intervals/Interval.h"
#include "capd/rounding/DoubleRounding.h"

namespace capd
{

typedef intervals::Interval<double, capd::rounding::DoubleRounding> DInterval;

}

#elif defined(__USE_FILIB__)

#include "capd/filib/Interval.h"
namespace capd
{

typedef ::capd::filib::Interval<double, ::filib::native_directed, ::filib::i_mode_normal >  DInterval;

}

#else

#error "No interval type selected!"

#endif

#ifndef __CAPD_DEFINE_INTERVAL__
#define __CAPD_DEFINE_INTERVAL__
namespace capd{
  typedef DInterval interval;
  typedef DInterval Interval;

}

#endif //__CAPD_DEFINE_INTERVAL__

#endif // CAPD_INTERVAL_LIB_H
