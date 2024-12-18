/////////////////////////////////////////////////////////////////////////////
/// @file LongDoubleInterval.h
///
/// @author Tomasz Kapela   @date 11-01-2006
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) Tomasz Kapela 2006
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details. 

#ifndef CAPD_INTERVAL_LONGDOUBLEINTERVAL_H
#define CAPD_INTERVAL_LONGDOUBLEINTERVAL_H

#include <cmath>

using std::log;

#include "capd/basicalg/doubleFun.h"
#include "capd/intervals/Interval.h"
#include "capd/rounding/DoubleRounding.h"

#if (CAPD_CPU_ARCH==CAPD_CPU_ARCH_X86)

typedef capd::intervals::Interval<long double, capd::rounding::DoubleRounding> LInterval;

//using namespace capd::intervals;
namespace capd{
	namespace intervals{
		typedef capd::intervals::Interval<long double, capd::rounding::DoubleRounding> LongDoubleInterval;
	}
}


#endif // check of architecture

#endif // CAPD_INTERVAL_LONGDOUBLEINTERVAL_H
