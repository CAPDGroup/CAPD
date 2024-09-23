//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file fields/lib.h
///
/// @author Tomasz Kapela   @date 2024-09-23
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) Tomasz Kapela 2024
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_FIELDS_LIB_H_
#define _CAPD_FIELDS_LIB_H_

#include "capd/fields/Complex.h"
#include "capd/intervals/lib.h"
namespace capd{
  typedef ::capd::fields::Complex<double> Complex;
  typedef ::capd::fields::Complex<interval> IComplex;
}


#endif // CAPD_FIELDS_LIB_H_
