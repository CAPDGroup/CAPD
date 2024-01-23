//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file IComplex.h
/// @deprecated
/// @author Tomasz Kapela   @date 2010-03-15
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) Tomasz Kapela 2010
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_INTERVAL_ICOMPLEX_H_
#define _CAPD_INTERVAL_ICOMPLEX_H_

#include "capd/fields/Complex.h"

namespace capd{
namespace intervals {

/// Definition for backward compatibility
/// @deprecated
template<typename T>
using IComplex = capd::fields::Complex<T>;
}}

#endif // _CAPD_INTERVAL_ICOMPLEX_H_
