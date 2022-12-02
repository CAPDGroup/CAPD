/////////////////////////////////////////////////////////////////////////////
/// @file PARISmithForm.cpp
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-06-20
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD::RedHom Group.
//
// This file constitutes a part of the CAPD::RedHom library,
// distributed under the terms of the GNU General Public License.
// Consult  http://redhom.ii.edu.pl/ for details.

#include "PARISmithForm.hpp"
#include <capd/vectalg/Matrix.hpp>
#include <capd/vectalg/Vector.hpp>

namespace capd
{
  namespace matrixAlgorithms
  {

    PARI_SMITH_FORM(short);
    PARI_SMITH_FORM(int);
    PARI_SMITH_FORM(long);
    typedef long long llong;
    PARI_SMITH_FORM(llong);
  }
}
