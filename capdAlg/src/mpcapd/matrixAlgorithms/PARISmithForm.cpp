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


#include "../../capd/matrixAlgorithms/PARISmithForm.hpp"
#include <capd/multiPrec/MpInt.h>
#include <capd/vectalg/Matrix.hpp>
#include <capd/vectalg/Vector.hpp>

namespace capd
{
  namespace matrixAlgorithms
  {

    using capd::multiPrec::MpInt;
    PARI_SMITH_FORM(MpInt);
  }
}
