/////////////////////////////////////////////////////////////////////////////
/// @file SmithFormFactory.cpp
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

#include <capd/multiPrec/MpInt.h>
#include "../../capd/matrixAlgorithms/SmithFormFactory.hpp"

#include <capd/vectalg/Matrix.hpp>
#include <capd/vectalg/Vector.hpp>

using capd::vectalg::Matrix;
using capd::multiPrec::MpInt;

namespace capd
{
  namespace matrixAlgorithms
  {

#define INSTANCE_MAT(type) \
    template \
    SmithForm<type>* SmithFormFactory::operator()<type>(type& B, \
    bool computeQ, bool computeQinv, bool computeR, bool computeRinv);

#define INSTANCE_DEF(type) \
    typedef Matrix<type,0,0> XXXX##type; \
    INSTANCE_MAT(XXXX##type)


    INSTANCE_DEF(MpInt);
  }
}
