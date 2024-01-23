/////////////////////////////////////////////////////////////////////////////
/// @file SmithFormFactory.cpp
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-06-19
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD::RedHom Group.
//
// This file constitutes a part of the CAPD::RedHom library,
// distributed under the terms of the GNU General Public License.
// Consult  http://redhom.ii.edu.pl/ for details.

#include "SmithFormFactory.hpp"

#include <capd/matrixAlgorithms/SmithFormFactory.h>
#include <capd/matrixAlgorithms/CAPDSmithForm.h>
#include <capd/vectalg/Matrix.hpp>
#include <capd/rings/Zp.h>
#include <capd/rings/Z2.h>

using capd::vectalg::Matrix;

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


    INSTANCE_DEF(short);
    INSTANCE_DEF(int);
    INSTANCE_DEF(long);
    typedef long long llong;
    INSTANCE_DEF(llong);
    INSTANCE_DEF(Zp);
    INSTANCE_DEF(Z2);
  }
}
