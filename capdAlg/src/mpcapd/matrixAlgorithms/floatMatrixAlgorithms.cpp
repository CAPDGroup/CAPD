
/////////////////////////////////////////////////////////////////////////////
/// @file mpcapd/matrixAlgorithms/floatMatrixAlgorithms.cpp
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include "capd/basicalg/minmax.h"
#include "capd/vectalg/mplib.h"
#include "capd/matrixAlgorithms/floatMatrixAlgorithms.hpp"
#include "capd/vectalg/Matrix_Interval.hpp"

#ifdef __HAVE_MPFR__
  using namespace capd;
  template void capd::matrixAlgorithms::orthonormalize<MpMatrix>(MpMatrix&);
  template void capd::matrixAlgorithms::QR_decompose<MpMatrix>(const MpMatrix&,MpMatrix&,MpMatrix&);
  template int capd::matrixAlgorithms::symMatrixDiagonalize<MpMatrix>(const MpMatrix&,MpMatrix&,MpFloat);
  template MpFloat capd::matrixAlgorithms::spectralRadiusOfSymMatrix<MpMatrix>(const MpMatrix&,MpFloat);
  template MpFloat capd::matrixAlgorithms::maxEigenValueOfSymMatrix<MpMatrix>(const MpMatrix&,MpFloat);
  template MpVector capd::matrixAlgorithms::gauss<MpMatrix,MpVector>(const MpMatrix&,const MpVector&);
  template void capd::matrixAlgorithms::croutDecomposition<MpMatrix>(const MpMatrix&,MpMatrix&,MpMatrix&);
  template MpMatrix capd::matrixAlgorithms::invLowerTriangleMatrix<MpMatrix>(const MpMatrix&);
  template MpMatrix capd::matrixAlgorithms::invUpperTriangleMatrix<MpMatrix>(const MpMatrix&);
  template MpMatrix capd::matrixAlgorithms::inverseMatrix<MpMatrix>(const MpMatrix&);
  template MpMatrix capd::matrixAlgorithms::gaussInverseMatrix<MpMatrix>(const MpMatrix&);
  template MpMatrix capd::matrixAlgorithms::matrixExp<MpMatrix>(const MpMatrix&,MpFloat);

  template void capd::matrixAlgorithms::orthonormalize<MpIMatrix>(MpIMatrix&);
  template void capd::matrixAlgorithms::orthonormalize<MpIMatrix>(MpIMatrix&,const MpIVector&);
  template void capd::matrixAlgorithms::QR_decompose<MpIMatrix>(const capd::MpIMatrix&,MpIMatrix&,MpIMatrix&);
  template int capd::matrixAlgorithms::symMatrixDiagonalize<MpIMatrix>(const capd::MpIMatrix&,MpIMatrix&,MpInterval);
  template MpInterval capd::matrixAlgorithms::spectralRadiusOfSymMatrix<MpIMatrix>(const capd::MpIMatrix&,MpInterval);
  template MpInterval capd::matrixAlgorithms:: maxEigenValueOfSymMatrix<MpIMatrix>(const capd::MpIMatrix&,MpInterval);
  template MpIVector capd::matrixAlgorithms::gauss<MpIMatrix,MpIVector>(const MpIMatrix&, const MpIVector&);
  template void capd::matrixAlgorithms::croutDecomposition<MpIMatrix>(const capd::MpIMatrix&,MpIMatrix&,MpIMatrix&);
  template capd::MpIMatrix capd::matrixAlgorithms::invLowerTriangleMatrix<MpIMatrix>(const capd::MpIMatrix&);
  template capd::MpIMatrix capd::matrixAlgorithms::invUpperTriangleMatrix<MpIMatrix>(const capd::MpIMatrix&);
  template capd::MpIMatrix capd::matrixAlgorithms::inverseMatrix<MpIMatrix>(const capd::MpIMatrix&);
  template capd::MpIMatrix capd::matrixAlgorithms::gaussInverseMatrix<MpIMatrix>(const capd::MpIMatrix&);
  template capd::MpIMatrix capd::matrixAlgorithms::krawczykInverse<MpIMatrix>(const capd::MpIMatrix&);
  template capd::MpIMatrix capd::matrixAlgorithms::matrixExp<MpIMatrix>(const MpIMatrix&,MpInterval);


#endif
