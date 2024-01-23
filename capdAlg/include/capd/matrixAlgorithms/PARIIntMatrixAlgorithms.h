/////////////////////////////////////////////////////////////////////////////
/// @file PARIIntMatrixAlgorithms.h
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-08-09
/////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library (capdAlg),
// distributed under the terms of the GNU General Public License.
// Consult http://capd.ii.uj.edu.pl and  http://redhom.ii.edu.pl/ for details.
/////////////////////////////////////////////////////////////////////////////

#ifndef CAPD_FILE_CAPDALG_MATRIXALGORITHMS_PARIINTMATRIXALGORITHMS_H
#define CAPD_FILE_CAPDALG_MATRIXALGORITHMS_PARIINTMATRIXALGORITHMS_H

// #include
#include "PARISmithForm.h"
#include "SolveLinearEquation.h"
#include "QuotientBaseMatrix.h"

namespace capd
{
  namespace matrixAlgorithms
  {

    struct PARIIntMatrixAlgorithms
    {
      typedef PARIIntMatrixAlgorithms This;

      template<typename Matrix>
      struct SmithForm
      {
        typedef matrixAlgorithms::PARISmithForm<typename Matrix::ScalarType> type;
      };

      template<typename Matrix>
      struct SolveLinearEquation
      {
        typedef matrixAlgorithms::SolveLinearEquation<Matrix, This> type;
      };

      template<typename Matrix>
      struct QuotientBaseMatrix
      {
        typedef matrixAlgorithms::QuotientBaseMatrix<Matrix, This> type;
      };


      template<typename Matrix>
      typename SmithForm<Matrix>::type
      smithForm(Matrix& B, bool computeQ, bool computeQinv, bool computeR, bool computeRinv)
      {
        typename SmithForm<Matrix>::type(B, computeQ, computeQinv, computeR, computeRinv);
      }

      template<typename Matrix>
      typename SolveLinearEquation<Matrix>::type
      solveLinearEquation(Matrix& m)
      {
        return typename SolveLinearEquation<Matrix>::type(m, *this);
      }

      template<typename Matrix>
      typename QuotientBaseMatrix<Matrix>::type
      quotientBaseMatrix(Matrix& A_W, Matrix& A_V)
      {
        return typename QuotientBaseMatrix<Matrix>::type(A_W, A_V, *this);
      }

    };

  }
}

#endif // CAPD_FILE_CAPDALG_MATRIXALGORITHMS_PARIINTMATRIXALGORITHMS_H
