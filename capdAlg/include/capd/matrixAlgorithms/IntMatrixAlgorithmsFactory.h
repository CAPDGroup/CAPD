/////////////////////////////////////////////////////////////////////////////
/// @file IntMatrixAlgorithmsFactory.h
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-08-17
/////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library (capdAlg),
// distributed under the terms of the GNU General Public License.
// Consult http://capd.ii.uj.edu.pl and  http://redhom.ii.edu.pl/ for details.
/////////////////////////////////////////////////////////////////////////////

#ifndef CAPD_FILE_CAPDALG_MATRIXALGORITHMS_INTMATRIXALGORITHMSFACTORY_H
#define CAPD_FILE_CAPDALG_MATRIXALGORITHMS_INTMATRIXALGORITHMSFACTORY_H

// #include
#include "SmithFormFactory.h"
#include <memory> //necessary for gcc 5.2 and unique_ptr

namespace capd
{
  namespace matrixAlgorithms
  {

    struct IntMatrixAlgorithmsFactory
    {
      typedef IntMatrixAlgorithmsFactory This;

      template<typename Matrix>
      struct SmithForm
      {
        struct SmithFormThroughFactory
        {
          typedef capd::matrixAlgorithms::SmithForm<Matrix> SmithFormType;
          typedef typename SmithFormType::MatrixQ MatrixQ;
          typedef typename SmithFormType::MatrixR MatrixR;


          SmithFormThroughFactory(Matrix& B, bool computeQ, bool computeQinv, bool computeR, bool computeRinv, bool usePARI)
          {
            SmithFormFactory factory(usePARI);
            _smithForm.reset(factory(B, computeQ, computeQinv, computeR, computeRinv));
          }

          void operator()()
          {
            (*_smithForm)();
          }

          const MatrixQ& getQ() const { return _smithForm->getQ(); }
          const MatrixQ& getQinv() const { return _smithForm->getQinv(); }
          const MatrixR& getR() const { return _smithForm->getR(); }
          const MatrixR& getRinv() const { return _smithForm->getRinv(); }

          const int& getT() const { return _smithForm->getT(); }
          const int& getS() const { return _smithForm->getS(); }

        private:
          std::unique_ptr<SmithFormType> _smithForm;
        };

        typedef SmithFormThroughFactory type;
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



      IntMatrixAlgorithmsFactory(bool usePARI):
        _usePARI(usePARI)
      {}


      template<typename Matrix>
      typename SmithForm<Matrix>::type
      smithForm(Matrix& B, bool computeQ, bool computeQinv, bool computeR, bool computeRinv)
      {
        return typename SmithForm<Matrix>::type(B, computeQ, computeQinv, computeR, computeRinv, _usePARI);
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


    private:
      bool _usePARI;
    };


  }
}

#endif // CAPD_FILE_CAPDALG_MATRIXALGORITHMS_INTMATRIXALGORITHMSFACTORY_H
