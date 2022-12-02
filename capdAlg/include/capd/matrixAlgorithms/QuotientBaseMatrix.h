/////////////////////////////////////////////////////////////////////////////
/// @file QuotientBaseMatrix.h
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

#ifndef CAPD_FILE_CAPDALG_MATRIXALGORITHMS_QUOTIENTBASEMATRIX_H
#define CAPD_FILE_CAPDALG_MATRIXALGORITHMS_QUOTIENTBASEMATRIX_H

// #include
#include <capd/auxil/Logger.h>
#include <vector>
#include <capd/auxil/RemoveConst.h>


namespace capd
{
  namespace matrixAlgorithms
  {

    template<class Matrix, typename IntMatrixAlgorithms>
    class QuotientBaseMatrix
    {

    public:
      typedef typename Matrix::ScalarType ScalarType;
      typedef typename auxil::RemoveConst<Matrix>::type InternalMatrix;
      typedef typename IntMatrixAlgorithms::template SmithForm<InternalMatrix>::type SmithForm;
      typedef std::vector<ScalarType> Orders;

      QuotientBaseMatrix(Matrix& A_W, Matrix& A_V, IntMatrixAlgorithms& intMatrixAlgorithms):
        A_W(A_W), A_V(A_V),
        _intMatrixAlgorithms(intMatrixAlgorithms)
      {}

      virtual ~QuotientBaseMatrix() {}

      void operator()()
      {
        typedef typename Matrix::ColumnVectorType VectorType;
        const int p=A_W.numberOfRows();
        const int m=A_W.numberOfColumns();
        const int n=A_V.numberOfColumns();
        const bool trace = (p < 100 && m < 100 && n < 100);

        CAPD_DEBUG("quotientBaseMatrix with sizes: " << p << " " << m << " " << n);

        if (trace) {
          CAPD_TRACE("quotientBaseMatrix args: " << cppReprezentation(A_W, "A_W", "TYPE") << " "
                     << cppReprezentation(A_V, "A_V", "TYPE"));
        }

        VectorType a;
        InternalMatrix B(m,n);
        typedef typename IntMatrixAlgorithms::template SolveLinearEquation<Matrix>::type Solver;
        Solver solver = _intMatrixAlgorithms.solveLinearEquation(A_W);

        for(int j=0;j<n;++j){
          solver(A_V.column(j),a);
          for(int i=0;i<m;++i) B[i][j]=a[i];
        }

        SmithForm smithForm = _intMatrixAlgorithms.smithForm(B, true, false, false, false);
        smithForm();

        const int s = smithForm.getS();
        _orders.resize(n - s);
        for(int i=s;i<n;++i) _orders[i-s]=B[i][i];

        InternalMatrix WQ=A_W * smithForm.getQ();
        A_U=MatrixSlice<InternalMatrix>(WQ,1,p,s+1,m);
      }
      const Matrix& pseudoBasis() const { return A_U; }

      template<typename V>
      void getOrders(V& v) const
      {
        v.resize(_orders.size());
        std::copy(_orders.begin(), _orders.end(), v.begin());
      }

    private:
      Matrix& A_W;
      Matrix& A_V;
      IntMatrixAlgorithms& _intMatrixAlgorithms;
      InternalMatrix A_U;
      Orders _orders;
    };


  }
}

#endif // CAPD_FILE_CAPDALG_MATRIXALGORITHMS_QUOTIENTBASEMATRIX_H
