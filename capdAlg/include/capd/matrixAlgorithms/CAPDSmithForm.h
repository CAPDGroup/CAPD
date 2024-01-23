/////////////////////////////////////////////////////////////////////////////
/// @file CAPDSmithForm.h
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-06-19
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.edu.pl/ for details.

#ifndef CAPD_FILE_CAPDSMITHFORM_H
#define CAPD_FILE_CAPDSMITHFORM_H

#include "MatrixOp.h"
#include "SmithForm.h"
#include <capd/auxil/Logger.h>

namespace capd
{
  namespace matrixAlgorithms
  {
    template<class MatrixT, typename Traits=SmithFormTraits<MatrixT> >
    class CAPDSmithForm: public SmithForm<MatrixT, Traits>
    {
      typedef SmithForm<MatrixT, Traits> Base;
    public:

      typedef typename Base::Matrix Matrix;
      typedef typename Base::MatrixQ MatrixQ;
      typedef typename Base::MatrixR MatrixR;

      CAPDSmithForm(Matrix& B, bool computeQ, bool computeQinv, bool computeR, bool computeRinv):
        Base(B, computeQ, computeQinv, computeR, computeRinv)
      {}

      virtual ~CAPDSmithForm() {}

      void operator()()
      {
        typedef typename MatrixT::ScalarType ScalarType;

        CAPD_DEBUG("SmithForm for matrix size: " << _m << "x" << _n);
        const bool trace = (_n <= 100 && _m <= 100);

        if (trace) {
          CAPD_TRACE("capdSmithForm args: " << cppReprezentation(_B, "B", "TYPE"));
        }

        _s=_t=0;
        while(nonZero(MatrixSlice<MatrixT>(_B, _t + 1, _m, _t + 1, _n))){
          _t++;
          partSmithForm(_B, _Q, _Qinv, _R, _Rinv, _t);
          if(_B(_t, _t)<ScalarType(0)) rowMultiply(_B, _Q, _Qinv, _t,-ScalarType(1));
          if(isInvertible(_B(_t, _t))) rowMultiply(_B, _Q, _Qinv, _t,inverse(ScalarType(_B(_t, _t))));
          if(_B(_t, _t)==ScalarType(1)) _s++;
        };

        if (trace) {
          CAPD_TRACE("capdSmithForm result: " << cppReprezentation(_B, "B", "TYPE"));
        }
      }

    private:
      using Base::_n;
      using Base::_m;
      using Base::_B;
      using Base::_Q;
      using Base::_Qinv;
      using Base::_R;
      using Base::_Rinv;
      using Base::_t;
      using Base::_s;

      //CAPD_CLASS_LOGGER;

      void moveMinNonzero(MatrixT& B, MatrixQ& Q,MatrixQ& Qinv,MatrixR& R,MatrixR& Rinv,int k){
        typedef typename MatrixT::ScalarType ScalarType;
        int m=B.numberOfRows();
        int n=B.numberOfColumns();
        int i,j;
        ScalarType s;
        smallestNonZero(MatrixSlice<MatrixT>(B,k,m,k,n),s,i,j);
        i+=k-1;
        j+=k-1;
        rowExchange(B,Q,Qinv,k,i);
        columnExchange(B,R,Rinv,k,j);
      }

      bool checkForDivisibility(MatrixT& B,int k,int& i,int& j,typename MatrixT::ScalarType &q){
        int m=B.numberOfRows();
        int n=B.numberOfColumns();
        for(i=k+1;i<=m;++i)
          for(j=k+1;j<=n;++j){
            q=B(i,j)/B(k,k);
            if(q*B(k,k)!=B(i,j)) return false;
          }
        return true;
      }

      void partSmithForm(MatrixT& B,MatrixQ& Q,MatrixQ& Qinv,MatrixR& R,MatrixR& Rinv,int k){
        typedef typename MatrixT::ScalarType ScalarType;
        int m=B.numberOfRows();
        int n=B.numberOfColumns();
        bool divisible=false;
        do{
          moveMinNonzero(B,Q,Qinv,R,Rinv,k);
          partRowReduce(B,Q,Qinv,k,k);
          if(nonZero(MatrixSlice<MatrixT>(B,k+1,m,k,k))) continue;
          partColumnReduce(B,R,Rinv,k,k);
          if(nonZero(MatrixSlice<MatrixT>(B,k,k,k+1,n))) continue;
          int i=0,j=0;
          ScalarType q=ScalarType(0);
          //divisible=checkForDivisibility(B,k,i,j,q);
          divisible = true;
          if(!divisible){
            rowAdd(B,Q,Qinv,i,k,ScalarType(1));
            columnAdd(B,R,Rinv,k,j,-q);
          }
        }while(!divisible);
      }

    };
  }

}

#endif // CAPD_FILE_CAPDSMITHFORM_H
