/////////////////////////////////////////////////////////////////////////////
/// @file PARISmithForm.h
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

#ifndef CAPD_FILE_PARISMITHFORM_H
#define CAPD_FILE_PARISMITHFORM_H

#include "MatrixOp.h"
#include "SmithForm.h"
#include <capd/auxil/Logger.h>

namespace capd
{
  namespace matrixAlgorithms
  {

    template<typename Scalar>
    struct GetPARISmithFormTraits
    {
      typedef capd::vectalg::Matrix<Scalar, 0, 0> Matrix;
      typedef SmithFormTraits<Matrix> type;
    };


    template<typename Scalar>
    class PARISmithForm: public SmithForm<capd::vectalg::Matrix<Scalar, 0, 0>,
                                          typename GetPARISmithFormTraits<Scalar>::type>
    {
      typedef SmithForm<capd::vectalg::Matrix<Scalar, 0, 0>,
                        typename GetPARISmithFormTraits<Scalar>::type> Base;
    public:
      typedef typename Base::Traits Traits;
      typedef typename Base::Matrix Matrix;

      PARISmithForm(Matrix& B, bool computeQ, bool computeQinv, bool computeR, bool computeRinv):
        Base(B, computeQ, computeQinv, computeR, computeRinv)
      {}

      virtual ~PARISmithForm() {}
      void operator()();

      static bool enabled();

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

    };
  }

}


#endif // CAPD_FILE_PARISMITHFORM_H
