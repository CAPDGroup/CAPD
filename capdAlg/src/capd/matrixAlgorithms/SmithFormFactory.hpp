/////////////////////////////////////////////////////////////////////////////
/// @file SmithFormFactory.hpp
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-06-20
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.edu.pl/ for details.

#ifndef CAPD_FILE_SMITHFORMFACTORY_HPP
#define CAPD_FILE_SMITHFORMFACTORY_HPP

#include <capd/matrixAlgorithms/SmithFormFactory.h>
#include <capd/matrixAlgorithms/CAPDSmithForm.h>
#include <capd/matrixAlgorithms/PARISmithForm.h>
#include <capd/vectalg/Matrix.hpp>

using capd::vectalg::Matrix;

namespace capd
{
  namespace matrixAlgorithms
  {

    template<typename MatrixT>
    struct _CreateSmithForm_
    {
    };

    template<typename Scalar>
    struct _CreateSmithForm_<Matrix<Scalar, 0, 0> >
    {
      typedef Matrix<Scalar, 0, 0> MatrixT;

      SmithForm<MatrixT>* operator()(MatrixT& B, bool computeQ, bool computeQinv, bool computeR, bool computeRinv, bool /*usePari*/)
      {
        return new CAPDSmithForm<MatrixT>(B, computeQ, computeQinv, computeR, computeRinv);
      }
    };

#define CREATE_SMITH_FORM(type) \
    template<> \
    struct _CreateSmithForm_<Matrix<type, 0, 0> > \
    {\
      typedef Matrix<type, 0, 0> MatrixT;\
      SmithForm<MatrixT>* operator()(MatrixT& B, bool computeQ, bool computeQinv, bool computeR, bool computeRinv, bool usePari) \
      {                                                                 \
        if (usePari && PARISmithForm<type>::enabled()) {                \
          return new PARISmithForm<type>(B, computeQ, computeQinv, computeR, computeRinv); \
        } else {                                                        \
          return new CAPDSmithForm<MatrixT>(B, computeQ, computeQinv, computeR, computeRinv); \
        }                                                               \
      }                                                                 \
    };

    CREATE_SMITH_FORM(short);
    CREATE_SMITH_FORM(int);
    CREATE_SMITH_FORM(long);
    typedef long long llong;
    CREATE_SMITH_FORM(llong);

    template<typename MatrixT>
    SmithForm<MatrixT>* SmithFormFactory::operator()(MatrixT& B, bool computeQ, bool computeQinv, bool computeR, bool computeRinv)
    {
      return _CreateSmithForm_<MatrixT>()(B, computeQ, computeQinv, computeR, computeRinv, _usePARI);
    }

  }
}

#endif // CAPD_FILE_SMITHFORMFACTORY_HPP
