/////////////////////////////////////////////////////////////////////////////
/// @file SmithForm.h
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

#ifndef CAPD_FILE_SMITHFORM_H
#define CAPD_FILE_SMITHFORM_H

namespace capd
{

  namespace matrixAlgorithms
  {

    template<class MatrixT>
    struct SmithFormTraits
    {
      typedef typename MatrixT::ScalarType ScalarType;
      typedef typename MatrixT::size_type size_type;
      typedef typename MatrixT::template rebind<ScalarType>::other Matrix; // remove const!
      typedef typename MatrixT::template rebind<ScalarType>::other MatrixQ;
      typedef typename MatrixT::template rebind<ScalarType>::other MatrixR;
    };

    template<class MatrixT, typename TraitsT=SmithFormTraits<MatrixT> >
    class SmithForm
    {
    public:
      typedef TraitsT Traits;
      typedef typename Traits::Matrix Matrix;
      typedef typename Traits::MatrixQ MatrixQ;
      typedef typename Traits::MatrixR MatrixR;

      SmithForm(Matrix& B, bool computeQ, bool computeQinv, bool computeR, bool computeRinv):
        _B(B),
        _m(_B.numberOfRows()), _n(_B.numberOfColumns()),
        _Q(computeQ ? MatrixQ::Identity(_m) : MatrixQ()),
        _Qinv(computeQinv ? MatrixQ::Identity(_m) : MatrixQ()),
        _R(computeR ? MatrixR::Identity(_n) : MatrixR()),
        _Rinv(computeRinv ? MatrixR::Identity(_n) : MatrixR()),
        _s(0), _t(0)
      {
      }

      const MatrixQ& getQ() const { return _Q; }
      const MatrixQ& getQinv() const { return _Qinv; }
      const MatrixR& getR() const { return _R; }
      const MatrixR& getRinv() const { return _Rinv; }

      const int& getT() const { return _t; }
      const int& getS() const { return _s; }

      virtual void operator()() = 0;
      virtual ~SmithForm() {}

    protected:
      MatrixT& _B;
      int _m, _n;
      MatrixQ _Q, _Qinv;
      MatrixR _R, _Rinv;
      int _s, _t;
    };

  }

}

#endif // CAPD_FILE_SMITHFORM_H
