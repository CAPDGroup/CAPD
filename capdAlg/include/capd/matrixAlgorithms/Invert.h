/////////////////////////////////////////////////////////////////////////////
/// @file Invert.h
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

#ifndef CAPD_FILE_INVERT_H
#define CAPD_FILE_INVERT_H

#include "MatrixOp.h"
#include "CAPDSmithForm.h"


namespace capd
{

  namespace matrixAlgorithms
  {
    template<typename matrix>
    struct IInvert
    {
      virtual bool operator()(matrix& result) const = 0;

      virtual matrix operator()() const
      {
	matrix tmp;
	if ((*this)(tmp)) {
	  return tmp;
	} else {
	  throw std::invalid_argument("Cannot invert the matrix");
	}
      }
    };

    template<typename matrix, typename Scalar=typename matrix::ScalarType>
    class Invert: public IInvert<matrix>
    {

    public:
      explicit Invert(const matrix& A):
	_m(A.numberOfRows()), _n(A.numberOfColumns()),
	_B(A), _Qinv(_m, _m), _R(_n, _n)
      {
	CAPDSmithForm<matrix> smithForm(_B, false, true, true, false);
	smithForm();

	_Qinv = smithForm.getQinv();
	_R = smithForm.getR();
      }

      bool operator()(matrix& Ainv) const
      {
	typedef typename matrix::ScalarType ScalarType;
	if(_m != _n) return false;

	for(int i=1; i <= _n; ++i) {
	  if(_B(i,i) != ScalarType(1)) {
	    return false;
	  }
	}
	Ainv = _R * _Qinv;

	return true;
      }

    private:
      int _m, _n;
      matrix _B;
      matrix _Qinv, _R;
    };

  }

}

#endif // CAPD_FILE_INVERT_H
