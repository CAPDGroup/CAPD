/// @addtogroup normalForms
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file normalForms.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_NORMALFORMS_NORMALFORMS_H_
#define _CAPD_NORMALFORMS_NORMALFORMS_H_

#include <complex>
#include <stdexcept>
#include <sstream>

#include "capd/diffAlgebra/Jet.h"
#include "capd/vectalg/Vector.hpp"
#include "capd/vectalg/Matrix.hpp"
#include "capd/vectalg/Multiindex.h"

namespace capd{
namespace normalForms{

// -------------------------------------------------------------------------- //

template<typename MatrixType, unsigned DEGREE>
void derivativesToSeries(capd::diffAlgebra::Jet< MatrixType,DEGREE >& s)
{
  typedef typename capd::diffAlgebra::Jet<MatrixType,DEGREE>::ScalarType ScalarType;
  for(typename MatrixType::size_type i=2;i<=s.degree();++i)
  {
    typename capd::diffAlgebra::Jet<MatrixType,DEGREE>::Multipointer mp = s.first(i);
    do{
     // cannot be just double - fails with Scalar=std::complex<Interval<MpFloat>>
     typename capd::TypeTraits<ScalarType>::Real fac = mp.factorial();
     s(mp) /= ScalarType(fac);
    }while(s.hasNext(mp));
  }
}

template<typename MatrixType, unsigned DEGREE>
void seriesToDerivatives(capd::diffAlgebra::Jet<MatrixType,DEGREE>& s)
{
  typedef typename capd::diffAlgebra::Jet<MatrixType,DEGREE>::ScalarType ScalarType;
  for(typename MatrixType::size_type i=2;i<=s.degree();++i)
  {
    typename capd::diffAlgebra::Jet<MatrixType,DEGREE>::Multipointer mp = s.first(i);
    do{
     typename capd::TypeTraits<ScalarType>::Real fac = mp.factorial();
     s(mp) *= ScalarType(fac);
    }while(s.hasNext(mp));
  }
}

}} // namespace capd::normalForms

#include "capd/normalForms/linearSubstitution.hpp"
#include "capd/normalForms/planarMaps.hpp"

#endif // _CAPD_NORMALFORMS_NORMALFORMS_H_

/// @}
