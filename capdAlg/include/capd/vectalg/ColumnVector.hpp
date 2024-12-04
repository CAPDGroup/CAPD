/// @addtogroup vectalg
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file ColumnVector.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details.

#ifndef _CAPD_VECTALG_COLUMNVECTOR_HPP_
#define _CAPD_VECTALG_COLUMNVECTOR_HPP_

#include <stdexcept>
#include <sstream>
#include <algorithm>

#include "capd/vectalg/ColumnVector.h"
#include "capd/vectalg/Vector.hpp"

namespace capd{
namespace vectalg{

template<typename Scalar,__size_type rows>
ColumnVector<Scalar,rows>::operator Vector<Scalar,rows>() const
{
  Vector<Scalar,rows> result(dimension(),true);
  const_iterator b=begin(), e=end();
  typename Vector<Scalar,rows>::iterator i=result.begin();
  while(b!=e)
  {
    *i = *b;
    ++i;
    ++b;
  }
  return result;
}


template<typename Scalar, __size_type rows>
std::string cppReprezentation(const ColumnVector<Scalar,rows> & A, const std::string& varName,
			      const std::string& typeName)
{
  std::stringstream out;
  out << "capd::vectalg::ColumnVector<" << typeName << ", " << rows << "> " << varName << "(";

  if (A.dimension() > 0) {
    out << "(" << typeName << "[" << A.dimension() << "])" << A;
  } else {
    out << "(__size_type)" << A.dimension();
  }
  out << ");";

  std::string str = out.str();
  std::replace(str.begin(), str.end(), '\n', ' ');

  return str;
}

}} // namespace capd::vectalg

#endif // _CAPD_VECTALG_COLUMNVECTOR_HPP_

/// @}
