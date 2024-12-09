/// @addtogroup diffIncl
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file DiffInclusion.hpp
///
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

/* Author: Tomasz Kapela, 2007 */

#ifndef CAPD_DIFFINCL_DIFFINCLUSION_HPP
#define CAPD_DIFFINCL_DIFFINCLUSION_HPP

#include <sstream>
#include <string>
#include <stdexcept>

#include "capd/dynsys/OdeSolver.hpp"
#include "capd/diffIncl/DiffInclusion.h"

namespace capd{
namespace diffIncl{


template <typename MapT, typename DynSysT>
DiffInclusion<MapT, DynSysT>::DiffInclusion(
           MultiMapType& diffIncl, 
           int order, 
           const NormType & norm
) 
  : m_norm ( norm.clone()),
    m_diffIncl(diffIncl),
    m_dynamicalSystem(diffIncl.getVectorField(), order) {
}

template <typename MapT, typename DynSysT>
DiffInclusion<MapT, DynSysT>::~DiffInclusion(){
  delete m_norm;
}
  

  
}} //namespace capd::diffIncl`

#endif // CAPD_DIFFINCL_DIFFINCLUSION_HPP

/// @}
