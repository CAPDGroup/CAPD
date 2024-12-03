/// @addtogroup map
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file mapLib.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details. 

#ifndef _CAPD_MAP_MAPLIB_H_ 
#define _CAPD_MAP_MAPLIB_H_ 

#include "capd/vectalg/vectalgLib.h"
#include "capd/map/Function.h"
#include "capd/map/Map.h"

typedef capd::map::Function<IVector> IFunction;
typedef capd::map::Map<IMatrix> IMap;

typedef capd::map::Function<DVector> DFunction;
typedef capd::map::Map<DMatrix> DMap;

typedef capd::map::Function<LDVector> LDFunction;
typedef capd::map::Map<LDMatrix> LDMap;

typedef capd::map::Function<IVectorMD> IFunctionMD;
typedef capd::map::Map<IMatrixMD> IMapMD;

typedef capd::map::Function<DVectorMD> DFunctionMD;
typedef capd::map::Map<DMatrixMD> DMapMD;

typedef capd::map::Function<LDVectorMD> LDFunctionMD;
typedef capd::map::Map<LDMatrixMD> LDMapMD;

#if (CAPD_CPU_ARCH==CAPD_CPU_ARCH_X86)

  typedef capd::map::Function<LIVector> LIFunction;
  typedef capd::map::Map<LIMatrix> LIMap;

  typedef capd::map::Function<LIVectorMD> LIFunctionMD;
  typedef capd::map::Map<LIMatrixMD> LIMapMD;

#endif



#endif // _CAPD_MAP_MAPLIB_H_ 

/// @}
