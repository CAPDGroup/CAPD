/// @addtogroup dynsys
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file dynsysLib.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details. 

/* Author: Daniel Wilczak 2005 */

#ifndef _CAPD_DYNSYS_DYNSYSLIB_H_ 
#define _CAPD_DYNSYS_DYNSYSLIB_H_ 

#include "capd/basicalg/factrial.h"
#include "capd/map/Map.h"
#include "capd/map/C2Map.h"
#include "capd/map/CnMap.h"
#include "capd/dynsys/SolverException.h"
#include "capd/dynsys/C2Solver.h"
#include "capd/dynsys/CnSolver.h"
#include "capd/dynsys/Linear2d.h"
#include "capd/dynsys/Linear3d.h"
#include "capd/dynsys/VLin3D.h"
#include "capd/dynsys/OdeNumTaylor.h"
#include "capd/vectalg/vectalgLib.h"
#include "capd/map/mapLib.h"

typedef capd::dynsys::DynSys<IMatrix> DynSys;
typedef capd::dynsys::Ode<IMatrix> Ode;
typedef capd::dynsys::OdeNum<IMatrix> OdeNum;
typedef capd::dynsys::OdeNumTaylor<IMatrix> OdeNumTaylor;
typedef capd::dynsys::VLin3D<IMatrix> VLin3D;
typedef capd::dynsys::Linear2d<IMatrix> linear2d;
typedef capd::dynsys::Linear3d<IMatrix> linear3d;
typedef capd::dynsys::Solver<IMap> Taylor;
typedef capd::dynsys::C2Solver<IC2Map> C2Taylor;
typedef capd::dynsys::CnSolver<ICnMap> CnTaylor;
typedef capd::dynsys::SolverException<IVector> TaylorException;


typedef capd::dynsys::DynSys<IMatrixMD> DynSysMD;
typedef capd::dynsys::Ode<IMatrixMD> OdeMD;
typedef capd::dynsys::OdeNum<IMatrixMD> OdeNumMD;
typedef capd::dynsys::OdeNumTaylor<IMatrixMD> OdeNumTaylorMD;
typedef capd::dynsys::VLin3D<IMatrixMD> VLin3DMD;
typedef capd::dynsys::Linear2d<IMatrixMD> linear2dMD;
typedef capd::dynsys::Linear3d<IMatrixMD> linear3dMD;
typedef capd::dynsys::Solver<IMapMD> TaylorMD;
typedef capd::dynsys::C2Solver<IC2MapMD> C2TaylorMD;
typedef capd::dynsys::CnSolver<ICnMapMD> CnTaylorMD;
typedef capd::dynsys::SolverException<IVectorMD> TaylorExceptionMD;

// classes for nonrigorous computations

typedef capd::dynsys::BasicSolver<DMap> BasicTaylor;
typedef capd::dynsys::BasicC2Solver<DC2Map> BasicC2Taylor;
typedef capd::dynsys::BasicCnSolver<DCnMap> BasicCnTaylor;

typedef capd::dynsys::BasicSolver<LDMap> LDTaylor;
typedef capd::dynsys::BasicC2Solver<LDC2Map> LDC2Taylor;
typedef capd::dynsys::BasicCnSolver<LDCnMap> LDCnTaylor;

typedef capd::dynsys::BasicSolver<DMapMD> BasicTaylorMD;
typedef capd::dynsys::BasicC2Solver<DC2MapMD> BasicC2TaylorMD;
typedef capd::dynsys::BasicCnSolver<DCnMapMD> BasicCnTaylorMD;

typedef capd::dynsys::BasicSolver<LDMapMD> LDTaylorMD;
typedef capd::dynsys::BasicC2Solver<LDC2MapMD> LDC2TaylorMD;
typedef capd::dynsys::BasicCnSolver<LDCnMapMD> LDCnTaylorMD;



#if (CAPD_CPU_ARCH==CAPD_CPU_ARCH_X86)

  typedef capd::dynsys::DynSys<LIMatrix> LIDynSys;
  typedef capd::dynsys::Ode<LIMatrix> LIOde;
  typedef capd::dynsys::OdeNum<LIMatrix> LIOdeNum;
  typedef capd::dynsys::OdeNumTaylor<LIMatrix> LIOdeNumTaylor;
  typedef capd::dynsys::VLin3D<LIMatrix> LIVLin3D;
  typedef capd::dynsys::Linear2d<LIMatrix> LIlinear2d;
  typedef capd::dynsys::Linear3d<LIMatrix> LIlinear3d;
  typedef capd::dynsys::Solver<LIMap> LITaylor;
  typedef capd::dynsys::C2Solver<LIC2Map> LIC2Taylor;
  typedef capd::dynsys::CnSolver<LICnMap> LICnTaylor;
  typedef capd::dynsys::SolverException<LIVector> LITaylorException;

  typedef capd::dynsys::DynSys<LIMatrixMD> LIDynSysMD;
  typedef capd::dynsys::Ode<LIMatrixMD> LIOdeMD;
  typedef capd::dynsys::OdeNum<LIMatrixMD> LIOdeNumMD;
  typedef capd::dynsys::OdeNumTaylor<LIMatrixMD> LIOdeNumTaylorMD;
  typedef capd::dynsys::VLin3D<LIMatrixMD> LIVLin3DMD;
  typedef capd::dynsys::Linear2d<LIMatrixMD> LIlinear2dMD;
  typedef capd::dynsys::Linear3d<LIMatrixMD> LIlinear3dMD;
  typedef capd::dynsys::Solver<LIMapMD> LITaylorMD;
  typedef capd::dynsys::C2Solver<LIC2MapMD> LIC2TaylorMD;
  typedef capd::dynsys::CnSolver<LICnMapMD> LICnTaylorMD;
  typedef capd::dynsys::SolverException<LIVectorMD> LITaylorExceptionMD;

#endif



#ifdef __HAVE_MPFR__

  typedef capd::dynsys::DynSys<MpIMatrix> MpDynSys;
  typedef capd::dynsys::Ode<MpIMatrix> MpOde;
  typedef capd::dynsys::OdeNum<MpIMatrix> MpOdeNum;
  typedef capd::dynsys::OdeNumTaylor<MpIMatrix> MpOdeNumTaylor;
  typedef capd::dynsys::VLin3D<MpIMatrix> MpVLin3D;
  typedef capd::dynsys::Linear2d<MpIMatrix> MpLinear2d;
  typedef capd::dynsys::Linear3d<MpIMatrix> MpLinear3d;
  typedef capd::dynsys::Solver<MpIMap> MpTaylor;
  typedef capd::dynsys::C2Solver<MpIC2Map> MpC2Taylor;
  typedef capd::dynsys::CnSolver<MpICnMap> MpCnTaylor;
  typedef capd::dynsys::SolverException<MpIVector> MpTaylorException;

  // classes for nonrigorous computations
  typedef capd::dynsys::BasicSolver<MpMap> MpBasicTaylor;
  typedef capd::dynsys::BasicC2Solver<MpC2Map> MpBasicC2Taylor;
  typedef capd::dynsys::BasicCnSolver<MpCnMap> MpBasicCnTaylor;



  typedef capd::dynsys::DynSys<MpIMatrixMD> MpDynSysMD;
  typedef capd::dynsys::Ode<MpIMatrixMD> MpOdeMD;
  typedef capd::dynsys::OdeNum<MpIMatrixMD> MpOdeNumMD;
  typedef capd::dynsys::OdeNumTaylor<MpIMatrixMD> MpOdeNumTaylorMD;
  typedef capd::dynsys::VLin3D<MpIMatrixMD> MpVLin3DMD;
  typedef capd::dynsys::Linear2d<MpIMatrixMD> MpLinear2dMD;
  typedef capd::dynsys::Linear3d<MpIMatrixMD> MpLinear3dMD;
  typedef capd::dynsys::Solver<MpIMapMD> MpTaylorMD;
  typedef capd::dynsys::C2Solver<MpIC2MapMD> MpC2TaylorMD;
  typedef capd::dynsys::CnSolver<MpICnMapMD> MpCnTaylorMD;
  typedef capd::dynsys::SolverException<MpIVectorMD> MpTaylorExceptionMD;

  // classes for nonrigorous computations
  typedef capd::dynsys::BasicSolver<MpMapMD> MpBasicTaylorMD;
  typedef capd::dynsys::BasicC2Solver<MpC2MapMD> MpBasicC2TaylorMD;
  typedef capd::dynsys::BasicCnSolver<MpCnMapMD> MpBasicCnTaylorMD;

#endif

#endif // _CAPD_DYNSYS_DYNSYSLIB_H_ 

/// @}
