/*
**  CXSC is a C++ library for eXtended Scientific Computing (V 2.5.4)
**
**  Copyright (C) 1990-2000 Institut fuer Angewandte Mathematik,
**                          Universitaet Karlsruhe, Germany
**            (C) 2000-2014 Wiss. Rechnen/Softwaretechnologie
**                          Universitaet Wuppertal, Germany   
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Library General Public
**  License as published by the Free Software Foundation; either
**  version 2 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Library General Public License for more details.
**
**  You should have received a copy of the GNU Library General Public
**  License along with this library; if not, write to the Free
**  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/* CVS $Id: compiler.h,v 1.23 2014/01/30 17:23:44 cxsc Exp $ */

#ifndef _CXSC_COMPILER_H_INCLUDED
#define _CXSC_COMPILER_H_INCLUDED
/* Dieses Header-File soll Compiler-Inkompatiblitaeten ueberwinden */

#include "o_syst.h"

// #ifndef __GNUG__
// #ifndef bool
//   typedef int bool;
//   #define false 0
//   #define true  1
// #endif
// #endif

//  (INTEL)
#if INTEL 
#define LOWREAL 0
#define HIGHREAL 1
#else
#define LOWREAL 1
#define HIGHREAL 0
#endif

#endif // _CXSC_COMPILER_H_INCLUDED

