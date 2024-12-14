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

/* CVS $Id: a_defs.h,v 1.21 2014/01/30 17:24:01 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : a_defs.h                              */
/*                                                              */
/****************************************************************/

#ifdef IEEE_HARDWARE
#if IBM_RT_C
#include <sys/FP.h>
#undef MININT
#define ROUND_DOWN              FP_DOWN
#define ROUND_NEAR              FP_NEAR
#define ROUND_UP                FP_UP
#else
#define ROUND_DOWN              (-1)
#define ROUND_NEAR              0
#define ROUND_UP                1
#endif
#endif

/* Rounding modes used by tenbyte arithemtic operations */
/*    DO NOT CHANGE THESE DEFINITIONS    */
#define DOWN            (-UP)
#define NEAREST         0
#define UP              1
#define CHOP            2





