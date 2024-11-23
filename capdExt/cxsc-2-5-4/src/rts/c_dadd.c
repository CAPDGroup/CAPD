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

/* CVS $Id: c_dadd.c,v 1.21 2014/01/30 17:24:05 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : c_dadd.c                              */
/*                                                              */
/*      Entries         : void c_dadd(cr,ci,a)                  */
/*                        dotprecision *cr,*ci;                 */
/*                        d_otpc a;                             */
/*                                                              */
/*      Arguments       : cr= dotprecision variable(real)       */
/*                        ci= dotprecision variable(imaginary)  */
/*                        a = complex dotprecision value        */
/*                                                              */
/*      Description     : Add complex to complex dotprecision   */
/*                        variable c = c+a                      */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
#endif

#ifdef LINT_ARGS
local void c_dadd(dotprecision *cr,dotprecision *ci,d_otpc a)
#else
local void c_dadd(cr,ci,a)

dotprecision *cr;
dotprecision *ci;
d_otpc a;
#endif
{
  E_TPUSH("c_dadd")

  d_otpc *ptemp= &a;

  d_vlcp((dotprecision*) &(ptemp)->RE);
  d_vlcp((dotprecision*) &(ptemp)->IM);

  d_dadd(cr,a.RE);
  d_dadd(ci,a.IM);

  ptemp= &a;
  d_free((dotprecision*) &(ptemp)->RE);
  d_free((dotprecision*) &(ptemp)->IM);

  E_TPOPP("c_dadd")
  return;
}





