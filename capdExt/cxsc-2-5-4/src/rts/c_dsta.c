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

/* CVS $Id: c_dsta.c,v 1.21 2014/01/30 17:24:05 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : c_dsta.c                              */
/*                                                              */
/*      Entries         : d_otpc c_dsta(cr,ci)                  */
/*                        dotprecision cr,ci;                   */
/*                                                              */
/*      Arguments       : cr= dotprecision variable(real)       */
/*                        ci= dotprecision variable(imaginary)  */
/*                                                              */
/*      Description     : Convert dotprecision to complex       */
/*                        dotprecision                          */
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
local d_otpc c_dsta (dotprecision cr,dotprecision ci)
#else
local d_otpc c_dsta ( cr,ci)

dotprecision cr;
dotprecision ci;
#endif
{
       d_otpc  res;
       d_otpc  *ptemp = &res;

       E_TPUSH("c_dsta")

       d_vlcp((dotprecision*) &cr);
       d_vlcp((dotprecision*) &ci);

       d_init((dotprecision*) &(ptemp)->RE);
       d_init((dotprecision*) &(ptemp)->IM);
  

        d_ass(&res.RE,cr);
        d_ass(&res.IM,ci);

        d_free((dotprecision*) &cr);
        d_free((dotprecision*) &ci);

        ptemp= &res;
        d_temp((dotprecision*) &(ptemp)->RE);
        d_temp((dotprecision*) &(ptemp)->IM);
        E_TPOPP("c_dsta")
        return res;
}





