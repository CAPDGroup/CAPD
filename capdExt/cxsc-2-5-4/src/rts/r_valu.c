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

/* CVS $Id: r_valu.c,v 1.21 2014/01/30 17:24:12 cxsc Exp $ */


/****************************************************************/
/*                                                              */
/*      Filename        : r_valu.c                              */
/*                                                              */
/*      Entries         : a_real r_valu(x)                      */
/*                                                              */
/*      Arguments       : a_intg x;                             */
/*                                                              */
/*                        E_CLS0 - signaling NaN                */
/*                        E_CLS1 - quiet NaN                    */
/*                        E_CLS2 - -infinity                    */
/*                        E_CLS3 - negative normalized nonzero  */
/*                        E_CLS4 - negative denormalized        */
/*                        E_CLS5 - -0                           */
/*                        E_CLS6 - +0                           */
/*                        E_CLS7 - positive denormalized        */
/*                        E_CLS8 - positive normalized nonzero  */
/*                        E_CLS9 - +infinity                    */
/*                                                              */
/*      Description     : Determine a real value belonging      */
/*                        to the specified class of values.     */
/*                                                              */
/*	1992-06-09 : use macro SET_SIGNAL			*/
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern a_real *r_minf;
extern a_real *r_mmax;
extern a_real *r_meps;
extern a_real *r_sero;
extern a_real *r_zero;
extern a_real *r_eps_;
extern a_real *r_max_;
extern a_real *r_pinf;
#endif

#ifdef LINT_ARGS
local a_real r_valu(a_intg x)
#else
local a_real r_valu(x)

a_intg x;
#endif
        {
        a_real res;
         
        E_TPUSH("r_valu")

        switch (x)               
           {
           case E_CLS0: /* signaling NaN */
                   ((a_btyp *)&res)[B_HPART] = SET_SIGNAL(EXPO_MASK);;
                   ((a_btyp *)&res)[B_LPART] = MAX_BASETYPE;
                   break;
           case E_CLS1: /* quiet NaN */
                   ((a_btyp *)&res)[B_HPART] = SET_QUIET(EXPO_MASK);
                   ((a_btyp *)&res)[B_LPART] = MAX_BASETYPE;
                   break;
           case E_CLS2: /* minus infinity */
                   R_ASSIGN(res,*r_minf);
                   break;
           case E_CLS3: /* minus normalized */
                   R_ASSIGN(res,*r_mmax);
                   break;
           case E_CLS4: /* minus denormalized */
                   R_ASSIGN(res,*r_meps);
                   break;
           case E_CLS5: /* minus zero */
                   R_ASSIGN(res,*r_sero);
                   break;
           case E_CLS6: /* plus zero */
                   R_ASSIGN(res,*r_zero);
                   break;
           case E_CLS7: /* plus denormalized */
                   R_ASSIGN(res,*r_eps_);
                   break;
           case E_CLS8: /* plus normalized */
                   R_ASSIGN(res,*r_max_);
                   break;
           case E_CLS9: /* plus infinity */
                   R_ASSIGN(res,*r_pinf);
                   break;
           }

        E_TPOPP("r_valu")
        return(res);
        }





