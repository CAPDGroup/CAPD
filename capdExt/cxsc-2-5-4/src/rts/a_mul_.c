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

/* CVS $Id: a_mul_.c,v 1.21 2014/01/30 17:24:02 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : a_mul_.c                              */
/*                                                              */
/*      Entries         : a_intg a_mul_(a,b)                    */
/*                        a_intg a;                             */
/*                        a_intg b;                             */
/*                                                              */
/*      Arguments       : a = factor of multiplication          */
/*                        b = factor of multiplication          */
/*                                                              */
/*      Description     : Multiplication of integer values      */
/*                                                              */
/*      Note            : 0 returned in case of error           */
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
local a_intg a_mul_(a_intg a,a_intg b)
#else
local a_intg a_mul_(a,b)

a_intg a;
a_intg b;
#endif
        {
        a_intg res,g;
        a_btyp f;

        E_TPUSH("a_mul_")

        if (a==0 || b==0) res = 0;
        else if (a==1) res = b;
        else if (b==1) res = a;
        else if (a!=MININT && b!=MININT)
           {
           res = ((a^b)<0) ? -1 : 1;
           g = (a<0) ? -a : a;
           f = (a_btyp)((b<0) ? -b : b);

           if (res==1)
              {
              res = (f & 1) ? g : 0;
              while ((f >>= 1)!=0)
                 {
                 if (MAXINT-g<g) break;
                 if (MAXINT-(g <<= 1)<res) break;
                 if (f & 1) res += g;
                 }
              }
           else
              {
              res = (f & 1) ? -g : 0;
              while ((f >>= 1)!=0)
                 {
                 if (MAXINT-g<g) break;
                 if (MININT+(g <<= 1)>res) break;
                 if (f & 1) res -= g;
                 }
              }
           if (f!=0)
              {
              e_trap(OVERFLOW,6,E_TMSG,15,E_TINT+E_TEXT(1),&a,
                                          E_TINT+E_TEXT(2),&b);
              res = 0;
              }
           }
        else
           {
           e_trap(OVERFLOW,6,E_TMSG,15,E_TINT+E_TEXT(1),&a,
                                       E_TINT+E_TEXT(2),&b);
           res = 0;
           }

        E_TPOPP("a_mul_")
        return(res);
        }





