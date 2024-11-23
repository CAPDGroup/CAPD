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

/* CVS $Id: b_addu.c,v 1.21 2014/01/30 17:24:02 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_addu.c                              */
/*                                                              */
/*      Entry           : void b_addu(a,b,cin,c,cout)           */
/*                        a_btyp a,b,*c;                        */
/*                        a_intg cin,*cout;                     */
/*                                                              */
/*      Arguments       : a = first operand of addition         */
/*                        b = second operand of addition        */
/*                        cin = carry in                        */
/*                        c = a+b+cin                           */
/*                        cout = carry out                      */
/*                                                              */
/*      Description     : Addition of a_btyp digits             */
/*                        with carry                            */
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
local void b_addu(a_btyp a,a_btyp b,a_intg cin,a_btyp *c,a_intg *cout)
#else
local void b_addu(a,b,cin,c,cout)

a_btyp a;
a_btyp b;
a_intg cin;
a_btyp *c;
a_intg *cout;
#endif
        {
#if C_P_1
        *c = a+b+cin;
        *cout = (~a<b || (*c==ZERO && cin)) ? 1 : 0;
#else
        /* a+b produce no carry                                 */
        if ( ~a>=b )
                {
                *c = a+b;

                /* add cin to (a+b)                             */
                if (cin && *c==MAX_BASETYPE ) {
                        *c    = ZERO;
                        *cout = 1;
                        }
                else
                        {
                        *c += cin;
                        *cout = 0;
                        }
                return;
                }

        /* a+b produces carry                                   */
        if (a & MSB)
            if (b & MSB) *c = ((a^MSB)+(b^MSB))+cin;
            else *c = (((a^MSB)+b)^MSB)+cin;
        else *c = (((b^MSB)+a)^MSB)+cin;
        *cout = 1;
#endif

        return;
        }





