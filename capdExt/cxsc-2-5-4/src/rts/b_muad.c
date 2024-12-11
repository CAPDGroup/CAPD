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

/* CVS $Id: b_muad.c,v 1.22 2014/01/30 17:24:04 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_muad.c                              */
/*                                                              */
/*      Entry           : void b_muad(a,b,lang)                 */
/*                        a_btyp a,b,*lang;                     */
/*                                                              */
/*      Arguments       : a = first operand for product         */
/*                        b = second operand for product        */
/*                        lang = a_btyp array (position for     */
/*                               least significant product part)*/
/*                                                              */
/*      Description     : Multiply a_btyp digits and add        */
/*                        product to a_btyp array               */
/*                                                              */
/*      Note            : Carry is propagated within array lang */
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
local void b_muad(a_btyp a,a_btyp b,a_btyp *lang)
#else
local void b_muad(a,b,lang)

a_btyp a;
a_btyp b;
a_btyp *lang;
#endif
        {
        a_btyp ahigh,alow,bhigh,blow,chigh,clow,d,dummy;
        a_intg cout;

        ahigh = GETHIGH(a);
        alow  = GETLOW(a);
        bhigh = GETHIGH(b);
        blow  = GETLOW(b);

        /* product of cross terms                               */
        dummy = ahigh*blow;
        d = MOVEHIGH(dummy);

#if C_P_1
        clow = alow*blow;
        cout = (~clow<d) ? 1 : 0;
        clow += d;
#else
        b_addu(alow*blow,d,(a_intg)0,&clow,&cout);
#endif
        d = GETHIGH(dummy);
        chigh = ahigh*bhigh+d+cout;
        dummy = alow*bhigh;
        d = MOVEHIGH(dummy);
#if C_P_1
        cout = (~clow<d) ? 1 : 0;
        clow += d;
#else
        b_addu(clow,d,(a_intg)0,&clow,&cout);
#endif
        d = GETHIGH(dummy);
        chigh += d+cout;

        /* add product                                          */
#if C_P_1
        cout = (~*lang<clow) ? 1 : 0;
        *lang += clow;
        clow = lang[-1]+chigh+cout;
        if (~lang[-1]<chigh || (clow==ZERO && cout))
                b_addc(lang-2);
        lang[-1] = clow;
#else
        b_addu(*lang,clow,(a_intg)0,lang,&cout);
        b_addu(lang[-1],chigh,cout,lang-1,&cout);
        if (cout) b_addc(lang-2);
#endif
        return;
        }





