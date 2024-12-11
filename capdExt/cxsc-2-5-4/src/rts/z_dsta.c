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

/* CVS $Id: z_dsta.c,v 1.21 2014/01/30 17:24:18 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : z_dsta.c                              */
/*                                                              */
/*      Entries         : d_otpz z_dsta(ZRI,ZII,ZRS,ZIS)        */
/*                        dotprecision ZRI,ZII,ZRS,ZIS;         */
/*                                                              */
/*      Arguments       : ZRI=dotprecision value(real lower)    */
/*                        ZII=dotprecision value(imag lower)    */
/*                        ZRS=dotprecision value(real upper)    */
/*                        ZIS=dotprecision value(imag upper)    */
/*                                                              */
/*      Description     : Convert complex interval dotprec.     */
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
local d_otpz z_dsta (dotprecision ZRI,dotprecision ZII,dotprecision ZRS,
                     dotprecision ZIS)
#else
local d_otpz z_dsta ( ZRI, ZII, ZRS, ZIS)

dotprecision ZRI;
dotprecision ZII;
dotprecision ZRS;
dotprecision ZIS;
#endif
{
        d_otpz  res;
        d_otpz *ptemp= &res;
        d_otpi *pctemp;

        E_TPUSH("z_dsta")

        d_vlcp(&ZRI);
        d_vlcp(&ZII);
        d_vlcp(&ZRS);
        d_vlcp(&ZIS);

        pctemp= &(ptemp)->RE;
        d_init((dotprecision*) &(pctemp)->INF);
        d_init((dotprecision*) &(pctemp)->SUP);

        pctemp= &(ptemp)->IM;
        d_init((dotprecision*) &(pctemp)->INF);
        d_init((dotprecision*) &(pctemp)->SUP);
        
        d_ass(&res.RE.INF,ZRI);
        d_ass(&res.IM.INF,ZII);
        d_ass(&res.RE.SUP,ZRS);
        d_ass(&res.IM.SUP,ZIS);

        d_free( &ZRI);
        d_free( &ZII);
        d_free( &ZRS);
        d_free( &ZIS);

        ptemp= &res;
        pctemp= &(ptemp)->RE;
        d_temp((dotprecision*) &(pctemp)->INF);
        d_temp((dotprecision*) &(pctemp)->SUP);

        pctemp= &(ptemp)->IM;
        d_temp((dotprecision*) &(pctemp)->INF);
        d_temp((dotprecision*) &(pctemp)->SUP);
        
        E_TPOPP("z_dsta")
        return res;
}





