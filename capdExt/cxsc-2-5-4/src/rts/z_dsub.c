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

/* CVS $Id: z_dsub.c,v 1.21 2014/01/30 17:24:18 cxsc Exp $ */


/****************************************************************/
/*                                                              */
/*      Filename        : z_dsub.c                              */
/*                                                              */
/*      Entries         : void z_dsub(ZRI,ZII,ZRS,ZIS,A)        */
/*                        dotprecision *ZRI,*ZII,*ZRS,*ZIS;     */
/*                        d_otpz A;                             */
/*                                                              */
/*      Arguments       : ZRI=dotprecision variable(real lower) */
/*                        ZII=dotprecision variable(imag lower) */
/*                        ZRS=dotprecision variable(real upper) */
/*                        ZIS=dotprecision variable(imag upper) */
/*                        A = complex interval dotprec. value   */
/*                                                              */
/*      Description     : SUB complex interval to dotprecision  */
/*                        variable                              */
/*                        Z = Z-A                               */
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
local void z_dsub(dotprecision *ZRI,dotprecision *ZII,dotprecision *ZRS,
                  dotprecision *ZIS,d_otpz A)
#else
local void z_dsub(ZRI,ZII,ZRS,ZIS,A)

dotprecision *ZRI;
dotprecision *ZII;
dotprecision *ZRS;
dotprecision *ZIS;
d_otpz A;
#endif

{
        d_otpz *ptemp= &A;
        d_otpi *pctemp= &(ptemp)->RE;

        E_TPUSH("z_dsub")
        d_vlcp((dotprecision*) &(pctemp)->INF);
        d_vlcp((dotprecision*) &(pctemp)->SUP);
        
        pctemp= &(ptemp)->IM;
        d_vlcp((dotprecision*) &(pctemp)->INF);
        d_vlcp((dotprecision*) &(pctemp)->SUP);
        

        d_dsub(ZRI,A.RE.SUP);
        d_dsub(ZII,A.IM.SUP);
        d_dsub(ZRS,A.RE.INF);
        d_dsub(ZIS,A.IM.INF);

        ptemp= &A;
        pctemp= &(ptemp)->RE;
        d_free((dotprecision*) &(pctemp)->INF);
        d_free((dotprecision*) &(pctemp)->SUP);
         
        pctemp= &(ptemp)->IM;
        d_free((dotprecision*) &(pctemp)->INF);
        d_free((dotprecision*) &(pctemp)->SUP);
         
          
        E_TPOPP("z_dsub")
        return;
}





