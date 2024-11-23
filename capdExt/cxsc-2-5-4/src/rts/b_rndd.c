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

/* CVS $Id: b_rndd.c,v 1.21 2014/01/30 17:24:04 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_rndd.c                              */
/*                                                              */
/*      Entries         : a_btyp b_rndd(lang,expo,vz)           */
/*                        a_btyp *lang;                         */
/*                        a_intg *expo;                         */
/*                        a_bool vz;                            */
/*                                                              */
/*      Arguments       : lang = result mantissa of length BSIZE*/
/*                        expo = exponent of result             */
/*                        vz = sign of result                   */
/*                                                              */
/*      Description     : Round downwards to IEEE number by     */
/*                        changing value of lang+D_U_RATIO-1    */
/*                                                              */
/*      Function value  : error code                            */
/*                        OVERFLOW - overflow occurred          */
/*                        INEXACT - inexact result              */
/*                                                              */
/*                   unnecessary e_sofo() removed               */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern a_bool e_efie;
extern a_bool e_efof;
extern a_bool e_ofie;
extern a_bool e_ofof;
#endif

#ifdef LINT_ARGS
local a_btyp b_rndd(a_btyp *lang,a_intg *expo,a_bool vz)
#else
local a_btyp b_rndd(lang,expo,vz)

a_btyp *lang;
a_intg *expo;
a_bool vz;
#endif
        {
        a_intg k;

        /* number is positive (no rounding required) */
        if (NOT(vz))
           {

           /* test for inexact data */
           for (k=D_U_RATIO;k<BSIZE;k++)
              {
              if (lang[k])
                 {
                 if (e_efie) return(INEXACT);
                 e_ofie = TRUE;
                 return(NO_ERROR);
                 }
              }
           return(NO_ERROR);
           }

        /* add correction if trailing mantissa is non-zero */
        for (k=D_U_RATIO;k<BSIZE;k++)
           {
           if (lang[k])
              {

              /* propagate carry */
              b_addc(lang+(D_U_RATIO-1));

              /* normalization (only msb is set) */
              if (SHFT_MASK & *lang)
                 {
                 *lang = HIDDEN_BIT;
                 (*expo)++;

                 /* overflow occurred */
                 if (*expo>EXPO_MAX)
                    {
                    if (e_efof)
                       {
                       *expo -= EXPO_ADJUST;
                       return(OVERFLOW);
                       }
                    e_ofof = TRUE;
                    }
                 }
              if (e_efie) return(INEXACT);
              e_ofie = TRUE;
              return(NO_ERROR);
              }
           }

        return(NO_ERROR);
        }





