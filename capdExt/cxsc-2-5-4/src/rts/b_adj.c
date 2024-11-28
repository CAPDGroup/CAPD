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

/* CVS $Id: b_adj.c,v 1.21 2014/01/30 17:24:02 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_adj.c                               */
/*                                                              */
/*      Entries         : a_btyp b_adj(lang,expo)               */
/*                        a_btyp *lang;                         */
/*                        a_intg *expo;                         */
/*                                                              */
/*      Arguments       : lang = result mantissa of length BSIZE*/
/*                        expo = exponent of number             */
/*                                                              */
/*      Function value  : exception code                        */
/*                        OVERFLOW = overflow detected          */
/*                        UNDERFLOW = underflow detected        */
/*                        INEXACT = shifting of last D_U_RATIO  */
/*                                  elements of lang causes     */
/*                                  loss of accuracy            */
/*                                                              */
/*      Description     : Adjust denormalized number and        */
/*                        test for inexact data.                */
/*                                                              */
/*      Notes           : no significant bits for IEEE double   */
/*                        format are lost.                      */
/*                                                              */
/*                   exponent set to EXPO_MIN if denormalized   */
/*                   changed test for underflow exception       */
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
extern a_bool e_efuf;
extern a_bool e_ofie;
extern a_bool e_ofuf;
extern a_bool e_ofof;
#endif

#ifdef LINT_ARGS
local a_btyp b_adj(a_btyp *lang,a_intg *expo)
#else
local a_btyp b_adj(lang,expo)

a_btyp *lang;
a_intg *expo;
#endif
        {
        a_intg i,k;
        a_btyp rc;

        E_TPUSH("b_adj")

        if (*expo>EXPO_MAX)
           {

           /* if overflow trap enabled, adjust exponent         */
           if (e_efof)
              {
              *expo -= EXPO_ADJUST;

              E_TPOPP("b_adj")
              return(OVERFLOW);
              }
           else
              {

              /* set all significant mantissa bits              */
              *lang = HIDDEN_BIT | MANT_HIGH;
              lang[1] = MAX_BASETYPE;
              lang[D_U_RATIO] = MSB;

              *expo = EXPO_MAX;
              e_ofof = TRUE;

              if (e_efie) return(INEXACT);
              e_ofie = TRUE;

              E_TPUSH("b_adj")
              return(NO_ERROR);
              }
           }
        else if (*expo<EXPO_MIN)
           {

           /* bits are lost due to denormalization */
           rc = NO_ERROR;
           for (i=D_U_RATIO;i<BSIZE;i++)
              {
              if (lang[i])
                 {
                 rc = UNDERFLOW;

                 /* if underflow trap enabled, adjust exponent */
                 if (e_efuf)
                    {
                    *expo += EXPO_ADJUST;

                    E_TPOPP("b_adj")
                    return(rc);
                    }

                 e_ofuf = TRUE;
                 break;
                 }
              }

           /* denormalized number will result or zero */
           k = (*expo<EXPO_MIN-MANTL) ? MANTL+1 : EXPO_MIN-*expo;
           b_shru(lang,BSIZE,k);

           if (rc==NO_ERROR)
              {
              for (i=D_U_RATIO;i<BSIZE;i++)
                 {
                 if (lang[i])
                    {
                    rc = UNDERFLOW;

                    /* if underflow trap enabled, adjust exponent */
                    if (e_efuf)
                       {
                       b_shlu(lang,BSIZE,k);
                       *expo += EXPO_ADJUST;

                       E_TPOPP("b_adj")
                       return(rc);
                       }

                    e_ofuf = TRUE;
                    break;
                    }
                 }
              }

           /* set least significant bit for the case that all bits */
           /* have been shifted out in case of detected UNDERFLOW  */
           else
              {
              lang[BSIZE-1] = LSB;
              }

           /* exponent set to minimum exponent */
           *expo = EXPO_MIN;
           }

        E_TPOPP("b_adj")
        return(NO_ERROR);
        }





