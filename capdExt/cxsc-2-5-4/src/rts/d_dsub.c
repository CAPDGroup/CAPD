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

/* CVS $Id: d_dsub.c,v 1.21 2014/01/30 17:24:06 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : d_dsub.c                              */
/*                                                              */
/*      Entries         : void d_dsub(a,b)                      */
/*                        d_otpr *a,b;                          */
/*                                                              */
/*      Arguments       : a = dotprecision variable             */
/*                        b = dotprecision value                */
/*                                                              */
/*      Description     : Subtraction of dotprecision values    */
/*                        a = a-b                               */
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
local void d_dsub(d_otpr *a,d_otpr b)
#else
local void d_dsub(a,b)

d_otpr *a;
d_otpr b;
#endif
        {
        a_intg aanf,aend,banf,bend,cout,i,j,a_ge_b;
        a_btyp *d,status;

        E_TPUSH("d_dsub")

        /* accu a is empty                                      */
        if ( (aanf = (a_intg)((*a)[A_BEGIN]))==0) {
           if (b[A_STATUS] & (A_PINFINITY+A_MINFINITY))
              {
              status = A_MZERO+A_PZERO;
              if (b[A_STATUS] & A_PINFINITY) status |= A_MINFINITY;
              if (b[A_STATUS] & A_MINFINITY) status |= A_PINFINITY;
              }
           else if (b[A_BEGIN]==ZERO)
              {
              status = (*a)[A_STATUS] & (A_MZERO+A_PZERO);
              if (b[A_STATUS] & A_MZERO) status |= A_PZERO;
              if (b[A_STATUS] & A_PZERO) status |= A_MZERO;
              }
           else
              status = A_MZERO+A_PZERO;
           d_ass(a,b);
           (*a)[A_SIGN] = 1-(*a)[A_SIGN];
           (*a)[A_STATUS] &= A_TEMPORARY+A_OWNTEMPORARY;
           (*a)[A_STATUS] |= status;
           }

        /* accu b is empty                                      */
        else if ( (banf = (a_intg)b[A_BEGIN])==0) {
           if (b[A_STATUS] & A_TEMPORARY) d_free(&b);
           }

        /* both accus are non-empty                             */
        else    {

           /* accus hold NaN or infinity                   */
           if (((*a)[A_STATUS] | b[A_STATUS]) &
               (A_QUIETNAN | A_MINFINITY | A_PINFINITY))
              {
              if ((*a)[A_STATUS] & A_QUIETNAN)
                 {
                 }
              else if (b[A_STATUS] & A_QUIETNAN)
                 {
                 (*a)[A_STATUS] |= A_QUIETNAN;
                 (*a)[A_LOWNAN] = b[A_LOWNAN];
                 }
              else if (((*a)[A_STATUS] & A_MINFINITY) &&
                       (b[A_STATUS] & A_MINFINITY))
                 {
                 e_trap(INV_OP+E_IEEE,6,E_TMSG,9,E_TDTP,a,E_TDTP,&b);
                 (*a)[A_STATUS] |= A_QUIETNAN;
                 (*a)[A_LOWNAN] = INV_OP;
                 }
              else if (((*a)[A_STATUS] & A_PINFINITY) &&
                       (b[A_STATUS] & A_PINFINITY))
                 {
                 e_trap(INV_OP+E_IEEE,6,E_TMSG,9,E_TDTP,a,E_TDTP,&b);
                 (*a)[A_STATUS] |= A_QUIETNAN;
                 (*a)[A_LOWNAN] = INV_OP;
                 }
              else if (b[A_STATUS] & A_QUIETNAN)
                 {
                 (*a)[A_STATUS] |= A_QUIETNAN;
                 (*a)[A_LOWNAN] = b[A_LOWNAN];
                 }
              else if (((*a)[A_STATUS] &
                       (A_PINFINITY | A_MINFINITY))==ZERO)
                 {
                 if (b[A_STATUS] & A_MINFINITY)
                    (*a)[A_STATUS] |= A_PINFINITY;
                 else
                    (*a)[A_STATUS] |= A_MINFINITY;
                 }
              if (b[A_STATUS] & A_TEMPORARY) d_free(&b);
              E_TPOPP("d_dsub")
              return;
              }

           /* determine A_BEGIN and A_END of result        */
           if (aanf>banf) (*a)[A_BEGIN] = (a_btyp)banf;
           if ( (aend = (*a)[A_END]) < (bend = b[A_END]) )
               (*a)[A_END] = bend;

           /* subtract dotprecision variables of different sign    */
           if ((*a)[A_SIGN]!=b[A_SIGN]) {
               if (b_addm(bend-banf+1,&((*a)[banf]),&b[banf])) {
                   b_addc(&((*a)[banf-1]));

                   /* carry extends used accu area         */
                   if ((*a)[(*a)[A_BEGIN]-1]) (*a)[A_BEGIN]--;
                   }
               }

           /* subtract dotprecision variables of same sign */
           else {

               /* compare dotprecision values              */
               if ( (a_ge_b = (aanf>banf) ? FALSE : TRUE)==TRUE) {
                   if (aanf==banf) {
                       for (i=aanf;i<=aend && i<=bend;i++)
                           if ((*a)[i]>b[i]) break;
                           else if ((*a)[i]<b[i]) {
                               a_ge_b = FALSE;
                               break;
                               }
                       if (i>aend && i<=bend) a_ge_b = FALSE;
                       }
                   }

               if (a_ge_b) {
                   if (b_subm(bend-banf+1,&((*a)[banf]),&b[banf]))
                       b_subc(&((*a)[banf-1]));
                   }
               else {
                   cout = 0;
                   for (i=bend;i>aend;i--) (*a)[i] = b[i];
                   for (i=aend;i>=aanf;i--)
                       b_subu(b[i],(*a)[i],cout,
                              &((*a)[i]),&cout);
                   for (j=aanf-1;j>=banf;j--) (*a)[j] = b[j];
                   (*a)[A_SIGN] = 1-(*a)[A_SIGN];
                   if (cout) b_subc(&((*a)[aanf-1]));
                   }

               /* eliminate leading zeros                  */
               d = &((*a)[(*a)[A_BEGIN]]);
               while (*d++==ZERO)
                   if (++((*a)[A_BEGIN])>(*a)[A_END]) {
                       (*a)[A_BEGIN] = (*a)[A_END] =
                                       (*a)[A_SIGN] = ZERO;
                       break;
                       }
               }

           /* eliminate trailing zeros                     */
           if ((*a)[A_BEGIN]) {
               d = &((*a)[(*a)[A_END]]);
               while (*d--==ZERO) (*a)[A_END]--;
               }

           /* free temporary dotprecison variables         */
           if (b[A_STATUS] & A_TEMPORARY) d_free(&b);
           }

        E_TPOPP("d_dsub")
        return;
        }





