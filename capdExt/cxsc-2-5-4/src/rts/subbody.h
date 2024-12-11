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

/* CVS $Id: subbody.h,v 1.22 2014/01/30 17:24:14 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : subbody.h                             */
/*                                                              */
/*      Description     : Subtraction of IEEE numbers.          */
/*                                                              */
/*      Notes           : included by r_subn.c,r_subd.c,r_subu.c*/
/*                                                              */
/*      Include         : body.h                                */
/*                                                              */
/*                   handle inf-inf correctly                   */
/*                   handle qNaN-real correctly                 */
/****************************************************************/

/* Note: at entry of this routine vzb holds sign of -b */

        if (expoa>EXPO_MAX || expob>EXPO_MAX)
                {

                /* a = infinity or NaN                          */
                if (expoa>EXPO_MAX) {

                   /* a = infinity                              */
                   if (MANT_INFINITY(manta)) {

                      /* b = infinity or NaN                    */
                      if (expob>EXPO_MAX) {

                         /* a = b = infinity                    */
                         if (MANT_INFINITY(mantb)) {
                            if (vza!=vzb)
                               e_trap(INV_OP+E_IEEE,8,E_TMSG,9,
                                      E_TDBL+E_TEXT(1),&a,
                                      E_TDBL+E_TEXT(2),&b,
                                      E_TDBL+E_TRES,&a);
                            E_TPOPP("subbody")
                            return(a);
                            }

                         /* a = infinity b = NaN                */
                         if (SIGNALING(mantb[0]))
                            e_trap(INV_OP+E_IEEE,8,E_TMSG,5,
                                   E_TDBL+E_TEXT(1),&a,
                                   E_TDBL+E_TEXT(2),&b,
                                   E_TDBL+E_TRES,&b);
                         E_TPOPP("subbody")
                         return(b);
                         }

                      /* subtract number from infinity => a     */
                      E_TPOPP("subbody")
                      return(a);
                      }

                   /* a = NaN                                   */
                   else {
                      if (SIGNALING(manta[0]))
                         e_trap(INV_OP+E_IEEE,8,E_TMSG,5,
                                E_TDBL+E_TEXT(1),&a,
                                E_TDBL+E_TEXT(2),&b,E_TDBL+E_TRES,&a);

                      /* a = qNaN     b = sNaN                  */
                      else if (expob>EXPO_MAX && !MANT_INFINITY(mantb) && 
                               SIGNALING(mantb[0]))
                         {
                         e_trap(INV_OP+E_IEEE,8,E_TMSG,5,
                                E_TDBL+E_TEXT(1),&a,
                                E_TDBL+E_TEXT(2),&b,E_TDBL+E_TRES,&b);
                         E_TPOPP("subbody")
                         return(b);
                         }

                      E_TPOPP("subbody")
                      return(a);
                      }
                   }

                /* b = infinity or NaN                          */
                else {

                   /* b = infinity                              */
                   if (MANT_INFINITY(mantb)) {
                      b_comp(&b,expob,mantb,vzb);
                      E_TPOPP("subbody")
                      return(b);
                      }

                   /* b = NaN                                   */
                   if (SIGNALING(mantb[0]))
                      e_trap(INV_OP+E_IEEE,8,E_TMSG,5,
                             E_TDBL+E_TEXT(1),&a,
                             E_TDBL+E_TEXT(2),&b,E_TDBL+E_TRES,&b);
                   E_TPOPP("subbody")
                   return(b);
                   }
                }

        /* check for one operand to be zero */
        if (zeroa)
           {

           /* special case (-0)-(-0),0-0 */
           if (zerob && vza!=vzb) vzb = vzz;
           
           b_comp(&b,expob,mantb,vzb);

           E_TPOPP("subbody")
           return(b);
           }
        else if (zerob)
           {
           E_TPOPP("subbody")
           return(a);
           }

        /* change sign of real number b (Now do addition !)     */
        R_ASSIGN(oldb,b);
        b_comp(&b,expob,mantb,vzb);

/* equivalent part for addition and subtraction                 */
#ifdef AIX
#include "/u/p88c/runtime/real/body.h"
#else
#include "body.h"
#endif





