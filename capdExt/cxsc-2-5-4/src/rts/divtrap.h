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

/* CVS $Id: divtrap.h,v 1.22 2014/01/30 17:24:06 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : divtrap.h                             */
/*                                                              */
/*      Description     : NaN or infinity operator in division  */
/*                                                              */
/*                   union q removed                            */
/****************************************************************/

        if ( expoa>EXPO_MAX || expob>EXPO_MAX )
                {

                /* a = infinity or NaN                          */
                if (expoa==EXPO_MAX+1) {

                   /* a = infinity                              */
                   if (MANT_INFINITY(manta)) {

                      /* b = infinity or NaN                    */
                      if (expob==EXPO_MAX+1) {

                         /* a = infinity b = NaN                */
                         if (!MANT_INFINITY(mantb)) {
                            if (SIGNALING(mantb[0]))
                               {
                               e_trap(INV_OP+E_IEEE,8,E_TMSG,5,
                                         E_TDBL+E_TEXT(1),&a,
                                         E_TDBL+E_TEXT(2),&b,
                                         E_TDBL+E_TRES,&b);
                               }
                            E_TPOPP("divtrap")
                            return(b);
                            }

                         /* a = b = infinity                    */
                         e_trap(INV_OP+E_IEEE,8,E_TMSG,4,
                                   E_TDBL+E_TEXT(1),&a,
                                   E_TDBL+E_TEXT(2),&b,
                                   E_TDBL+E_TRES,&a);
                         E_TPOPP("divtrap")
                         return(a);
                         }

                      /* infinity divided by number             */
                      if (vza!=vzc) b_comp(&a,expoa,manta,vzc);
                      E_TPOPP("divtrap")
                      return(a);
                      }

                   /* a = NaN                                   */
                   else {
                      if (SIGNALING(manta[0]))
                         {
                         e_trap(INV_OP+E_IEEE,8,E_TMSG,5,
                                   E_TDBL+E_TEXT(1),&a,
                                   E_TDBL+E_TEXT(2),&b,E_TDBL+E_TRES,&a);
                         }
                      else if (expob==EXPO_MAX+1 && SIGNALING(mantb[0]))
                         {
                         e_trap(INV_OP+E_IEEE,8,E_TMSG,5,
                                   E_TDBL+E_TEXT(1),&a,
                                   E_TDBL+E_TEXT(2),&b,E_TDBL+E_TRES,&b);
                         }
                      E_TPOPP("divtrap")
                      return(a);
                      }
                   }

                /* b = infinity or NaN                          */
                else {

                   /* b = infinity                              */
                   if (MANT_INFINITY(mantb))
                      {
                      /* Generate +/- zero                      */
                      R_ASSIGN(a,((vzc) ? *r_sero : *r_zero));
                      E_TPOPP("divtrap")
                      return(a);
                      }

                   /* b = NaN                                   */
                   if (SIGNALING(mantb[0]))
                      {
                      e_trap(INV_OP+E_IEEE,8,E_TMSG,5,
                                E_TDBL+E_TEXT(1),&a,
                                E_TDBL+E_TEXT(2),&b,E_TDBL+E_TRES,&b);
                      }
                   E_TPOPP("divtrap")
                   return(b);
                   }
                }





