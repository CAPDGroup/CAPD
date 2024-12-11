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

/* CVS $Id: s_ins2.c,v 1.21 2014/01/30 17:24:14 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : s_ins2.c                              */
/*                                                              */
/*      Entries         : a_VOID s_ins2(s,rs,re)                */
/*                        s_etof s;                             */
/*                        a_intg rs,re;                         */
/*                                                              */
/*      Arguments       : s = set                               */
/*                        rs = start of element range           */
/*                        re = end of element range             */
/*                                                              */
/*      Function value  : pointer to updated set                */
/*                                                              */
/*      Description     : Add element range rs-re to set s.     */
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
local a_VOID s_ins2(s_etof s,a_intg rs,a_intg re)
#else
local a_VOID s_ins2(s,rs,re)

s_etof s;
a_intg rs;
a_intg re;
#endif
        {
        register a_intg i,k;
        static unsigned char bits[] =
           { 0xFF,0x7F,0x3F,0x1F,0x0F,0x07,0x03,0x01 };

        E_TPUSH("s_ins2")

        if (rs<0 || re>s_SIZE*BITS_PER_CHAR-1)
           e_trap(INDEX_RANGE,4,E_TINT+E_TEXT(5),&rs,
                                E_TINT+E_TEXT(6),&re);
        else
           {
           if ((i = rs>>LOG_BITS_PER_CHAR)<(k = re>>LOG_BITS_PER_CHAR))
              {
              s[i] |= bits[rs&MOD_BITS_PER_CHAR];
              while (++i<k) s[i] = 0xFF;
              s[i] |= ~(bits[re&MOD_BITS_PER_CHAR]>>1);
              }
           else
              {
              s[i] |= bits[rs&MOD_BITS_PER_CHAR] & 
                      ~(bits[re&MOD_BITS_PER_CHAR]>>1);
              }
           }

        E_TPOPP("s_ins2")
        return((a_VOID)&s[0]);
        }





