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

/* CVS $Id: a_ival.c,v 1.22 2014/01/30 17:24:02 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : a_ival.c                              */
/*                                                              */
/*      Entries         : a_intg a_ival(s)                      */
/*                        s_trng s;                             */
/*                                                              */
/*      Arguments       : s = input string                      */
/*                                                              */
/*      Description     : Convert integer literal to integer    */
/*                                                              */
/*      Note            : 0 returned in case of error           */
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
local a_intg a_ival(s_trng s)
#else
local a_intg a_ival(s)

s_trng s;
#endif
        {
        a_intg dig,res,len;
        a_bool sign;
        a_char ch;

        E_TPUSH("a_ival")

        len = res = 0;
        sign = FALSE;

        /* ignore leading blanks                             */
        while (len<s.clen)
           if (s.ptr[len]==' ') len++;
           else break;

        if (len==s.clen)
           {
           e_trap(I_O_ERROR,2,E_TMSG,62);
           }
        else
           {

           /* check for sign                                 */
           if (s.ptr[len]=='+') len++;
           else if (s.ptr[len]=='-')
              {
              len++;
              sign = TRUE;
              }

           if (len==s.clen)
              {
              e_trap(I_O_ERROR,4,E_TMSG,63,E_TSTG+E_TEXT(10),&s);
              }

           /* handle first decimal digit                     */
           else if (!isdigit((int)s.ptr[len]))
              {
              ch = s.ptr[len];
              e_trap(I_O_ERROR,4,E_TMSG,21,E_TCHR,&ch);
              }

           /* convert decimal string to integer              */
           else
              {
              res = s.ptr[len]-'0';

              while (++len<s.clen)
              if (isdigit((int)s.ptr[len]))
                 {
                 dig = s.ptr[len]-'0';

                 /* check for integer overflow               */
                 if ((MAXINT-dig)/10<res)
                    {
                    if (sign && (MININT+dig)/10==-res)
                       res = MININT;
                    else
                       {
                       if (sign) res = MININT;
                       else res = MAXINT;
                       while (++len<s.clen)
                          if (!isdigit((int)s.ptr[len])) break;
                       e_trap(OVERFLOW,4,E_TMSG,15,E_TSTG+E_TEXT(10),&s);
                       break;
                       }
                    }

                 /* add digit                                */
                 else
                    res = res*10+dig;
                 }
              else
                 break;

              if (sign && res!=MININT) res = -res;
              }
           }

        if (len>BUFFERSIZE && isdigit((int)s.ptr[len]))
           e_trap(I_O_BUFFER,2,E_TMSG,56);

        if (s.tmp) s_free(&s);

        E_TPOPP("a_ival")
        return(res);
        }





