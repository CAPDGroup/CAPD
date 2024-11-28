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

/* CVS $Id: f_rhex.c,v 1.22 2014/01/30 17:24:07 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : f_rhex.c                              */
/*                                                              */
/*      Entry           : void f_rhex(desc,r,mode)              */
/*                        f_text *desc;                         */
/*                        a_real *r;                            */
/*                        a_char mode;                          */
/*                                                              */
/*      Arguments       : desc   - device descriptor            */
/*                        r      - real variable                */
/*                        mode   - read mode                    */
/*                                 'x' = hexadecimal            */
/*                                                              */
/*      Description     : read real value in special format     */
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
local void f_rhex(f_text *desc,a_real *r,a_char mode)
#else
local void f_rhex(desc,r,mode)

f_text *desc;
a_real *r;
a_char mode;
#endif
        {
        int i;
        a_btyp b = ZERO;

        E_TPUSH("f_rhex")

        if (b_text(desc,TRUE))
           {
           while (desc->eof==FALSE && desc->win.ch[0]==' ')
              f_getc(desc);

           switch(mode)
              {
              case 'X' :
              case 'x' :
                 for (i=0;i<16;i++)
                    {
                    b <<= 4;
                    if (desc->eof==TRUE)
                       {
                       e_trap(I_O_ERROR,4,E_TMSG,20,
                              E_TSTR+E_TEXT(8),desc->name);
                       break;
                       }
                    else if (desc->eoln==TRUE)
                       {
                       e_trap(I_O_ERROR,4,E_TMSG,53,
                              E_TSTR+E_TEXT(8),desc->name);
                       break;
                       }
                    else if (isdigit((int)desc->win.ch[0]))
                       b += desc->win.ch[0]-'0';
                    else if (isalpha((int)desc->win.ch[0]))
                       b += toupper(desc->win.ch[0])-('A'-10);
                    else
                       {
                       e_trap(I_O_ERROR,4,E_TMSG,52,
                              E_TCHR+E_TEXT(10),&desc->win.ch[0]);
                       break;
                       }

                    f_getc(desc);

                    if (i==7)
                       {
                       ((a_btyp *)r)[B_HPART] = b;
                       b = ZERO;
                       }
                    else if (i==15) ((a_btyp *)r)[B_LPART] = b;
                    }
                 break;

              default :
                 e_trap(I_O_ERROR,4,E_TMSG,51,E_TCHR,&mode);
              }
           }

        E_TPOPP("f_rhex")
        return;
        }





