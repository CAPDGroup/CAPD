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

/* CVS $Id: f_rdr1.c,v 1.21 2014/01/30 17:24:07 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : f_rdr1.c                              */
/*                                                              */
/*      Entry           : void f_rdr1(desc,r)                   */
/*                        f_text *desc;                         */
/*                        a_real *r;                            */
/*                                                              */
/*      Arguments       : desc   - device descriptor            */
/*                        r      - real variable                */
/*                                                              */
/*      Description     : perform PASCAL read(real).            */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern a_bool f_ppcc;
extern a_bool e_ofie;
#endif

#ifdef LINT_ARGS
local void f_rdr1(f_text *desc,a_real *r)
#else
local void f_rdr1(desc,r)

f_text *desc;
a_real *r;
#endif
        {
        a_bool ie_flag;

        E_TPUSH("f_rdr1")

        if (b_text(desc,TRUE))
           {
           ie_flag = e_ofie;
           e_ofie = FALSE;
           r_read(desc->fp,r,(a_intg)0,(int)desc->win.ch[0]);
           f_getc(desc);
           if (e_ofie)
              {
              if (f_ppcc)
                 e_trap(NO_ERROR+E_EMSG+E_ECNT,6,E_TMSG,68,E_TMSG,60,
                        E_TDBL,r);
              e_ofie = TRUE;
              }
           else if (ie_flag) e_ofie = TRUE;
           }

        E_TPOPP("f_rdr1")
        return;
        }





