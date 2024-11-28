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

/* CVS $Id: s_whex.c,v 1.21 2014/01/30 17:24:14 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : s_whex.c                              */
/*                                                              */
/*      Entry           : s_trng s_whex( r, mode)               */
/*                        a_real r;                             */
/*                        a_char mode;                          */
/*                                                              */
/*      Arguments       : r      - real value                   */
/*                                                              */
/*      Description     : write   real value to hex string      */
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

static char lc[] = "0123456789abcdef";
static char uc[] = "0123456789ABCDEF";

#ifdef LINT_ARGS
local s_trng s_whex(a_real r, a_char mode)
#else
local s_trng s_whex(r, mode)
a_real r;
a_char mode;
#endif
        {
        size_t i;
	   s_trng s;
	   unsigned int c;
        a_char *p;

        E_TPUSH("s_whex")

        if (mode!='x' && mode!='X')
              e_trap(INV_ARG  ,4,E_TMSG,57,E_TCHR,&mode);
	   else
        {
              s_init(&s,(size_t)2*sizeof(r));
		    if (s.ptr != NULL)
		    {
		    	s.clen=s.alen;
		    	s.tmp= TRUE;
#if INTEL
              	p = ((a_char *)&r)+(sizeof(r)-1);
              	for (i=0;i<sizeof(r);i++)
		    	{
		    	  	 c=(unsigned int) ((*p & 0xf0)>>4);
        			 s.ptr[2*i] = (mode=='x') ? lc[c] : uc[c];
        			 c=(unsigned int) (*p-- & 0x0f);
        			 s.ptr[2*i+1] = (mode=='x') ? lc[c] : uc[c];
#else
              	p = (a_char *)&r;
              	for (i=0;i<sizeof(r);i++)
		    	{
		    	  	 c=(unsigned int) ((*p & 0xf0)>>4);
        			 s.ptr[2*i] = (mode=='x') ? lc[c] : uc[c];
        			 c= (unsigned int)(*p++ & 0x0f);
        			 s.ptr[2*i+1] = (mode=='x') ? lc[c] : uc[c];
#endif
              	}
              }
        }
        E_TPOPP("s_whex")
        return s;
}





