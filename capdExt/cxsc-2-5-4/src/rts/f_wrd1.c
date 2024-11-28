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

/* CVS $Id: f_wrd1.c,v 1.21 2014/01/30 17:24:08 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : f_wrd1.c                              */
/*                                                              */
/*      Entry           : void f_wrd1(desc,r,mode)              */
/*                        f_text *desc;                         */
/*                        d_otpr r;                             */
/*                        a_char mode ;                         */
/*                                                              */
/*      Arguments       : desc   - device descriptor            */
/*                        r      - real value                   */
/*                                                              */
/*      Description     : performs accu hex dump                */
/*                        The first 6 words contain status      */
/*                        information. 64 Bit of guard digits   */
/*                        are build by A_LOWnan and a_Start.    */
/*                        The accu is left aligned. On the      */
/*                        right side there are 28 bits unused   */
/*                        f_wrd1 writes status information      */
/*                        In respect of mode the output is:     */
/*                        'x' : prints all words in the accu in */
/*                              lower case or upper case if 'X' */
/*                        'z' : prints only nonzero words in the*/ 
/*                              accu in lower case or upper if  */
/*                              'Z' is specified.               */
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
local void f_wrd1(f_text *desc,d_otpr r, a_char mode )
#else
local void f_wrd1(desc,r,mode )

f_text *desc;
d_otpr r;
a_char mode ;
#endif
{
        E_TPUSH("f_wrd1")


        if (b_text(desc,FALSE))
	   {
		 f_wrid(desc->fp,r,(char)mode);
	   }
        E_TPOPP("f_wrd1");
	   return ;
}





