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

/* CVS $Id: f_op88.c,v 1.21 2014/01/30 17:24:07 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : f_op88.c                              */
/*                                                              */
/*      Entry           : a_intg f_op88(desc,name,level)        */
/*                        f_text *desc;                         */
/*                        s_trng name;                          */
/*                        a_intg level;                         */
/*                                                              */
/*      Arguments       : desc   - descriptor of opened file    */
/*                        name   - filename to be searched for  */
/*                        level  - search level                 */
/*                                                              */
/*      Return value    : level on which files is found         */
/*                           0 = not found                      */
/*                          10 = found in sysdir                */
/*                          20 = found in PXSC_SYS              */
/*                          50 = found in usrdir                */
/*                          60 = found in PXSC_USR              */
/*                          90 = found in current directory     */
/*                                                              */
/*      Description     : Open a system file for reading.       */
/*                        File must be a text file.             */
/*                        Filename must not exceed 'f_fnsz-1'.  */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#define local static
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#undef local
#define local
#endif

#define not_found_level          0

#ifdef LINT_ARGS
local a_intg f_op88(f_text *desc,s_trng name,a_intg level)
#else
local a_intg f_op88(desc,name,level)

f_text *desc;
s_trng name;
a_intg level;
#endif
 
        {
        E_TPUSH("f_op88")

        if (desc->text==FALSE)
           {
           e_trap(I_O_ERROR,2,E_TMSG,19);
           level = not_found_level;
           }
        else if (name.clen==0)
           {
           e_trap(I_O_ERROR,2,E_TMSG,62);
           level = not_found_level;
           }
        else if (name.clen>=f_fnsz)
           {
           e_trap(I_O_ERROR,4,E_TMSG,30,E_TSTG+E_TEXT(8),&name);
           }
        else if ((level = b_op88(desc,name,level))>0)
           f_getc(desc);

        if (name.tmp) s_free(&name);

        E_TPOPP("f_op88")
        return(level);
        }





