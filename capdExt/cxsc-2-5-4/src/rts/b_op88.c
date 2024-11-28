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

/* CVS $Id: b_op88.c,v 1.21 2014/01/30 17:24:04 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_op88.c                              */
/*                                                              */
/*      Entry           : a_intg b_op88(desc,name,level)        */
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
/*      Description     : Open a system text file for reading.  */
/*                                                              */
/*      Note            : File must be a text file. (unchecked) */
/*                        Temporary string 'name' not freed.    */
/*                        No get on opened file performed.      */
/*                        Filename not longer than 'f_fnsz-1'.  */
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
extern char f_name[];
extern char o_sdir[];
extern char o_udir[];
extern char *o_text[];
#endif

#define cd_level                90
#define PXSC_USR_level          60
#define usrdir_level            50
#define PXSC_SYS_level          20
#define sysdir_level            10
#define not_found_level          0

#ifdef LINT_ARGS
local a_intg b_op88(f_text *desc,s_trng name,a_intg level)
#else
local a_intg b_op88(desc,name,level)

f_text *desc;
s_trng name;
a_intg level;
#endif

        {
        char *p,*filename;
        size_t len,len1,size;

        E_TPUSH("b_op88")

        /* copy filename to top of buffer and add '\0'  */
        size = (f_fnsz-1)-name.clen;
        filename = &f_name[size];
        filename[name.clen] = '\0';
        memcpy(filename,name.ptr,name.clen);

        desc->fp = NULL;

        /* search in current directory */
        if (level>=cd_level)
           {
           if ((desc->fp = fopen(filename,"r"))!=NULL)
              {
              len = 0;
              level = cd_level;
              goto set_flags;
              }
           }

        /* search in user directory */
        if (level>=PXSC_USR_level)
           {
           if ((p = getenv(o_text[4]))!=NULL)
              {
              if ((len = strlen(p))<=size)
                 {
                 memcpy(filename-len,p,len);
                 if ((desc->fp = fopen(filename-len,"r"))!=NULL)
                    {
                    level = PXSC_USR_level;
                    goto set_flags;
                    }
                 }
              }
           }

        /* search in user directory */
        if (level>=usrdir_level)
           {
           if ((len = strlen(o_udir))<=size)
              {
              memcpy(filename-len,o_udir,len);
              if ((p = getenv(o_text[3]))!=NULL)
                 len1 = strlen(p);
              else
                 len1 = 0;

              if ((len += len1)<=size)
                 {
                 if (len1) memcpy(filename-len,p,len1);
                 if ((desc->fp = fopen(filename-len,"r"))!=NULL)
                    {
                    level = usrdir_level;
                    goto set_flags;
                    }
                 }
              }
           }

        /* search in system directory */
        if (level>=PXSC_SYS_level)
           {
           if ((p = getenv(o_text[5]))!=NULL)
              {
              if ((len = strlen(p))<=size)
                 {
                 memcpy(filename-len,p,len);
                 if ((desc->fp = fopen(filename-len,"r"))!=NULL)
                    {
                    level = PXSC_SYS_level;
                    goto set_flags;
                    }
                 }
              }
           }

        /* search in system directory */
        if (level>=sysdir_level)
           {
           if ((len = strlen(o_sdir))<=size)
              {
              memcpy(filename-len,o_sdir,len);
              if ((desc->fp = fopen(filename-len,"r"))!=NULL)
                 {
                 level = sysdir_level;
                 goto set_flags;
                 }
              }
           }

        level = not_found_level;

set_flags:
        desc->outf = desc->stdi = desc->stdo = desc->temp = FALSE;
        if (desc->fp==NULL)
           {
           desc->eof = desc->eoln = desc->err = TRUE;
           desc->infl = desc->asgd = FALSE;
           level = not_found_level;
           }
        else
           {
           desc->eof = desc->eoln = desc->err = FALSE;
           desc->infl = desc->asgd = TRUE;
           memcpy(desc->name,filename,name.clen);
           desc->name[name.clen] = '\0';
           }

        E_TPOPP("b_op88")
        return(level);
        }





