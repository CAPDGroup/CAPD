/*
**  CXSC is a C++ library for eXtended Scientific Computing (V 2.5.4)
**
**  Copyright (C) 1990-2000 Institut fuer Angewandte Mathematik,
**                          Universitaet Karlsruhe, Germany
**            (C) 2000-2012 Wiss. Rechnen/Softwaretechnologie
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

/* CVS $Id: a_time.c,v 1.24 2014/02/27 15:49:21 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : a_time.c                              */
/*                                                              */
/*      Entries         : a_intg a_gtim()                       */
/*                      : void   a_itim()                       */
/*                                                              */
/*      Description     : Gets System CPU Time from clock       */
/*                        in ms sec, may wrap around on some    */
/*                        systems. Init timer, get timer        */
/*                        recognizes SUN systems                */
/*                                                              */
/*      Function value  : time between a_itim and a_gtim call   */
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

a_intg a_tvar; /* this var holds initial timer value from a_itim */


#if SUN4_OS4_C+DEC_ALPHA_C+IBM_LINUX_C+IBM_EMX_C+WINDOWS_X86_32
/* timer for sun4 os4, OSF/1  systems */
/* holds rusage structure */
#if WINDOWS_X86_32
#include <time.h>
struct timeval{
	long tv_sec;
	long tv_usec;
};
#else
#include <sys/time.h>
#include <sys/resource.h>
#endif

#if IBM_EMX_C+WINDOWS_X86_32
/* emx 0.9a does not support getrusage and rusage   */
/* but struct timeval and timezone and gettimeofday */
struct rusage {
	struct timeval ru_utime;
};

int getrusage(who, ruse)
int who;
struct rusage * ruse;
{
#if WINDOWS_X86_32
	return 0;
#else
	return gettimeofday(&(ruse->ru_utime),NULL);
#endif
}
#endif

#ifdef LINT_ARGS
local void   a_itim(void)
#else
local void   a_itim()
#endif
        {

        struct rusage buf;
	   int ret;

        E_TPUSH("a_itim")

        ret=getrusage(0,&buf); 
	   a_tvar = (a_intg)(buf.ru_utime.tv_sec*1000 + 
					 buf.ru_utime.tv_usec / 1000 ); 

        E_TPOPP("a_itim")
        return;
        }

#ifdef LINT_ARGS
local a_intg a_gtim(void)
#else
local a_intg a_gtim()
#endif
        {
        register a_intg res;
        struct rusage buf;

        E_TPUSH("a_gtim")

        res=getrusage(0,&buf); 
	   res = (a_intg)(buf.ru_utime.tv_sec*1000 - a_tvar + 
	                  buf.ru_utime.tv_usec / 1000 ); 

        E_TPOPP("a_gtim")
        return(res);
        }
 
#else
/* timer for non sun systems */
#ifdef LINT_ARGS
local void   a_itim(void)
#else
local void   a_itim()
#endif
        {

        E_TPUSH("a_itim")

        a_tvar = a_clck();

        E_TPOPP("a_itim")
        return;
        }

#ifdef LINT_ARGS
local a_intg a_gtim(void)
#else
local a_intg a_gtim()
#endif
        {
        register a_intg res;

        E_TPUSH("a_gtim")

        res = (a_clck()-a_tvar) / 1000;

        E_TPOPP("a_gtim")
        return(res);
        }
#endif







