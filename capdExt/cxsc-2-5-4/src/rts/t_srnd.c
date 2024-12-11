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

/* CVS $Id: t_srnd.c,v 1.22 2014/01/30 17:24:17 cxsc Exp $ */

/*****************************************************************
**                                                              **
**  t_srnd.c   06.02.91 Baumhof                                 **
**                                                              **
**  Rundungsauswahl                                             **
**                                                              **
**  Format:                                                     **
**      int  t_srnd (int rnd);  alle viere                      **
**      int  t_grnd ;           aktuelle Rundung abfragen       **
**                                                              **
**      Date    : 1992-07-21                                    **
*****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
extern int b_rflg;

#ifdef IEEE_HARDWARE
#if SUN4_OS4_C+SUN4_GNU_C+SUN4_CPP_C
/* #define __CB_SUN4_OS4_HARDWARE */
#endif
#endif


#ifdef IEEE_HARDWARE

#ifdef __CB_SUN4_OS4_HARDWARE

#define __T_SRND_DEFINED

extern unsigned int e_srq_, e_srqn, e_srqd, e_srqu, e_srqc;

#ifdef LINT_ARGS
int t_srnd(int rnd)
#else
int t_srnd(rnd)
int rnd;
#endif
{
    if (rnd==DOWN || rnd==UP || rnd==NEAREST || rnd==CHOP)
    {
        b_rflg = rnd;
	switch (rnd)
	{
	    case NEAREST:	e_srq_ = e_srqn; break;
	    case UP:		e_srq_ = e_srqu; break;
	    case DOWN:		e_srq_ = e_srqd; break;
	    case CHOP:		e_srq_ = e_srqc; break;
	}
        return(0);
    }
    else
    {
        b_rflg = NEAREST;
	e_srq_ = e_srqn;
        return(1);
    }
}

#endif	/* __CB_SUN4_OS4_HARDWARE */

#if IBM_LINUX_C+IBM_EMX_C

#define __T_SRND_DEFINED

extern unsigned int e_cwe_, e_cwen, e_cwed, e_cweu, e_cwec;

#ifdef LINT_ARGS
int t_srnd(int rnd)
#else
int t_srnd(rnd)
int rnd;
#endif
{
    if (rnd==DOWN || rnd==UP || rnd==NEAREST || rnd==CHOP)
    {
        b_rflg = rnd;
	switch (rnd)
	{
	    case NEAREST:	e_cwe_ = e_cwen; break;
	    case UP:		e_cwe_ = e_cweu; break;
	    case DOWN:		e_cwe_ = e_cwed; break;
	    case CHOP:		e_cwe_ = e_cwec; break;
	}
        return(0);
    }
    else
    {
        b_rflg = NEAREST;
	e_cwe_ = e_cwen;
        return(1);
    }
}

#endif	/* IBM_LINUX_C */

#endif	/* IEEE_HARDWARE */

#ifndef __T_SRND_DEFINED

#ifdef LINT_ARGS
int t_srnd(int rnd)
#else
int t_srnd(rnd)
int rnd;
#endif
{
    if (rnd==DOWN || rnd==UP || rnd==NEAREST || rnd==CHOP)
    {
        b_rflg = rnd;
        return(0);
    }
    else
    {
        b_rflg = NEAREST;
        return(1);
    }
}

#endif	/* __T_SRND_DEFINED */


#ifdef LINT_ARGS
int t_grnd(void)
#else
int t_grnd()
#endif
{
    return(b_rflg);
}





