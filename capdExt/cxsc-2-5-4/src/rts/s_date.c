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

/* CVS $Id: s_date.c,v 1.22 2014/01/30 17:24:13 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : s_date.c                              */
/*                                                              */
/*      Entries         : s_trng s_date(fmt)                    */
/*                        s_trng fmt                            */
/*                                                              */
/*      Arguments       : fmt  - format string                  */
/*                                                              */
/*      Description     : Return date and time from system      */
/*                        fmt-'' ctime result                   */
/*                        fmt='...' format string for strftime  */
/*                        function                              */
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
#if GNU_C+IBM_EMX_C+HP_9000_C+SUN4_GNU_C+SUN4_OS5_GNU_C
#ifdef CLOCKS_PER_SEC
#undef CLOCKS_PER_SEC
#endif
#include <time.h>
#endif
#if SUN4_GNU_C
/* only defined when c++ flag defined in 2.5.8 */
/* extern long unsigned int strftime (char *, long unsigned int , const char *, const struct tm *); */
extern long time (long *);
#endif
#ifdef LINT_ARGS
local s_trng s_date(s_trng fmt)
#else
local s_trng s_date(fmt)
s_trng fmt;
#endif
{
        time_t tsys;     /* date/time stucture */
        s_trng datestr;  /* p string structure */
        char *p;         /* help str           */
        size_t len;

        E_TPUSH("s_date")
        /*
	    * get time/date from system
	    */
        time(&tsys);

        if ( fmt.clen == 0 )        
	   {
	   	   	/*
	    		 * empty format string
	    		 * only simple date output
	    		 */
        		p = (char *) ctime(&tsys);

			/* skip cr , len should be 25*/
         	 	len = (p!=NULL) ? (strlen(p)-1) : 0; 
        }
	   else
	   {
			/* 
			 * formatted date/time requested
			 */
			 len = fmt.clen *2 +50; /* estimation for strlen */
			 if ((p = (char*) malloc(len * sizeof(char)))==NULL) len = 0;
			 else
			 {
				struct tm *tmnow;
				char help;
				tmnow = (struct tm*) gmtime( &tsys);
#ifdef HEAP_CHECK
b_geth((a_char *)&p,(a_char *)p,(a_char *)"s_date");
#endif
				/* real length of formatted string */
				help=fmt.ptr[fmt.clen];
				fmt.ptr[fmt.clen] = '\0';
			 	len = (size_t) strftime(p,len,(char*)fmt.ptr,tmnow);
				fmt.ptr[fmt.clen] = help;
			 }
	   }
        s_init(&datestr, len);

	   if (datestr.ptr!=NULL)
        {    /* this implies len > 0 */
        		memcpy((char*)datestr.ptr,p,len);
              	datestr.ptr[len] = '\0';
              	datestr.clen = len;
        }
	   else datestr.clen = 0;

	   /*
	    * free allocated buffer
	    */
	   if (fmt.clen!=0 && p!= NULL) 
	   {
#ifdef HEAP_CHECK
b_freh((a_char *)&p,(a_char *)p,(a_char *)"s_date");
#endif
			free(p);
	   }
	   /*
	    * clear temporary object
	    */
	   if (fmt.tmp) s_free(&fmt);
        E_TPOPP("s_date")
        return datestr;
}





