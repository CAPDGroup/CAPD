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

/* CVS $Id: e_tprt.c,v 1.22 2014/01/30 17:24:06 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : e_tprt.c                              */
/*                                                              */
/*      Entry           : void e_tprt(e_argc,e_argv)            */
/*                        int e_argc;                           */
/*                        va_list e_argv;                       */
/*                                                              */
/*      Arguments       : e_argc - number of arguments          */
/*                        e_argv - reference to arguments       */
/*                                                              */
/*      Description     : Print variable argument list.         */
/*                                                              */
/*      Note            : e_argc must be even.                  */
/*                        e_argv is a list of pairs representing*/
/*                        integer values for the data type      */
/*                        specification and pointers to the     */
/*                        objects.                              */
/*                        Displaying stops if one result        */
/*                        argument is detected or the list ends.*/
/*                                                              */
/*                   do not display number if E_TRES is set.    */
/*                   TRAP_REAL uses r_writ() to display real.   */
/*                   variable p must be of type a_intg          */
/*                   f_errr replaces stderr                     */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern char *e_head;
extern f_text f_errr;
extern int e_rtyp;
extern a_VOID e_rptr;
#endif

#ifdef LINT_ARGS
local void e_tprt(int e_argc,va_list e_argv)
#else
local void e_tprt(e_argc,e_argv)
int e_argc;
va_list e_argv;
#endif
        {
        int i,k,res,type,msg;    /* !!! must be int !!! */

        e_rtyp = 0;
        e_rptr = NULL;

        for (i=0;i<e_argc;)
           {
           /* possible pointer allignment problems              */
           type = e_ref(int);

           /* Process message */
           if (type==E_TMSG)
              {
              /* possible pointer allignment problems           */
              type = e_ref(int);
              e_tmsg(type);
              e_argc -= 2;
              continue;
              }

           /* display header if not result argument */
           if ((res = type & E_TRES)==0)
              fprintf(f_errr.fp,"%s",e_head);

           /* display argument number */
           if ((msg=(type & E_TMSG))+res==0)
              fprintf(f_errr.fp,"%1d. ",i/2+1);

           /* display object according to specified type */
           switch (type & ~(E_TMSG+E_TRES))
              {

              case E_TCHR:
                 /* possible pointer allignment problems        */
                 e_rptr = e_ref(char *);
                 if (res) break;
                 if (msg) e_tmsg(msg);
                 else fprintf(f_errr.fp,"       char   : ");
                 fprintf(f_errr.fp,"0x%2.2x   '%c'\n",
                         *(char *)e_rptr,*(char *)e_rptr);
                 break;

              case E_TDBL:
                 /* possible pointer allignment problems        */
#if DEC_ALPHA_C
			  /* double has special alignment other than a_real */
                 e_rptr = (char *)e_ref(double *);
#else
                 e_rptr = (char *)e_ref(a_real *);
#endif
                 if (res) break;
                 if (msg) e_tmsg(msg);
                 else fprintf(f_errr.fp,"       real   : 0x");
                 for (k=0;k<D_U_RATIO;k++)
                    fprintf(f_errr.fp,"%8.8lx", (unsigned long)
                  ((a_btyp *)e_rptr)[B_HPART+((B_LPART-B_HPART)/
                                     (D_U_RATIO-1))*k]);
#ifdef TRAP_REAL
                 fprintf(f_errr.fp,"   ");
                 r_writ(f_errr.fp,*(a_real *)e_rptr,16,0,0);
#endif
                 fprintf(f_errr.fp,"\n");
                 break;

              case E_TLNG:
                 /* possible pointer allignment problems        */
                 e_rptr = (char *)e_ref(a_long *);
                 if (res) break;
                 if (msg) e_tmsg(msg);
                 else fprintf(f_errr.fp,"       long   : 0x");
                 for (k=0;k<sizeof(a_long);k++)
fprintf(f_errr.fp,"%2.2x",((unsigned char *)e_rptr)[k]);
                 fprintf(f_errr.fp,"\n");
                 break;

              case E_TDTP:
                 /* possible pointer allignment problems        */
                 e_rptr = (char *)e_ref(dotprecision *);
                 if (res) break;
                 if (msg) e_tmsg(msg);
                 else fprintf(f_errr.fp,"        dot   : ");
                 if ((*(dotprecision *)e_rptr)[A_STATUS] &
                     A_QUIETNAN)
                    fprintf(f_errr.fp,"qNaN = %8.8lx\n",(unsigned long)
                            (*(dotprecision *)e_rptr)[A_LOWNAN]);
                 else if ((*(dotprecision *)e_rptr)[A_STATUS] &
                          A_PINFINITY)
                    fprintf(f_errr.fp,"+infinity\n");
                 else if ((*(dotprecision *)e_rptr)[A_STATUS] &
                          A_MINFINITY)
                    fprintf(f_errr.fp,"-infinity\n");
                 else if ((*(dotprecision *)e_rptr)[A_BEGIN]==ZERO)
                      fprintf(f_errr.fp,"0\n");
                 else
                    {
                    fprintf(f_errr.fp,"sign='%c' ",
                       ((*(dotprecision *)e_rptr)[A_SIGN])
                        ? '-' : '+');
                    fprintf(f_errr.fp,"begin=%-2ld ",(long)
                       (*(dotprecision *)e_rptr)[A_BEGIN]);
                    fprintf(f_errr.fp,"end=%-2ld ",(long)
                       (*(dotprecision *)e_rptr)[A_END]);
                    for (k=0;k<=(*(dotprecision *)e_rptr)[A_END]-
                         (*(dotprecision *)e_rptr)[A_BEGIN];k++)
                       {
                       if ((k & 7)==0) fprintf(f_errr.fp,"\n%s   ",e_head);
                       fprintf(f_errr.fp,"%8.8lx ",(unsigned long)
                          (*(dotprecision *)e_rptr)
                          [(*(dotprecision *)e_rptr)[A_BEGIN]+k]);
                       }
                    fprintf(f_errr.fp,"\n");
                    }
                 break;

              case E_TINT:
                 /* possible pointer allignment problems        */
                 e_rptr = (char *)e_ref(a_intg *);
                 if (res) break;
                 if (msg) e_tmsg(msg);
                 else fprintf(f_errr.fp,"   integer    : ");
                 fprintf(f_errr.fp,"0x%8.8lx   ",(long int) *(a_intg *)e_rptr);
                 fprintf(f_errr.fp,"%10ld\n",(long int) *(a_intg *)e_rptr);
                 break;

              case E_TMLT:
                 /* possible pointer allignment problems        */
                 e_rptr = (char *)e_ref(multiprecision *);
                 if (res) break;
                 if (msg) e_tmsg(msg);
                 else fprintf(f_errr.fp,"     dynamic  : ");
                 if ((*(multiprecision *)e_rptr)->z)
                    fprintf(f_errr.fp,"0\n");
                 else
                    {
                    fprintf(f_errr.fp,"sign='%c' ",
                            ((*(multiprecision *)e_rptr)->s) ? '-' : '+');
                    fprintf(f_errr.fp,"exp=%-10ld ",(long)
                            (*(multiprecision *)e_rptr)->e);
                    fprintf(f_errr.fp,"len=%-10ld ",(long)
                            (*(multiprecision *)e_rptr)->l);
                    for (k=0;k<(*(multiprecision *)e_rptr)->l;k++)
                       {
                       if ((k & 7)==0) fprintf(f_errr.fp,"\n%s   ",e_head);
                       fprintf(f_errr.fp,"%8.8lx ", (unsigned long)
                               (*(multiprecision *)e_rptr)->m[k]);
                       }
                    fprintf(f_errr.fp,"\n");
                    }
                 break;

              case E_TSTR:
                 /* possible pointer allignment problems        */
                 e_rptr = e_ref(char *);
                 if (res) break;
                 if (msg) e_tmsg(msg);
                 else fprintf(f_errr.fp,"       string : ");
                 fprintf(f_errr.fp,"'%s'\n",(char *)e_rptr);
                 break;
/*
              case E_TULT:
                 e_rptr = (char *)e_ref(ultraprecision *);
                 if (res) break;
                 if (msg) e_tmsg(msg);
                 else fprintf(f_errr.fp,"       ultra  : ");
                 if ((*(ultraprecision *)e_rptr)->z)
                    fprintf(f_errr.fp,"0\n");
                 else
                    {
                    fprintf(f_errr.fp,"sign='%c' ",
                    ((*(ultraprecision *)e_rptr)->s) ? '-' : '+');
                    fprintf(f_errr.fp,"exp=%-10d ",
                    (*(ultraprecision *)e_rptr)->e);
                    fprintf(f_errr.fp,"len=%-10d ",
                    (*(ultraprecision *)e_rptr)->l);
                    for (k=0;k<(*(ultraprecision *)e_rptr)->l;k++)
                       {
                       if ((k & 7)==0)
                          fprintf(f_errr.fp,"\n%s   ",e_head);
                       fprintf(f_errr.fp,"%08.8lx ", (unsigned long)
                       (*(ultraprecision *)e_rptr)->m[k]);
                       }
                    fprintf(f_errr.fp,"\n");
                    }
                 break;
*/
              case E_TSTG:
                 /* possible pointer allignment problems        */
                 e_rptr = (char *)e_ref(s_trng *);
                 if (res) break;
                 if (msg) e_tmsg(msg);
                 else fprintf(f_errr.fp,"       string : ");
                 fprintf(f_errr.fp,"'");
                 for (k=0;k<((s_trng *)e_rptr)->clen;k++)
                    fprintf(f_errr.fp,"'%c'\n",
                    ((s_trng *)e_rptr)->ptr[k]);
                 fprintf(f_errr.fp,"'\n");
                 break;
              }
           i += 2;

           /* Display stops if result */
           if (res)
              {
              e_rtyp = type ^ E_TRES;
              return;
              }
           }

        return;
        }





