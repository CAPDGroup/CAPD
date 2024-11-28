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

/* CVS $Id: b_popt.c,v 1.22 2014/01/30 17:24:04 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_popt.c                              */
/*                                                              */
/*      Entry           : int b_popt(device,option)             */
/*                        FILE *device;                         */
/*                        char *option;                         */
/*                                                              */
/*      Function value  : 0 - option is unknown / nothing done  */
/*                        1 - option found and processed/cont.  */
/*                        2 - option found and processed/exit   */
/*                                                              */
/*      Arguments       : device - output device                */
/*                        option - name of option               */
/*                                                              */
/*      Description     : Check for runtime option and process  */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#include "/u/p88c/runtime/o_revs.h"
#else
#include "o_defs.h"
#include "o_revs.h"
#endif
#define local
#ifdef COPYIB_ENABLE
extern a_bool f_ppib;
#endif
#ifdef DIRECT_ENABLE
extern char o_sdir[];
extern char o_udir[];
#endif
#ifdef SIGNED_ENABLE
extern a_bool f_ppsz;
#endif
#ifdef NORMALIZE_ENABLE
extern a_bool f_ppdn;
#endif
extern a_bool f_ppmt;
#ifdef STATTB_ENABLE
extern a_bool f_pptb;
#endif
#ifdef TRACEC_ENABLE
extern a_bool f_pptr;
#endif
extern a_bool f_ppcc;
extern a_bool f_pppl;
extern char *f_pplt[];
extern a_bool f_pptf;
extern char **f_argv;
extern char *o_text[];
#endif

#ifdef LINT_ARGS
local int b_popt(FILE *device,char *option)
#else
local int b_popt(device,option)

FILE *device;
char *option;
#endif

        {
        int i;     /* !!! must be int !!! */
#ifdef OPTION_ENABLE
        int ch,k;     /* !!! must be int !!! */
#endif
        a_bool key_found;
        char *flag,*key;
        size_t size;

        /* display active options only */
        key = o_text[12];
        size = strlen(option);
        while ((flag=strchr(++key,(int)'|'))!=NULL)
              {
              if (flag-key==size)
                 {
                 if (memcmp(key,option,size)==0)
                    {
                    fprintf(device,"%s(%s)(%s)(%s)(%s)(%s) ",
                       o_text[26],o_text[27],
                       o_text[29],o_text[30],
                       o_text[28],o_text[31]);
                    fprintf(device,"%s ",o_text[14]);
#ifdef COPYIB_ENABLE
                    fprintf(device,"%s ",o_text[15]);
#endif
#ifdef OPTION_ENABLE
                    fprintf(device,"%s(:(key))(@file) ",o_text[16]);
#endif
#ifdef NORMALIZE_ENABLE
                    fprintf(device,"%s ",o_text[17]);
#endif
                    fprintf(device,"%s ",o_text[18]);
                    fprintf(device,"%s ",o_text[13]);
#ifdef DIRECT_ENABLE
                    fprintf(device,"%s ",o_text[19]);
                    fprintf(device,"%s ",o_text[20]);
#endif
#ifdef SIGNED_ENABLE
                    fprintf(device,"%s ",o_text[21]);
#endif
#ifdef TRACEC_ENABLE
                    fprintf(device,"%s ",o_text[22]);
#endif
#ifdef STATTB_ENABLE
                    fprintf(device,"%s ",o_text[23]);
#endif
                    fprintf(device,"%s ",o_text[24]);
                    fprintf(device,"%s\n",o_text[25]);
                    return(2);
                    }
                 }
              key = flag;
              }

        /* display program parameter list only */
        if (strcmp(option,o_text[18])==0)
              {
              f_pppl = FALSE;
              for(i=0;i<FPPLT && f_pplt[i]!=NULL;i++)
                 fprintf(device,"%s ",f_pplt[i]);
              if (i>0) fprintf(device,"\n");
              return(2);
              }

#ifdef SIGNED_ENABLE
        /* signed zero is output */
        else if (strcmp(option,o_text[21])==0)
              {
              f_ppsz = NOT(f_ppsz);
              }
#endif

        /* display version number */
        else if (strcmp(option,o_text[25])==0)
              {
              fprintf(device,"%s",PRODUCT_TEXT);
              fprintf(device,"%s %s %s %s %s\n",VERSION_TEXT, REVISION_TEXT,
		            VERSION_SYS_TEXT, VERSION_EXT_TEXT, VERSION_ARI_TEXT );
              fprintf(device,"%s",COPYRIGHT_TEXT);
              }

        /* prompt for missing command line arguments */
        else if (strcmp(option,o_text[13])==0)
              {
              f_ppmt = NOT(f_ppmt);
              }

        /* generate temporary files instead of prompting */
        else if (strcmp(option,o_text[24])==0)
              {
              f_pptf = NOT(f_pptf);
              }

#ifdef TRACEC_ENABLE
        /* trace function calls pushed to trace back stack */
        else if (strcmp(option,o_text[22])==0)
              {
              f_pptr = NOT(f_pptr);
              }
#endif

#ifdef STATTB_ENABLE
        /* generate "static" trace back and no           */
        /* runtime function trace                        */
        else if (strcmp(option,o_text[23])==0)
              {
              f_pptb = NOT(f_pptb);
              }
#endif

        /* display message if inexact real conversion    */
        else if (strcmp(option,o_text[14])==0)
              {
              f_ppcc = NOT(f_ppcc);
              }

#ifdef NORMALIZE_ENABLE
        /* generate normalized successor/predecessor of 0.0 */
        else if (strcmp(option,o_text[17])==0)
              {
              f_ppdn = NOT(f_ppdn);
              }
#endif

#ifdef COPYIB_ENABLE
        /* copy index bounds of dynamic arrays           */
        else if (strcmp(option,o_text[15])==0)
              {
              f_ppib = NOT(f_ppib);
              }
#endif

#ifdef DIRECT_ENABLE
        /* display or overwrite variable sysdir          */
        else if (strncmp(option,o_text[19],strlen(o_text[19]))==0)
              {
              flag = option+strlen(o_text[19]);
              if (*flag==':' || *flag=='\0')
                 {
                 if (*flag==':') flag++;
                 size = strlen(flag);
                 if (size>0)
                    {
                    if (size>=DIRNAME_LENGTH)
                       fprintf(device,"SYSTEM DIRECTORY name too long.\n");
                    else
                       {
                       strncpy(o_sdir,flag,DIRNAME_LENGTH-1);
                       o_sdir[size] = '\0';
                       }
                    }
                 else
                    {
                    fprintf(device,"SYSTEM DIRECTORY = %s\n",o_sdir);
                    return(2);
                    }
                 }
              }

        /* display or overwrite variable usrdir          */
        else if (strncmp(option,o_text[20],strlen(o_text[20]))==0)
              {
              flag = option+strlen(o_text[20]);
              if (*flag==':' || *flag=='\0')
                 {
                 if (*flag==':') flag++;
                 size = strlen(flag);
                 if (size>0)
                    {
                    if (size>=DIRNAME_LENGTH)
                       fprintf(device,"USER DIRECTORY name too long.\n");
                    else
                       {
                       strncpy(o_udir,flag,DIRNAME_LENGTH-1);
                       o_udir[size] = '\0';
                       }
                    }
                 else
                    {
                    fprintf(device,"USER DIRECTORY = %s\n",o_udir);
                    return(2);
                    }
                 }
              }
#endif

        /* activate IEEE trap handling                   */
        else if (strncmp(option,o_text[26],strlen(o_text[26]))==0)
              {
              flag = option+strlen(o_text[26]);
              while (*flag)
                 {
                 if (*flag==o_text[27][0])
                    {
		    if (e_dz_e()) e_rdze(); else e_sdze();
                    }
                 else if (*flag==o_text[28][0])
                    { 
		    if (e_io_e()) e_rioe(); else e_sioe();
                    }
                 else if (*flag==o_text[29][0])
                    { 
		    if (e_of_e()) e_rofe(); else e_sofe();
                    }
                 else if (*flag==o_text[30][0])
                    { 
		    if (e_uf_e()) e_rufe(); else e_sufe();
                    }
                 else if (*flag==o_text[31][0])
                    { 
		    if (e_ie_e()) e_riee(); else e_siee();
                    }
                 flag++;
                 }
              }

#ifdef OPTION_ENABLE
        /* display info file                             */
        else if (strncmp(option,o_text[16],strlen(o_text[16]))==0)
              {
              f_text desc;
              s_trng path;

              desc.text = TRUE;
              desc.ellen = 1;

              path.suba = TRUE;
              path.fix = path.tmp = FALSE;

              key = option+strlen(o_text[16]);
              if ((flag = strchr(key,(int)'@'))!=NULL)
                 {
                 size = flag-key;
                 path.ptr = flag+1;
                 path.clen = strlen(path.ptr);
                 }
              else
                 {
                 size = strlen(key);
                 path.ptr = o_text[7];
                 path.clen = strlen(o_text[7]);
                 }
              path.alen = path.clen+1;

              if (b_op88(&desc,path,100)==0)
                 {
                 fprintf(device,
                         "Information file \"%s\" not available.\n",
                         path.ptr);
                 }
              else
                 {

                 /* introduction and keywords displayed if no key is given */
                 if (*key==':') 
                    {
                    key++;
                    size--;
                    }
                 key_found = (size) ? FALSE : TRUE;

                 /* count of keywords */
                 i = 0;

                 /* scan characters of info file */
                 ch = fgetc(desc.fp);
                 while (ch!=EOF)
                    {

                    /* display character */
                    if (key_found) fprintf(device,"%c",ch);

                    /* if newline then check for new keyword */
                    if (ch=='\n')
                       {

                       /* end of keyword display after introduction */
                       if (size==0 && i>0) key_found = FALSE;

                       /* new keyword starts */
                       if ((ch = fgetc(desc.fp))==':')
                          {
                          i++;

                          /* display keyword only */
                          if (size==0) key_found = TRUE;

                          /* check for specified keyword */
                          else
                             {

                             /* only first 'size' characters must match */
                             k = 0;
                             while ((ch = fgetc(desc.fp))==key[k]) k++;

                             /* keyword found */
                             key_found =
                                (key[k]=='\0' || key[k]=='@')
                                            ? TRUE : FALSE;

                             /* skip info file until end of line */
                             while (ch!='\n' && ch!=EOF)
                                ch = fgetc(desc.fp);

                             if (ch!=EOF) ch = fgetc(desc.fp);
                             }
                          }
                       }
                    else ch = fgetc(desc.fp);
                    }

                 /* close info file     */
                 fclose(desc.fp);
                 }
              return(2);
              }
#endif

        else
              {
              return(0);
              }

        return(1);
        }





