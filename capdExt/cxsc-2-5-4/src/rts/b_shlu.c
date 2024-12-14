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

/* CVS $Id: b_shlu.c,v 1.21 2014/01/30 17:24:04 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_shlu.c                              */
/*                                                              */
/*      Entries         : void b_shlu(lang,laenge,dist);        */
/*                        a_btyp *lang;                         */
/*                        a_intg laenge,dist;                   */
/*                                                              */
/*      Arguments       : lang = array of length laenge         */
/*                        laenge = length of array lang         */
/*                        dist = number of bits to be shifted   */
/*                                                              */
/*      Description     : Shift a_btyp array lang of length     */
/*                        laenge by dist positions to left      */
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
local void b_shlu(a_btyp *lang,a_intg laenge,a_intg dist)
#else
local void b_shlu(lang,laenge,dist)

a_btyp *lang;
a_intg laenge;
a_intg dist;
#endif
        {
        a_intg i,bigdist;

        /* shift by blocks of B_LENGTH                          */
        if ( (bigdist = B_ASHR(dist,LOG_B_LENGTH))!=0 ) {
                for ( i=0; i<laenge-bigdist; i++ )
                        lang[i] = lang[i+bigdist];
                for ( i=laenge-1; i>=laenge-bigdist && i>=0; i-- )
                    lang[i] = ZERO;
                dist &= (B_LENGTH-1);
                }

        /* shift by block of length dist                        */
        if (dist) {
            for ( i=0; i<laenge-bigdist-1; i++ )
                lang[i] = (lang[i]<<dist) |
                          (lang[i+1]>>(B_LENGTH-dist));
            lang[i] <<= dist;
            }

        return;
        }





