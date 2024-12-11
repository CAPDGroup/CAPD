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

/* CVS $Id: RtsFunc.h,v 1.22 2014/01/30 17:23:43 cxsc Exp $ */

#ifndef _CXSC_RTSFUNC_H_INCLUDED
#define _CXSC_RTSFUNC_H_INCLUDED

/* Aus dotlib.hpp */

extern "C" {

/*extern char *dm;
extern char *dmhlp;


int d_init_dm (void);
*/
#define dotprecision Dotprecision

/* ---------------------------------------------------------------- */

  /*---- D_OUT.C -- mr 30.09.1990 ----------------------------------*/
  void d_out(a_intg*, char*, a_intg*, a_intg*, dotprecision);

  /*---- D_OUTP.C - mr 30.09.1990 ----------------------------------*/
  void d_outp(char*, dotprecision, a_intg, a_intg, a_intg, a_intg*);

  /*---- D_SCAN.C - mr 19.10.1990 ----------------------------------*/
  char* d_scan (char*, a_intg*, a_intg*, char*, a_intg*, a_intg*);

  /*---- D_SCANI.C - mr 19.10.1990 ---------------------------------*/
  void d_scani(dotprecision, char*, a_intg*, a_intg*, a_intg*);

  /*---- D_SCANF.C - mr 19.10.1990 ---------------------------------*/
  a_intg d_scanf(dotprecision, char*, a_intg*, a_intg*, a_intg*, a_intg);

  /*---- D_SCANP.C - mr 20.10.1990 --------------------------------*/
  char* d_scanp(dotprecision, char*, a_intg, a_intg*);

/* ---------------------------------------------------------------- */

// Fuer dotio:
void     b_rnd (a_intg rnd,char *buffer,a_intg digits,a_intg pos,a_intg *bdp,a_intg *dexpo); 
void     b_outi(a_intg *digits,char *buffer,a_intg *bdp,a_intg *dexpo,dotprecision c);      
void     b_outf(a_intg *digits,char *buffer,a_intg *bdp,a_intg *dexpo,dotprecision c);  

// Fuer realio:
a_bool   b_deko(a_real x,a_intg *expo,a_btyp *mant,a_bool *vz);
void     b_out (a_btyp *mant,a_intg expo,a_intg digits,char *buffer,a_intg *bdp,a_intg *dexpo);
void     b_rnd (a_intg rnd,char *buffer,a_intg digits,a_intg pos,a_intg *bdp,a_intg *dexpo);
                
// Fuer l_interv.cpp
void     b_shr1(a_btyp *lang,a_intg laenge);



#undef dotprecision
}

#endif // _CXSC_RTSFUNC_H_INCLUDED
