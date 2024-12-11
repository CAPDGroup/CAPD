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

/* CVS $Id: dot_defs.hpp,v 1.25 2014/01/30 17:23:45 cxsc Exp $ */

/**********************************************************************
 * (C) 1993 University of Karlsruhe
 **********************************************************************/
#ifndef _CXSC_DOT_DEFS_HPP_INCLUDED
#define _CXSC_DOT_DEFS_HPP_INCLUDED

namespace cxsc {

/* ----------------------------------------------------------------- */
/* Maximale dezimale Laenge des ganzzahligen Anteils eines Akkus ist */
/*    log(2)/log(10) * maximale binaere Laenge                       */

#define A_I_DIGITS  ((a_intg)((A_I_LENGTH * B_LENGTH * 30103L) / 100000L))

/* ---------------------------------------------------------------- */
/* Maximale dezimale Laenge des gebrochenen Anteils eines Akkus ist */
/*     maximale binaere Laenge                                      */

#define A_F_DIGITS  ((A_F_LENGTH * B_LENGTH) + 1)

/* ---------------------------------------------------------------- */
/* Anzahl der Digits fuer den Exponenten                            */

#define A_E_DIGITS  4
#define A_E_MAX     10000

/* ---------------------------------------------------------------- */
/* Laenge eines Strings zur Aufnahme eines Akku                     */

#define A_DIGITS   (1 + A_I_DIGITS + 1 + A_F_DIGITS + 2 + A_E_DIGITS + 1 + 20)

#if _WIN32
extern __declspec(thread) char *dm;
extern __declspec(thread) char *dmhlp;
#elif __APPLE__ && !CXSC_FORCE_TLS
extern char *dm;
extern char *dmhlp;
#else
extern __thread char *dm;
extern __thread char *dmhlp;
#endif


int d_init_dm (void);

/* ---------------------------------------------------------------- */

  /*---- D_OUT.C -- mr 30.09.1990 ----------------------------------*/
  void d_out(a_intg*, char*, a_intg*, a_intg*, Dotprecision);

  /*---- D_OUTP.C - mr 30.09.1990 ----------------------------------*/
  void d_outp(char*, Dotprecision, a_intg, a_intg, a_intg, a_intg*);

  /*---- D_SCAN.C - mr 19.10.1990 ----------------------------------*/
  char* d_scan (char*, a_intg*, a_intg*, char*, a_intg*, a_intg*);

  /*---- D_SCANI.C - mr 19.10.1990 ---------------------------------*/
  void d_scani(Dotprecision, char*, a_intg*, a_intg*, a_intg*);

  /*---- D_SCANF.C - mr 19.10.1990 ---------------------------------*/
  a_intg d_scanf(Dotprecision, char*, a_intg*, a_intg*, a_intg*, a_intg);

  /*---- D_SCANP.C - mr 20.10.1990 --------------------------------*/
  char* d_scanp(Dotprecision, char*, a_intg, a_intg*);

/* ---------------------------------------------------------------- */

} // namespace cxsc 

#endif // _CXSC_DOT_DEFS_HPP_INCLUDED
