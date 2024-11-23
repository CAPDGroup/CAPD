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

/* CVS $Id: RtsTyp.h,v 1.22 2014/01/30 17:23:43 cxsc Exp $ */

#ifndef _CXSC_RTSTYP_H_INCLUDED
#define _CXSC_RTSTYP_H_INCLUDED
/* Dieses Headerfile definiert die grundsaetzlichen Datentypen des RTS

Die Header-Dateien des RTS muessen im Include-Pfad liegen!

*/

#include <stddef.h> /* Fuer Definition von "size_t" in p88rts */
#include <stdio.h>  /* Fuer Definition von "FILE *" in p88rts ... (?!? Wer braucht das ?!?)*/

extern "C" { /* Da die folgenden Header zu einer C-Libary gehoeren */

#include <p88rts.h>

#define dotprecision Dotprecision
#include <r_defs.h> /* fuer b_defs */
#include <b_defs.h> /* fuer d_defs */
#include <d_defs.h> /* fuer Buffersize eines Dotprecision */
#include <o_type.h> /* fuer TByte */
#undef  dotprecision

}
#endif // _CXSC_RTSTYP_H_INCLUDED
