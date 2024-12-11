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

/* CVS $Id: rts_real.hpp,v 1.24 2014/01/30 17:23:48 cxsc Exp $ */


#ifndef _CXSC_RTS_REAL_HPP_INCLUDED
#define _CXSC_RTS_REAL_HPP_INCLUDED

/*	Deklarationen der Funktionen des Laufzeitsystems, die fuer die
	gerichtete Rundung der Klasse "real" benoetigt werden, und
	Definition dazu benoetigter Konvertierungsoperatoren.
	(Declaration of functions of the runtime system, used for the
	directed rounding in the class "real", and definitions of needed
	conversion operators.)
*/

// Deklaration von Hilfsfunktionen des Laufzeitsystems
// Der Typ a_real ist in p88rts.h des Laufzeitsystems als 64-Bit-Zahl definiert
// (Help functions of the runtime system declarated here)

#include "RtsTyp.h"

namespace cxsc {

// Verknuepfungen mit Rundung nach oben bzw. unten
// (Operators with rounding upwards or downwards)
extern "C" {
#if ( OPT80387 )
	a_real r_addd80387 (a_real a, a_real b); // speziell fuer 80387-Koprozessor
	a_real r_addu80387 (a_real a, a_real b); // (especially for 80387)
	a_real r_subd80387 (a_real a, a_real b);
	a_real r_subu80387 (a_real a, a_real b);
	a_real r_muld80387 (a_real a, a_real b);
	a_real r_mulu80387 (a_real a, a_real b);
	a_real r_divd80387 (a_real a, a_real b);
	a_real r_divu80387 (a_real a, a_real b);
#else
	a_real r_addd (a_real a, a_real b);
	a_real r_addu (a_real a, a_real b);
	a_real r_subd (a_real a, a_real b);
	a_real r_subu (a_real a, a_real b);
	a_real r_muld (a_real a, a_real b);
	a_real r_mulu (a_real a, a_real b);
	a_real r_divd (a_real a, a_real b);
	a_real r_divu (a_real a, a_real b);
#endif

#if HP_9000_CPP+SUN4_CPP_C
	void r_lfsr();
#endif
}

// namespace real
//{

// Operatoren fuer Umwandlungen zwischen real und a_real
// (operators for conversions between real and a_real)

class real;

inline a_real _a_real(const real &x)
{	return *((const a_real *)(&x));
}

// As a_real are doubles this is defined in real.hpp
// inline real _real(const a_real &x)
// {	return *((const real *)(&x));
// }

//} // namespace real

} // namespace cxsc 

#endif // _CXSC_RTS_REAL_HPP_INCLUDED

