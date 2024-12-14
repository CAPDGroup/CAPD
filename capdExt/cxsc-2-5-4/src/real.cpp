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

/* CVS $Id: real.cpp,v 1.28 2014/01/30 17:23:48 cxsc Exp $ */

#include "real.hpp"
#include "ioflags.hpp"

namespace cxsc {

//----------------------------------------------------------------------------
// MakeHexReal - erstellt aus den binaer angegebenen Einzelteilen einer
//               IEEE 64-bit Gleitkommazahl eine solche, hierbei ist
//                 sign         Flag 'zu erzeugener real ist negativ'
//                 expo         Exponent (11 Bit 0x000 - 0x7FF)
//                 manthigh     obere  20 Bit der Mantisse (Bit 32 - 53)
//                 mantlow      untere 32 Bit der Mantisse (Bit  0 - 31)
//
//             - die Mantisse ist normalisiert dargestellt, d.h.
//               implizit wird noch ein Bit 1 als Vorkommastelle
//               vorangestellt.
//               Eine Ausnahme bilden die Zahlen mit Exponent = 0,
//               oder Exponent = 0x7FF (alle Exponentenbits = 1)
//
/*!
\param sign The sign of the number
\param expo The exponent of the number
\param manthigh The high byte of the mantissa
\param mantlow The low byte of the mantissa
\return The IEEE 64-bit floating point number
*/
const real& MakeHexReal(int sign, unsigned int expo, a_btyp manthigh, a_btyp mantlow)
{
  static real a;
  ((a_btyp*)&a)[LOWREAL]   = mantlow,
  ((a_btyp*)&a)[HIGHREAL]  = manthigh & 0xFFFFFL,
  ((a_btyp*)&a)[HIGHREAL] |= ((a_btyp) (expo & 0x7FF)) << 20,
  ((a_btyp*)&a)[HIGHREAL] |= (sign ? 0x80000000L : 0x00000000);
  return a;
}

const real MinReal      = MakeHexReal(0, 0x001, 0x00000L, 0x00000000L);
const real minreal      = MakeHexReal(0, 0x000, 0x00000L, 0x00000001L);
// Blomquist, 26.09.02; minreal = smallest positive real number.
const real MaxReal      = MakeHexReal(0, 0x7FE, 0xFFFFFL, 0xFFFFFFFFL);
const real Infinity     = MakeHexReal(0, 0x7FF, 0x00000L, 0x00000000L);
const real SignalingNaN = MakeHexReal(1, 0x7FF, 0x80000L, 0x00000000L);
const real QuietNaN     = MakeHexReal(0, 0x7FF, 0x00000L, 0x00000001L);
const real Epsilon      = power(2,-53);
const real Factor       = power(2, 27) + 1;


// The following constants are roundet to the nearest machine nunmber:

const real Pi_real = 7074237752028440.0 / 2251799813685248.0;  // Pi
const real Pi2_real = 7074237752028440.0/1125899906842624.0;   // 2*Pi
const real Pi3_real = 5305678314021330.0/562949953421312.0;    // 3*Pi
const real Pid2_real = 7074237752028440.0/4503599627370496.0;  // Pi/2
const real Pid3_real = 4716158501352294.0/4503599627370496.0;  // Pi/3
const real Pid4_real = 7074237752028440.0/9007199254740992.0;  // Pi/4
const real Pir_real = 5734161139222659.0/18014398509481984.0;  // 1/Pi
const real Pi2r_real = 5734161139222659.0/36028797018963968.0; // 1/(2*Pi)
const real Pip2_real = 5556093337880030.0/562949953421312.0;   // Pi^2
const real SqrtPi_real = 7982422502469483.0/4503599627370496.0;// sqrt(Pi)
const real Sqrt2Pi_real = 5644425081792262.0/2251799813685248.0; // sqrt(2Pi)
const real SqrtPir_real = 5081767996463981.0/9007199254740992.0; // 1/sqrt(Pi)
const real Sqrt2Pir_real = 7186705221432913.0/18014398509481984.0; // 1/sqrt(2Pi)
const real Sqrt2_real = 6369051672525773.0/4503599627370496.0; // sqrt(2)
const real Sqrt5_real = 5035177455121576.0 / 2251799813685248.0; // sqrt(5)
const real Sqrt7_real = 5957702309312746.0 / 2251799813685248.0; // sqrt(7)
const real Sqrt2r_real = 6369051672525773.0/9007199254740992.0;// 1/sqrt(2)
const real Sqrt3_real = 7800463371553962.0/4503599627370496.0; // sqrt(3)
const real Sqrt3d2_real = 7800463371553962.0/9007199254740992.0;  // sqrt(3)/2
const real Sqrt3r_real = 5200308914369308.0/9007199254740992.0;// 1/sqrt(3)
const real Ln2_real = 6243314768165359.0 / 9007199254740992.0; // ln(2)
const real Ln2r_real = 6497320848556798.0 / 4503599627370496.0; // 1/ln(2)
const real Ln10_real = 5184960683398422.0 / 2251799813685248.0; // ln(10)
const real Ln10r_real = 7823553867474190.0/18014398509481984.0; // 1/ln(10)
const real LnPi_real = 5155405087351229.0 / 4503599627370496.0; // ln(Pi)
const real Ln2Pi_real = 8277062471433909.0/4503599627370496.0;  // ln(2Pi)
const real E_real = 6121026514868073.0 / 2251799813685248.0;    // e
const real Er_real = 6627126856707896.0 / 18014398509481984.0;  // 1/e
const real Ep2_real = 8319337573440942.0 / 1125899906842624.0;  // e^2
const real Ep2r_real = 4875967449235916.0/36028797018963968.0;  // 1/e^2
const real EpPi_real = 6513525919879994.0/281474976710656.0;    // e^(Pi)
const real Ep2Pi_real = 4710234414611993.0/8796093022208.0;     // e^(2Pi)
const real EpPid2_real = 5416116035097439.0/1125899906842624.0; // e^(Pi/2)
const real EpPid4_real = 4938827609611434.0/2251799813685248.0; // e^(Pi/4)

} // namespace cxsc
