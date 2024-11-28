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

/* CVS $Id: oldinterval.cpp,v 1.23 2014/01/30 17:23:48 cxsc Exp $ */

//----------------------------------------------------------------------------
// File "Interval.CPP"  vom 22.05.1990
//----------------------------------------------------------------------------
//
// Realisierung des skalaren Rechentyps interval
//
//----------------------------------------------------------------------------

#include "dotlib.hpp"

#include "errors.hpp"
#include "interval.hpp"

namespace cxsc {

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------

#if (OPT80387)

#define RCDown  0x0400
#define RCUp    0x0800
#define RCNear  0x0000
#define RCChop  0x0C00

#define RCMask  (RCDown OR RCUp OR RCNear OR RCChop)

  // -------------------------------------------------------------
  //    fnclex                  // clear exceptions, no check
  //    fstcw   bx              // store to ctrlwrd
  //    and     bx,RCMask       // clear RC
  //    or      bx,mode \
  //    fldcw   bx              // load changed ctrlwrd

#define setrndmode(mode) \
    asm { \
          fnclex ; \
          fstcw word ptr ctrlwrd ; \
          mov   bx, word ptr ctrlwrd ; \
          and   bx, NOT RCMask ; \
          or    bx, mode ; \
          mov   word ptr ctrlwrd, bx ; \
          fldcw word ptr ctrlwrd \
    }

  // -------------------------------------------------------------
  //    fnclex                  // clear exceptions, no check
  //    fstcw   save            // store to ctrlwrd
  //    and     save, RCMask    // clear all but RC

#define getrndmode(save) \
    asm { \
      fnclex ;               \
      fstcw word ptr save ;  \
      and word ptr save, RCMask   \
    }


// double name##80387 (real& a, real& b) { \
// erg = _double(a) operator _double(b); \

#define Co80387(name,rndmode,operator) \
    double name##80387 (double a, double b) { \
      int rnd, ctrlwrd; \
      double erg; \
      getrndmode (rnd); \
      setrndmode (rndmode); \
      erg = a operator b; \
      setrndmode (word ptr rnd); \
      return erg; \
    }

  Co80387 (r_addd, RCDown, +)
  Co80387 (r_addu, RCUp,   +)
  Co80387 (r_subd, RCDown, -)
  Co80387 (r_subu, RCUp,   -)
  Co80387 (r_muld, RCDown, *)
  Co80387 (r_mulu, RCUp,   *)
  Co80387 (r_divd, RCDown, /)
  Co80387 (r_divu, RCUp,   /)
#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Explizite Typumwandlungsfunktionen

interval _interval(const real& a) {
  return interval (a,a);
}
interval _interval(const real& a, const real& b) {
  return interval (a,b);
}

interval _unchecked_interval(const real& a, const real& b) {
  interval tmp;
  tmp.inf = a;
  tmp.sup = b;
  return tmp;
}

//----------------------------------------------------------------------------
// Konstruktoren fuer complex's
// Zuweisungsoperatoren

  interval::interval(const real& a, const real& b) {
    if (a > b) errmon (ERR_INTERVAL(EMPTY));
    ((a_btyp*)&(((double*)this)[0]))[LOWREAL]  = ((a_btyp*)&a)[LOWREAL],
    ((a_btyp*)&(((double*)this)[0]))[HIGHREAL] = ((a_btyp*)&a)[HIGHREAL];
    ((a_btyp*)&(((double*)this)[1]))[LOWREAL]  = ((a_btyp*)&b)[LOWREAL],
    ((a_btyp*)&(((double*)this)[1]))[HIGHREAL] = ((a_btyp*)&b)[HIGHREAL];
  }
  interval::interval(const interval& a) {
    ((a_btyp*)&(((double*)this)[0]))[LOWREAL]  = ((a_btyp*)&a.inf)[LOWREAL],
    ((a_btyp*)&(((double*)this)[0]))[HIGHREAL] = ((a_btyp*)&a.inf)[HIGHREAL];
    ((a_btyp*)&(((double*)this)[1]))[LOWREAL]  = ((a_btyp*)&a.sup)[LOWREAL],
    ((a_btyp*)&(((double*)this)[1]))[HIGHREAL] = ((a_btyp*)&a.sup)[HIGHREAL];
  }

  interval& interval::operator= (const real& a) {
    this->inf = this->sup = a;
    return *this;
  }
  interval& interval::operator= (interval& a) {
    this->inf = a.inf;
    this->sup = a.sup;
    return *this;
  }

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

  //-----------------------------------------------------
  // Eingabeoperatoren auf Streams und Strings

  char*    operator>> (char *s,    interval& a)
  {
    s = skipwhitespacessinglechar (s, '[');
    s = s >> SaveOpt >> RndDown >> a.inf;
    s = skipwhitespacessinglechar (s, ',');
    s = s >> RndUp >> a.sup >> RestoreOpt;
    s = skipwhitespaces (s);
    if (*s == ']') s++;

    if (a.inf > a.sup) {
      errmon (ERR_INTERVAL(EMPTY));
    }
    return s;
  }

  //-----------------------------------------------------

  std::istream& operator>> (std::istream& s, interval& a)
  {
    char c;

    skipeolnflag = inpdotflag = TRUE;
    c = skipwhitespacessinglechar (s, '[');
    if (inpdotflag) s.putback(c);

    s >> SaveOpt >> RndDown >> a.inf;

    skipeolnflag = inpdotflag = TRUE;
    c = skipwhitespacessinglechar (s, ',');
    if (inpdotflag) s.putback(c);

    s >> RndUp >> a.sup >> RestoreOpt;

    if (!waseolnflag) {
      skipeolnflag = FALSE, inpdotflag = TRUE;
      c = skipwhitespaces (s);
      if (inpdotflag && c != ']') s.putback(c);
    }

    if (a.inf > a.sup) {
      errmon (ERR_INTERVAL(EMPTY));
    }
    return s;
  }

  //-----------------------------------------------------
  // Ausgabeoperatoren auf Streams und Strings

  char* operator<< (char* s, interval& a)
  {
    sprintf (s, "[");
    sprintf (&s[strlen(s)] << SaveOpt << RndDown << a.inf, ",");
    sprintf (&s[strlen(s)] << RndUp << a.sup << RestoreOpt, "]");
    return &s[strlen(s)];
  }

  //----------------------------------------------------------

  std::ostream& operator<< (std::ostream& s, interval& a)
  {
    s << '[' << SaveOpt << RndDown << a.inf << ','
      << RndUp << a.sup << RestoreOpt << ']';
    return s;
  }

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

  //----------------------------------------------------------
  // Operatoren : interval 'operator' interval
  // Standard Rechen-Operatoren fuer interval's

  interval operator-  (interval& a) {
    return interval (-a.sup, -a.inf);
  }
  interval operator+  (interval& a, interval& b) {
    return interval (
      addd(a.inf,b.inf),
      addu(a.sup,b.sup));
  }
  interval operator-  (interval& a, interval& b) {
    return interval (
      subd(a.inf,b.sup),
      subu(a.sup,b.inf));
  }

  //------------------------------------------------------------------------
        //  a * b  = ueber Entscheidungstabelle :
  //
  //                      bi,bs >= 0       bi < 0, bs >=0       bi,bs < 0
  // ----------------+------------------+------------------+----------------
  // ai,as >= 0      I   ai*bi, as*bs   I   as*bi, as*bs   I   as*bi, ai*bs
  // ----------------+------------------+------------------+----------------
  // ai < 0, as >= 0 I   ai*bs, as*bs   I      ....        I   as*bi, ai*bi
  // ----------------+------------------+------------------+----------------
  // ai,as < 0       I   ai*bs, as*bi   I   ai*bs, ai*bi   I   as*bs, ai*bi
  // ----------------+------------------+------------------+----------------
  //
  //  .... :  min(ai*bs, as*bi), max(ai*ai, as*as)
  //
  interval operator*  (interval& a, interval& b)
  {
    interval tmp;

    if (sign(a.inf) >= 0) {             // 1. Zeile der Entscheidungstabelle
      if (sign(b.inf) >= 0)
      {                         //     1. Spalte: [ai*bi, as*bs]
        tmp.inf = muld(a.inf,b.inf);
        tmp.sup = mulu(a.sup,b.sup);
      }
      else if (sign(b.sup) >= 0) {    //     2. Spalte: [as*bi, as*bs]
        tmp.inf = muld(a.sup,b.inf);
        tmp.sup = mulu(a.sup,b.sup);
      }
      else {                            //     3. Spalte: [as*bi, ai*bs]
        tmp.inf = muld(a.sup,b.inf);
        tmp.sup = mulu(a.inf,b.sup);
      }
    }
    else if (sign(a.sup) >= 0) {      // 2. Zeile der Entscheidungstabelle
      if (sign(b.inf) >= 0) {           //     1. Spalte: [ai*bs, as*bs]
        tmp.inf = muld(a.inf,b.sup);
        tmp.sup = mulu(a.sup,b.sup);
      }
      else if (sign(b.sup) >= 0) {    //   2. Spalte: [min(ai*bs, as*bi),
        real hlp;                       //                 max(ai*ai, as*as)]

        tmp.inf = muld(a.inf,b.sup);
        hlp     = muld(a.sup,b.inf);
        if (hlp < tmp.inf) tmp.inf = hlp;

        tmp.sup = mulu(a.inf,b.inf);
        hlp     = mulu(a.sup,b.sup);
        if (hlp > tmp.sup) tmp.sup = hlp;  //hier stand bis 30.04.91 tmp.inf!
      }
      else {                            //     3. Spalte: [as*bi, ai*bi]
        tmp.inf = muld(a.sup,b.inf);
        tmp.sup = mulu(a.inf,b.inf);
      }
    }
    else {                              // 3. Zeile der Entscheidungstabelle
      if (sign(b.inf) >= 0) {           //     1. Spalte: [ai*bs, as*bi]
        tmp.inf = muld(a.inf,b.sup);
        tmp.sup = mulu(a.sup,b.inf);
      }
      else if (sign(b.sup) >= 0) {    //   2. Spalte: [ai*bs, ai*bi]
        tmp.inf = muld(a.inf,b.sup);
        tmp.sup = mulu(a.inf,b.inf);
      }
      else {                            //     3. Spalte: [as*bs, ai*bi]
        tmp.inf = muld(a.sup,b.sup);
        tmp.sup = mulu(a.inf,b.inf);
      }
    }

    return tmp;
  }

  //------------------------------------------------------------------------
        //  a / b  = ueber Entscheidungstabelle :
  //
  //                      bi,bs > 0          bi,bs < 0
  // ----------------+------------------+------------------
  // ai,as >= 0      I   ai/bs, as/bi   I   as/bs, ai/bi
  // ----------------+------------------+------------------
  // ai < 0, as >= 0 I   ai/bi, as/bi   I   as/bs, ai/bs
  // ----------------+------------------+------------------
  // ai,as < 0       I   ai/bi, as/bs   I   as/bi, ai/bs
  // ----------------+------------------+------------------
  //
  interval operator/  (interval& a, interval& b)
  {
    interval tmp;

    if (!b) {
      errmon (ERR_INTERVAL(DIVZERO));
      tmp.inf = (!a) ? QuietNaN : Infinity;  // mr????
      tmp.sup = tmp.inf;
      return tmp;
    }

    if (sign(a.inf) >= 0) {     // 1. Zeile der Entscheidungstabelle
      if (sign(b.inf) > 0) {            //     1. Spalte: [ai/bs, as/bi]
        tmp.inf = divd(a.inf,b.sup);
        tmp.sup = divu(a.sup,b.inf);
      }
      else {                            //     2. Spalte: [as/bs, ai/bi]
        tmp.inf = divd(a.sup,b.sup);
        tmp.sup = divu(a.inf,b.inf);
      }
    }
    else if (sign(a.sup) >= 0) {      // 2. Zeile der Entscheidungstabelle
      if (sign(b.inf) > 0) {          //     1. Spalte: [ai/bi, as/bi]
        tmp.inf = divd(a.inf,b.inf);
        tmp.sup = divu(a.sup,b.inf);
      }
      else {                            //     2. Spalte: [as/bs, ai/bs]
        tmp.inf = divd(a.sup,b.sup);
        tmp.sup = divu(a.inf,b.sup);
      }
    }
    else {                              // 3. Zeile der Entscheidungstabelle
      if (sign(b.inf) > 0) {          //     1. Spalte: [ai/bi, as/bs]
        tmp.inf = divd(a.inf,b.inf);
        tmp.sup = divu(a.sup,b.sup);
      }
      else {                            //     2. Spalte: [as/bi, ai/bs]
        tmp.inf = divd(a.sup,b.inf);
        tmp.sup = divu(a.inf,b.sup);
      }
    }

    return tmp;
  }

  //----------------------------------------------------------

  interval& operator+= (interval& a, interval& b) { return a = a + b; }
  interval& operator-= (interval& a, interval& b) { return a = a - b; }
  interval& operator*= (interval& a, interval& b) { return a = a * b; }
  interval& operator/= (interval& a, interval& b) { return a = a / b; }

  // -----------------------------------------------------
  // Standard Vergleichsoperatoren fuer interval's

  int operator!  (interval& a) {
    return (a.inf <= CXSC_Zero && a.sup >= CXSC_Zero);
  }
  int operator== (interval& a, interval& b) {
    return (a.inf == b.inf && a.sup == b.sup);
  }
  int operator!= (interval& a, interval& b) {
    return (a.inf != b.inf || a.sup != b.sup);
  }

  // -----------------------------------------------------
  // Standard Mengenvergleichsoperatoren fuer interval's

  int operator<  (interval& a, interval& b) {
    return (a != b) && (a <= b);
  }
  int operator<= (interval& a, interval& b) {
    return (a.inf >= b.inf && a.sup <= b.sup);
  }
  int operator>= (interval& a, interval& b) {
    return (b <= a);
  }
  int operator>  (interval& a, interval& b) {
    return (a != b) && (b <= a);
  }

  // -----------------------------------------------------
  // Standard Mengenoperatoren fuer interval's

  //------------------------------------------------------------------
  // Schnitt zweier interval:
  //    falls der Schnitt leer ist:
        //       ==> Rueckgabe des 'leeren' interval [1,-1]
  //           Bem. Im Konstruktor interval() wird eine Fehlermeldung
  //                generiert
  //    sonst
  //       ==> interval (max(a.inf,b.inf), min (a.sup, b.sup));
  //
  interval operator& (interval& a, interval& b) {
    if (a.sup < b.inf || a.inf > b.sup) {
      errmon(ERR_INTERVAL(EMPTY));
      return _unchecked_interval(1.0, -1.0);
    }
    if (a.inf > b.inf) {
      if (a.sup < b.sup) return a;                      // [a.inf, a.sup]
      return interval (a.inf, b.sup);                   // [a.inf, b.sup]
    }
    else {
      if (a.sup < b.sup) {
        return interval (b.inf, a.sup);                 // [b.inf, a.sup]
      }
      return b;                                         // [b.inf, b.sup]
    }
  }

  //----------------------------------------------------------------------
        // Konvexe Huelle:
  //    interval (min(a.inf,b.inf), max (a.sup, b.sup));
  //
  interval operator| (interval& a, interval& b) {
    if (a.inf > b.inf) {
      if (a.sup > b.sup) {
        return interval (b.inf, a.sup);                 // [b.inf, a.sup]
      }
      return b;                                         // [b.inf, b.sup]
    }
    else {
      if (a.sup > b.sup) return a;                      // [a.inf, a.sup]
      return interval (a.inf, b.sup);                   // [a.inf, b.sup]
    }
  }

  interval& operator&= (interval& a, interval& b) { return a = a & b; }
  interval& operator|= (interval& a, interval& b) { return a = a | b; }


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  // Operatoren : interval 'operator' real, bzw. real 'operator' interval
  // Standard Rechen-Operatoren fuer interval's

  interval operator+  (const real& a, interval& b) {
    return interval (
      addd(a,b.inf),
      addu(a,b.sup));
  }
  interval operator+  (interval& a, const real& b) {
    return interval (
       addd(a.inf,b),
       addu(a.sup,b));
  }
  interval operator-  (const real& a, interval& b) {
    return interval (
      subd(a,b.sup),
      subu(a,b.inf));
  }
  interval operator-  (interval& a, const real& b) {
    return interval (
      subd(a.inf,b),
      subu(a.sup,b));
  }

  interval operator*  (const real& a, interval& b)
  {
    interval tmp;

    if (sign(a) == 0) tmp = _interval(0.0,0.0);
    else if (sign(a) > 0) {
      tmp.inf = muld(a,b.inf);
      tmp.sup = mulu(a,b.sup);
    }
    else if (sign(a) < 0) {
      tmp.inf = muld(a,b.sup);
      tmp.sup = mulu(a,b.inf);
    }
    return tmp;
  }

  interval operator*  (interval& a, const real& b)
  {
    interval tmp;

    if (sign(b) == 0) tmp = _interval(0.0,0.0);
    else if (sign(b) > 0) {
      tmp.inf = muld(a.inf,b);
      tmp.sup = mulu(a.sup,b);
    }
    else if (sign(b) < 0) {
      tmp.inf = muld(a.sup,b);
      tmp.sup = mulu(a.inf,b);
    }
    return tmp;
  }

  interval operator/  (const real& a, interval& b) {
    return (_interval(a) / b);
  }

  interval operator/  (interval& a, const real& b) {
    return (a / _interval(b));
  }

  interval& operator+= (interval& a, const real& b) { return a = a + b; }
  interval& operator-= (interval& a, const real& b) { return a = a - b; }
  interval& operator*= (interval& a, const real& b) { return a = a * b; }
  interval& operator/= (interval& a, const real& b) { return a = a / b; }


  // -----------------------------------------------------
  // Standard Vergleichsoperatoren fuer interval's

  int operator==  (const real& a, interval& b) {
    return interval(a,a) == b;
  }
  int operator==  (interval& a, const real& b) {
    return a == interval(b,b);
  }
  int operator!=  (const real& a, interval& b) {
    return interval(a,a) != b;
  }
  int operator!=  (interval& a, const real& b) {
    return a != interval(b,b);
  }

  // --------------------------------------------------------------
  // Standard Mengenvergleichsoperatoren fuer interval's mit double's

  int operator<  (const real& a, interval& b) {
    return interval(a,a) <  b;
  }
  int operator<= (const real& a, interval& b) {
    return interval(a,a) <= b;
  }
  int operator<= (interval &a, const real& b) {
    return a <= interval(b,b);
  }
  int operator>= (const real& a, interval& b) {
    return interval(a,a) >= b;
  }
  int operator>= (interval &a, const real& b) {
    return a >= interval(b,b);
  }
  int operator>  (interval &a, const real& b) {
    return a > interval(b,b);
  }

  // --------------------------------------------------------------
  // Standard Mengenoperatoren fuer interval's mit double's

  interval operator& (const real& a, interval& b) {
    return interval(a,a) & b;
  }
  interval operator& (interval &a, const real& b) {
    return a & interval(b,b);
  }
  interval operator| (const real& a, interval& b) {
    return interval(a,a) | b;
  }
  interval operator| (interval &a, const real& b) {
    return a | interval(b,b);
  }
  interval& operator&= (interval& a, const real& b) {
    return a = a & interval(b,b);
  }
  interval& operator|= (interval& a, const real& b) {
    return a = a | interval(b,b);
  }

  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
  // Diverse Funktionen

  real& Inf (interval& a) { return a.inf; }
  real& Sup (interval& a) { return a.sup; }

  interval& SetInf (interval& a, const real& b) {
    if (b > a.sup) {
      errmon(ERR_INTERVAL(EMPTY));
    }
    a.inf = b; return a;
  }
  interval& SetSup (interval& a, const real& b) {
    if (a.inf > b) {
      errmon(ERR_INTERVAL(EMPTY));
    }
    a.sup = b; return a;
  }
  interval& UncheckedSetInf (interval& a, const real& b) {
    a.inf = b; return a;
  }
  interval& UncheckedSetSup (interval& a, const real& b) {
    a.sup = b; return a;
  }

  int IsEmpty (interval& a)   { return (a.inf  > a.sup); }

  interval abs  (interval& a)
  {
    real h1  = abs(a.inf);
    real h2  = abs(a.sup);

    if (IsEmpty(a)) return a;
    if (!a)         return interval(0.0, (h1 > h2) ? h1 : h2);
    if (h1 > h2)    return interval(h2, h1);
    return interval(h1, h2);
  }

  real     mid  (interval& a)
  {
    // --------------------------------------------------------
    //   dotakku[4] = a.inf + a.sup

    dotprecision dot(a.inf);
    dot += a.sup;

    if (dot != 0.0) {
     // Division nur bei ungleich 0.0
      Dotprecision ptr = *dot.ptr();

      // --------------------------------------------------------
      //  Dividiere dotakku[0] durch 2, mittels 1 Rechtsshift

      ptr[(a_intg)++ptr[A_END]] = ZERO;
      b_shr1 (
  &ptr[(a_intg)ptr[A_BEGIN]], (a_intg)(ptr[A_END]-ptr[A_BEGIN]+1));
      if (ptr[(a_intg)ptr[A_END]]   == ZERO) --ptr[A_END];
      if (ptr[(a_intg)ptr[A_BEGIN]] == ZERO) ++ptr[A_BEGIN];

      // --------------------------------------------------------
    }

    return rnd(dot);
  }

  real diam (interval& a)
  {
    return (subu(a.sup,a.inf));
  }

} // namespace cxsc
