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

/* CVS $Id: real.inl,v 1.33 2014/01/30 17:23:48 cxsc Exp $ */

#ifdef _CXSC_REAL_HPP_INCLUDED

#include "fi_lib.hpp"
#include <cmath>

namespace cxsc {
extern "C"
{
   a_real r_comp(a_real,a_intg);
   a_real r_mant(a_real);
   a_intg r_expo(a_real);
}

inline double _double(const real &a) throw()
{ 
   return a.w; 
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::real::real(const double &a)
*/
inline real   _real(const double &a) throw()
{ 
   return real(a); 
}

inline real pred(const real& r) throw()
{
   real ret;
   ret = fi_lib::q_pred(r.w);
   return ret;
}

inline real succ(const real& r) throw()
{
   real ret;
   ret = fi_lib::q_succ(r.w);
   return ret;
}

inline a_intg expo(const real &r) throw()
{
   return r_expo(*(a_real*)&r.w);
}

inline real comp(const real &r, a_intg i) throw()
{
   return real(r_comp(*(a_real*)&r.w,i));
}

inline real mant(const real &r) throw()
{
   return real(r_mant(*(a_real*)&r.w));
}

inline real operator -(const real &a) throw()       { return -a.w; }
inline real operator +(const real &a) throw()       { return a; }

inline real operator +(const real &a,const real &b) throw () { return a.w+b.w; }
inline real operator -(const real &a,const real &b) throw () { return a.w-b.w; }
inline real operator *(const real &a,const real &b) throw () { return a.w*b.w; }
inline real operator /(const real &a,const real &b) throw () { return a.w/b.w; }

inline real& operator +=(real &a, const real &b) throw () { a.w+=b.w; return a; }
inline real& operator -=(real &a, const real &b) throw () { a.w-=b.w; return a; }
inline real& operator *=(real &a, const real &b) throw () { a.w*=b.w; return a; }
inline real& operator /=(real &a, const real &b) throw () { a.w/=b.w; return a; }

inline bool operator!  (const real& a)                throw () { return !a.w; }
inline bool operator== (const real& a, const real& b) throw () { return a.w==b.w; }
inline bool operator!= (const real& a, const real& b) throw () { return a.w!=b.w; }
inline bool operator<  (const real& a, const real& b) throw () { return a.w<b.w;  }
inline bool operator<= (const real& a, const real& b) throw () { return a.w<=b.w; }
inline bool operator>= (const real& a, const real& b) throw () { return a.w>=b.w; }
inline bool operator>  (const real& a, const real& b) throw () { return a.w>b.w;  }

inline bool operator== (const real& a, const int & b) throw () { return a.w==b; }
inline bool operator!= (const real& a, const int & b) throw () { return a.w!=b; }
inline bool operator== (const int & a, const real& b) throw () { return a==b.w; }
inline bool operator!= (const int & a, const real& b) throw () { return a!=b.w; }
inline bool operator== (const real& a, const long & b) throw () { return a.w==b; }
inline bool operator!= (const real& a, const long & b) throw () { return a.w!=b; }
inline bool operator== (const long & a, const real& b) throw () { return a==b.w; }
inline bool operator!= (const long & a, const real& b) throw () { return a!=b.w; }
inline bool operator== (const real& a, const float & b) throw () { return a.w==b; }
inline bool operator!= (const real& a, const float & b) throw () { return a.w!=b; }
inline bool operator== (const float & a, const real& b) throw () { return a==b.w; }
inline bool operator!= (const float & a, const real& b) throw () { return a!=b.w; }
inline bool operator== (const real& a, const double & b) throw () { return a.w==b; }
inline bool operator!= (const real& a, const double & b) throw () { return a.w!=b; }
inline bool operator== (const double & a, const real& b) throw () { return a==b.w; }
inline bool operator!= (const double & a, const real& b) throw () { return a!=b.w; }


inline real abs(const real &a) throw()
{
  /* if(a.w>=0)
      return a.w;
   return -a.w;*/
  return fabs(a.w);
}

inline int sign(const real &a) throw()
{
   if(a.w>0)
      return 1;
   else
      if(a.w)
         return -1;
      else 
         return 0;
}   

inline void times2pown(real& r, const int n) // Blomquist 01.10.02.
{   // times2pown(r,n): The reference parameter r gets the new value:  r*2^n; 
    // r*2^n in the normalized range --> r*2^n is the exact value.
    // r*2^n in the denormalized range --> r*2^n is in general not exact.
    int j = expo(r) + n;
    if (j >= -1021) { r = comp(mant(r),j); } // r is normalized up to now !
    else 
    {   // result now in the denormalized range:
	j += 1021;
	r = comp(mant(r), -1021); // r is still normalized
	if (j<-53) r = 0;
        else r = r * comp(0.5,j+1);
    } 
} // end of times2pown(...)

inline real pow2n(const int n) throw() // Blomquist 03.10.02.
{   // Fast and exact calculation of 2^n; -1074 <= n <= +1023;
    return comp(0.5,n+1);
}

inline real max(const real & a, const real & b) { return a>b?a:b; }
inline real min(const real & a, const real & b) { return a<b?a:b; }
inline real Max(const real & a, const real & b) { return a>b?a:b; } 


// ------- Verknuepfungen mit nach oben bzw. unten gerichteter Rundung ------
// ------- (Operators with rounding direction upwards or downwards) ---------

} // namespace cxsc

#include "rts_real.hpp"

#if ROUND_ASM
#include <r_ari.h>
#endif

#if ROUND_C99_SAVE+ROUND_C99_QUICK
#include <fenv.h>
#endif

#ifdef _MSC_VER
#include <float.h>
#pragma fenv_access (on)
	static inline void setround(int r) {
		unsigned int control_word;
		_controlfp_s(&control_word,0,0);
		switch(r) {
			case -1:
				_controlfp_s(&control_word,_RC_DOWN,_MCW_RC);
				break;
			case 0:
				_controlfp_s(&control_word,_RC_NEAR,_MCW_RC);
				break;
			case 1:
				_controlfp_s(&control_word,_RC_UP,_MCW_RC);
				break;
			case 2: 
				_controlfp_s(&control_word,_RC_CHOP,_MCW_RC);
				break;
			default:
				_controlfp_s(&control_word,_RC_NEAR,_MCW_RC);
		}
	}

#endif

#if ROUND_C96_SAVE+ROUND_C96_QUICK
#include <fenv96.h>
#define fegetround fegetround96
#define fesetround fesetround96
#endif


namespace cxsc {

inline real addup(const real &x, const real &y)
{
	#if ROUND_ASM
		double ret;
		ret = ra_addu(x.w,y.w);
		return real(ret);	        	
	#elif ROUND_C99_SAVE+ROUND_C96_SAVE
                int mode;
		mode=fegetround();
                fesetround(FE_UPWARD);
		double ret;
		ret = x.w + y.w;
                fesetround(mode);
		return real(ret);
	#elif ROUND_C99_QUICK+ROUND_C96_QUICK
	        fesetround(FE_UPWARD);
                return( x.w + y.w);
    #elif _MSC_VER
        unsigned int cw;
		unsigned int mask = 0xFFFFF;
		_controlfp_s(&cw,0,0);
        setround(1);
		double ret;
		ret = x.w + y.w;
		_controlfp_s(&cw,cw,mask);
		return real(ret);
	#else
	        a_real ret;
		ret = r_addu(_a_real(x.w), _a_real(y.w));
        	return _real(ret);
	#endif
}

inline real adddown(const real &x, const real &y)
{
	#if ROUND_ASM
		double ret;
		ret = ra_addd(x.w,y.w);
		return real(ret);
	#elif ROUND_C99_SAVE+ROUND_C96_SAVE
                int mode;
		mode=fegetround();
                fesetround(FE_DOWNWARD);
		double ret;
		ret = x.w + y.w;
                fesetround(mode);
		return real(ret);
	#elif ROUND_C99_QUICK+ROUND_C96_QUICK
	        fesetround(FE_DOWNWARD);
                return( x.w + y.w);
    #elif _MSC_VER
        unsigned int cw;
		unsigned int mask = 0xFFFFF;
		_controlfp_s(&cw,0,0);
        setround(-1);
		double ret;
		ret = x.w + y.w;
		_controlfp_s(&cw,cw,mask);
		return real(ret);
	#else
	a_real ret;
		ret = r_addd(_a_real(x.w), _a_real(y.w));
        	return _real(ret);
	#endif
}

inline real subup(const real &x, const real &y)
{
	#if ROUND_ASM
	        double ret;
	        ret = ra_subu(x.w,y.w);
	    	return real(ret);
	#elif ROUND_C99_SAVE+ROUND_C96_SAVE
                int mode;
		mode=fegetround();
                fesetround(FE_UPWARD);
		double ret;
		ret = x.w - y.w;
                fesetround(mode);
		return real(ret);
	#elif ROUND_C99_QUICK+ROUND_C96_QUICK
	        fesetround(FE_UPWARD);
                return( x.w - y.w);
    #elif _MSC_VER
        unsigned int cw;
		unsigned int mask = 0xFFFFF;
		_controlfp_s(&cw,0,0);
        setround(1);
		double ret;
		ret = x.w - y.w;
		_controlfp_s(&cw,cw,mask);
		return real(ret);
	#else
	        a_real ret;
		ret = r_subu(_a_real(x.w), _a_real(y.w));
        	return _real(ret);
	#endif
}

inline real subdown(const real &x, const real &y)
{
	#if ROUND_ASM
		double ret;
		ret = ra_subd(x.w,y.w);
        	return real(ret);		        	
	#elif ROUND_C99_SAVE+ROUND_C96_SAVE
                int mode;
		mode=fegetround();
                fesetround(FE_DOWNWARD);
		double ret;
		ret = x.w - y.w;
                fesetround(mode);
		return real(ret);
	#elif ROUND_C99_QUICK+ROUND_C96_QUICK
	        fesetround(FE_DOWNWARD);
                return( x.w - y.w);
    #elif _MSC_VER
        unsigned int cw;
		unsigned int mask = 0xFFFFF;
		_controlfp_s(&cw,0,0);
        setround(-1);
		double ret;
		ret = x.w - y.w;
		_controlfp_s(&cw,cw,mask);
		return real(ret);
	#else
	        a_real ret;
		ret = r_subd(_a_real(x.w), _a_real(y.w));
        	return _real(ret);
	#endif
}

inline real multup(const real &x, const real &y)
{
	#if ROUND_ASM
	        double ret;
	        ret = ra_mulu(x.w,y.w);
        	return real(ret);
	#elif ROUND_C99_SAVE+ROUND_C96_SAVE
                int mode;
		mode=fegetround();
                fesetround(FE_UPWARD);
		double ret;
		ret = x.w * y.w;
                fesetround(mode);
		return real(ret);
	#elif ROUND_C99_QUICK+ROUND_C96_QUICK
	        fesetround(FE_UPWARD);
                return( x.w * y.w);
    #elif _MSC_VER
        unsigned int cw;
		unsigned int mask = 0xFFFFF;
		_controlfp_s(&cw,0,0);
        setround(1);
		double ret;
		ret = x.w * y.w;
		_controlfp_s(&cw,cw,mask);
		return real(ret);
	#else
	        a_real ret;
		ret = r_mulu(_a_real(x.w), _a_real(y.w));
        	return _real(ret);
	#endif
}

inline real multdown(const real &x, const real &y)
{
	#if ROUND_ASM
	        double ret;
	        ret = ra_muld(x.w,y.w);
	    	return real(ret);
	#elif ROUND_C99_SAVE+ROUND_C96_SAVE
                int mode;
		mode=fegetround();
                fesetround(FE_DOWNWARD);
		double ret;
		ret = x.w * y.w;
                fesetround(mode);
		return real(ret);
	#elif ROUND_C99_QUICK+ROUND_C96_QUICK
	        fesetround(FE_DOWNWARD);
                return( x.w * y.w);
    #elif _MSC_VER
        unsigned int cw;
		unsigned int mask = 0xFFFFF;
		_controlfp_s(&cw,0,0);
        setround(-1);
		double ret;
		ret = x.w * y.w;
		_controlfp_s(&cw,cw,mask);
		return real(ret);
	#else
	        a_real ret;
		ret = r_muld(_a_real(x.w), _a_real(y.w));
        	return _real(ret);
	#endif
}

inline real divup(const real &x, const real &y)
{
	#if ROUND_ASM
		double ret;
	        ret = ra_divu(x.w,y.w);
        	return real(ret);
	#elif ROUND_C99_SAVE+ROUND_C96_SAVE
                int mode;
		mode=fegetround();
                fesetround(FE_UPWARD);
		double ret;
		ret = x.w / y.w;
                fesetround(mode);
		return real(ret);
	#elif ROUND_C99_QUICK+ROUND_C96_QUICK
	        fesetround(FE_UPWARD);
                return( x.w / y.w);
    #elif _MSC_VER
        unsigned int cw;
		unsigned int mask = 0xFFFFF;
		_controlfp_s(&cw,0,0);
        setround(1);
		double ret;
		ret = x.w / y.w;
		_controlfp_s(&cw,cw,mask);
		return real(ret);
	#else
	        a_real ret;
		ret = r_divu(_a_real(x.w), _a_real(y.w));
        	return _real(ret);
	#endif
}

inline real divdown(const real &x, const real &y)
{
	#if ROUND_ASM
		double ret;
	        ret = ra_divd(x.w,y.w);
        	return real(ret);
	#elif ROUND_C99_SAVE+ROUND_C96_SAVE
                int mode;
		mode=fegetround();
                fesetround(FE_DOWNWARD);
		double ret;
		ret = x.w / y.w;
                fesetround(mode);
		return real(ret);
	#elif ROUND_C99_QUICK+ROUND_C96_QUICK
	        fesetround(FE_DOWNWARD);
                return( x.w / y.w);
    #elif _MSC_VER
        unsigned int cw;
		unsigned int mask = 0xFFFFF;
		_controlfp_s(&cw,0,0);
        setround(-1);
		double ret;
		ret = x.w / y.w;
		_controlfp_s(&cw,cw,mask);
		return real(ret);
	#else
	        a_real ret;
		ret = r_divd(_a_real(x.w), _a_real(y.w));
        	return _real(ret);
	#endif
}

//----------------------------------------------------------------------------
// IsInfinity - prueft ob die uebergebene IEEE 64-Bit Gleitkommazahl
//              'Unendlich' darstellt, dies ist der Fall, falls:
//                alle Bits des Exponenten 1 sind
//                alle Bits der Mantisse 0 sind
//
inline bool IsInfinity (const real &a) 
{
   if ((((a_btyp*)&a)[HIGHREAL] & 0x7FF00000L) != 0x7FF00000L) 
      return false;
   if ((((a_btyp*)&a)[HIGHREAL] & 0x000FFFFFL) != 0L) 
      return false;
   if (((a_btyp*)&a)[LOWREAL] != 0L) 
      return false;
  return true;  // a ist 'unendlich'
}

//----------------------------------------------------------------------------
// IsNan      - prueft ob die uebergebene IEEE 64-Bit Gleitkommazahl
//              eine ungueltige Zahl (Not a number) darstellt,
//              dies ist der Fall, falls:
//                alle Bits des Exponenten 1 sind
//                und nicht alle Bits der Mantisse 0 sind
//
inline bool IsQuietNaN (const real &a) 
{
   if ((((a_btyp*)&a)[HIGHREAL] & 0x7FF00000L) != 0x7FF00000L) 
      return false;

   if ((((a_btyp*)&a)[HIGHREAL] & 0x000FFFFFL) == 0L) 
   {
      if (((a_btyp*)&a)[LOWREAL] == 0L) 
         return false;
   }
   return true;
}

//----------------------------------------------------------------------------
// IsSignalingNaN - prueft ob die uebergebene IEEE 64-Bit Gleitkommazahl
//              undefiniert ist (z.B. Berechnung der Wurzel einer negativen
//              Zahl), dies ist der Fall, falls
//                das Vorzeichenbit 1 ist
//                alle Bits des Exponenten 1 sind
//                das hoechste Bit der Mantisse 1 ist
//                und alle anderen Bits der Mantisse 0 sind
//
inline bool IsSignalingNaN (const real &a) 
{
   if ((((a_btyp*)&a)[HIGHREAL] & 0xFFF00000L) != 0xFFF00000L) 
      return false;
   if ((((a_btyp*)&a)[HIGHREAL] & 0x000FFFFFL) != 0x00080000L) 
      return false;
   if (((a_btyp*)&a)[LOWREAL] != 0L) 
      return false;
   return true;
}

} // namespace cxsc

#endif // _CXSC_REAL_HPP_INCLUDED
