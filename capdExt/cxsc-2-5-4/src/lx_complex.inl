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

/* CVS $Id: lx_complex.inl,v 1.8 2014/01/30 17:23:47 cxsc Exp $ */

/*
**  F. Blomquist, University of Wuppertal, 19.09.2007;
*/

namespace cxsc {
// --------------------------------------------------------------------------
// ------- Inline functions and operators related to type lx_complex --------
// --------------------------------------------------------------------------

inline lx_real Re(const lx_complex &a)
{ return a.re; }

inline lx_real Im(const lx_complex &a)
{ return a.im; }

inline int StagPrec(const lx_complex &a) throw() 
{ return StagPrec(a.re); }

inline real expoRe(const lx_complex &a) throw()
{ return expo(a.re); }

inline real expoIm(const lx_complex &a) throw()
{ return expo(a.im); }

inline l_real lr_partRe(const lx_complex &a) throw()
{ return lr_part(a.re); }

inline l_real lr_partIm(const lx_complex &a) throw()
{ return lr_part(a.im); }

inline lx_complex & SetRe(lx_complex &a, const lx_real &b)
{ a.re = b; return a; } // The real part of a is substituted by b.
inline lx_complex & SetRe(lx_complex &a, const l_real &b)
{ a.re = b; return a; } // The real part of a is substituted by b.
inline lx_complex & SetRe(lx_complex &a, const real &b)
{ a.re = b; return a; } // The real part of a is substituted by b.

inline lx_complex & SetIm(lx_complex &a, const lx_real &b)
{ a.im = b; return a; } // The imaginary part of a is substituted by b.
inline lx_complex & SetIm(lx_complex &a, const l_real &b)
{ a.im = b; return a; } // The imaginary part of a is substituted by b.
inline lx_complex & SetIm(lx_complex &a, const real &b)
{ a.im = b; return a; } // The imaginary part of a is substituted by b.

inline lx_complex conj(const lx_complex& a) throw()
{ return lx_complex(a.re, -a.im); }

inline bool operator ! (const lx_complex& a) throw()
{ return !a.re && !a.im; }

inline bool operator == (const lx_complex &a, const lx_complex &b) 
		throw() { return (a.re == b.re && a.im == b.im); }

inline bool operator == (const lx_complex &a, const l_complex &b) 
		throw() { return (a.re == Re(b) && a.im == Im(b)); }
inline bool operator == (const lx_complex &a, const complex &b) 
		throw() { return (a.re == Re(b) && a.im == Im(b)); }
inline bool operator == (const l_complex &a, const lx_complex &b) 
		throw() { return b == a; }
inline bool operator == (const complex &a, const lx_complex &b) 
		throw() { return (b == a); }

inline bool operator == (const lx_complex &a, const lx_real &b) throw()
{ return a.re == b && a.im == 0.0; }
inline bool operator == (const lx_complex &a, const l_real &b) throw()
{ return a.re == b && a.im == 0.0; }
inline bool operator == (const lx_complex &a, const real &b) throw()
{ return a.re == b && a.im == 0.0; }

inline bool operator == (const lx_real &a, const lx_complex &b) throw()
{ return a == b.re && b.im == 0.0; }
inline bool operator == (const l_real &a, const lx_complex &b) throw()
{ return a == b.re && b.im == 0.0; }
inline bool operator == (const real &a,   const lx_complex &b) throw()
{ return a == b.re && b.im == 0.0; }

inline bool operator != (const lx_complex &a, const lx_complex &b) throw()
{ return !(a == b); }

inline bool operator != (const lx_complex &a, const l_complex &b) throw()
{ return !(a == b); }
inline bool operator != (const lx_complex &a, const complex   &b) throw()
{ return !(a == b); }
inline bool operator != (const l_complex &a, const lx_complex &b) throw()
{ return !(a == b); }
inline bool operator != (const complex   &a, const lx_complex &b) throw()
{ return !(a == b); }

inline bool operator != (const lx_complex &a, const lx_real &b) throw()
{ return !(a == b); }
inline bool operator != (const lx_complex &a, const l_real &b) throw()
{ return !(a == b); }
inline bool operator != (const lx_complex &a, const real &b)   throw()
{ return !(a == b); }
inline bool operator != (const lx_real &a, const lx_complex &b) throw()
{ return !(a == b); }
inline bool operator != (const l_real &a, const lx_complex &b) throw()
{ return !(a == b); }
inline bool operator != (const real &a,   const lx_complex &b) throw()
{ return !(a == b); }

inline lx_complex operator + (const lx_complex &a) throw()
{ return a; }
inline lx_complex operator - (const lx_complex &a) throw()
{ return lx_complex(-a.re,-a.im); }

inline lx_complex operator + (const lx_complex& a, const lx_complex& b) throw()
{ return lx_complex(a.re+b.re,a.im+b.im); }

inline lx_complex operator + (const lx_complex& a, const l_complex& b) throw()
{ return a + lx_complex(b); }
inline lx_complex operator + (const lx_complex& a, const complex& b) throw()
{ return a + lx_complex(b); }
inline lx_complex operator + (const l_complex& a, const lx_complex& b) throw()
{ return lx_complex(a) + b; }
inline lx_complex operator + (const complex& a, const lx_complex& b) throw()
{ return lx_complex(a) + b; }
inline lx_complex operator + (const lx_complex& a, const lx_real& b) throw()
{ return lx_complex(a.re + b, Im(a)); }
inline lx_complex operator + (const lx_real& a, const lx_complex& b) throw()
{ return lx_complex(b.re + a, Im(b)); }
inline lx_complex operator + (const lx_complex& a, const l_real& b) throw()
{ return lx_complex(a.re + b, Im(a)); }
inline lx_complex operator + (const l_real& a, const lx_complex& b) throw()
{ return lx_complex(b.re + a, Im(b)); }
inline lx_complex operator + (const lx_complex& a, const real& b) throw()
{ return lx_complex(a.re + b, Im(a)); }
inline lx_complex operator + (const real& a, const lx_complex& b) throw()
{ return lx_complex(b.re + a, Im(b)); }

inline lx_complex & operator +=(lx_complex& a, const lx_complex& b) throw()
{ return a = a+b; }
inline lx_complex & operator +=(lx_complex& a, const l_complex& b) throw()
{ return a = a+b; }
inline lx_complex & operator +=(lx_complex& a, const complex& b) throw()
{ return a = a+b; }
inline lx_complex & operator +=(lx_complex& a, const lx_real& b) throw()
{ return a = a+b; }
inline lx_complex & operator +=(lx_complex& a, const l_real& b) throw()
{ return a = a+b; }
inline lx_complex & operator +=(lx_complex& a, const real& b) throw()
{ return a = a+b; }

inline lx_complex operator - (const lx_complex& a, const lx_complex& b) throw()
{ return a + (-b); }
inline lx_complex operator - (const lx_complex& a, const l_complex& b) throw()
{ return a + (-b); }
inline lx_complex operator - (const lx_complex& a, const complex& b) throw()
{ return a + (-b); }
inline lx_complex operator - (const l_complex& a, const lx_complex& b) throw()
{ return a + (-b); }
inline lx_complex operator - (const complex& a, const lx_complex& b) throw()
{ return a + (-b); }
inline lx_complex operator - (const lx_complex& a, const lx_real& b) throw()
{ return a + (-b); }
inline lx_complex operator - (const lx_complex& a, const l_real& b) throw()
{ return a + (-b); }
inline lx_complex operator - (const lx_complex& a, const real& b) throw()
{ return a + (-b); }
inline lx_complex operator - (const lx_real& a, const lx_complex& b) throw()
{ return a + (-b); }
inline lx_complex operator - (const l_real& a, const lx_complex& b) throw()
{ return a + (-b); }
inline lx_complex operator - (const real& a, const lx_complex& b) throw()
{ return a + (-b); }

inline lx_complex & operator -=(lx_complex& a, const lx_complex& b) throw()
{ return a = a-b; }
inline lx_complex & operator -=(lx_complex& a, const l_complex& b) throw()
{ return a = a-b; }
inline lx_complex & operator -=(lx_complex& a, const complex& b) throw()
{ return a = a-b; }
inline lx_complex & operator -=(lx_complex& a, const lx_real& b) throw()
{ return a = a-b; }
inline lx_complex & operator -=(lx_complex& a, const l_real& b) throw()
{ return a = a-b; }
inline lx_complex & operator -=(lx_complex& a, const real& b) throw()
{ return a = a-b; }


inline lx_complex operator * (const lx_complex& a, const lx_complex& b) throw()
{
	lx_real x,y;
	
	x = a.re*b.re - a.im*b.im;
	y = a.im*b.re + a.re*b.im;
	
	return lx_complex(x,y);
}
inline lx_complex operator * (const lx_complex& a, const l_complex& b) throw()
{ return a*lx_complex(b); }
inline lx_complex operator * (const lx_complex& a, const complex& b) throw()
{ return a*lx_complex(b); }
inline lx_complex operator * (const l_complex& a, const lx_complex& b) throw()
{ return lx_complex(a)*b; }
inline lx_complex operator * (const complex& a, const lx_complex& b) throw()
{ return lx_complex(a)*b; }
inline lx_complex operator * (const lx_complex& a, const lx_real& b) throw()
{ return a*lx_complex(b); }
inline lx_complex operator * (const lx_complex& a, const l_real& b) throw()
{ return a*lx_complex(b); }
inline lx_complex operator * (const lx_complex& a, const real& b) throw()
{ return a*lx_complex(b); }
inline lx_complex operator * (const lx_real& a, const lx_complex& b) throw()
{ return lx_complex(a)*b; }
inline lx_complex operator * (const l_real& a, const lx_complex& b) throw()
{ return lx_complex(a)*b; }
inline lx_complex operator * (const real& a, const lx_complex& b) throw()
{ return lx_complex(a)*b; }

inline lx_complex & operator *=(lx_complex& a, const lx_complex& b) throw()
{ return a = a*b; }
inline lx_complex & operator *=(lx_complex& a, const l_complex& b) throw()
{ return a = a*b; }
inline lx_complex & operator *=(lx_complex& a, const complex& b) throw()
{ return a = a*b; }
inline lx_complex & operator *=(lx_complex& a, const lx_real& b) throw()
{ return a = a*b; }
inline lx_complex & operator *=(lx_complex& a, const l_real& b) throw()
{ return a = a*b; }
inline lx_complex & operator *=(lx_complex& a, const real& b) throw()
{ return a = a*b; }


inline lx_complex operator / (const lx_complex& a, const lx_complex& b) throw()
{
	lx_real x,y,Ne;
		
	Ne = b.re*b.re + b.im*b.im;
	x = (a.re*b.re + a.im*b.im) / Ne;
	y = (a.im*b.re - a.re*b.im) / Ne;
	return lx_complex(x,y);
}
inline lx_complex operator / (const lx_complex& a, const l_complex& b) throw()
{ return a/lx_complex(b); }
inline lx_complex operator / (const lx_complex& a, const complex& b) throw()
{ return a/lx_complex(b); }
inline lx_complex operator / (const l_complex& a, const lx_complex& b) throw()
{ return lx_complex(a)/b; }
inline lx_complex operator / (const complex& a, const lx_complex& b) throw()
{ return lx_complex(a)/b; }
inline lx_complex operator / (const lx_complex& a, const lx_real& b) throw()
{ return a/lx_complex(b); }
inline lx_complex operator / (const lx_complex& a, const l_real& b) throw()
{ return a/lx_complex(b); }
inline lx_complex operator / (const lx_complex& a, const real& b) throw()
{ return a/lx_complex(b); }
inline lx_complex operator / (const lx_real& a, const lx_complex& b) throw()
{ return lx_complex(a)/b; }
inline lx_complex operator / (const l_real& a, const lx_complex& b) throw()
{ return lx_complex(a)/b; }
inline lx_complex operator / (const real& a, const lx_complex& b) throw()
{ return lx_complex(a)/b; }

inline lx_complex & operator /=(lx_complex& a, const lx_complex& b) throw()
{ return a = a/b; }
inline lx_complex & operator /=(lx_complex& a, const l_complex& b) throw()
{ return a = a/b; }
inline lx_complex & operator /=(lx_complex& a, const complex& b) throw()
{ return a = a/b; }
inline lx_complex & operator /=(lx_complex& a, const lx_real& b) throw()
{ return a = a/b; }
inline lx_complex & operator /=(lx_complex& a, const l_real& b) throw()
{ return a = a/b; }
inline lx_complex & operator /=(lx_complex& a, const real& b) throw()
{ return a = a/b; }

// --------------------------- Output ---------------------------------

inline std::ostream& operator << (std::ostream& s, const lx_complex& a) 
		throw()
// A value a of type lx_complex is written to the output channel.
// The output has the form:  { ? , ? }
{     
	s << '('
			<< a.re
			<< " , "
			<< a.im
			<< ')';
	return s;
}

inline std::string & operator << (std::string &s, const lx_complex& a) throw()
// The value of a variable a of type lx_complex is copied to a string s.
// s has the form:  ({2**(...)*...} , {2**(...)*...})
{
	string str;
	s += "(";
	str = "";
	str << a.re;
	s += str;
	s += " , ";
	str = "";
	str << a.im;
	s += str;
	s += ")";
	return s;
}

} // end namespace cxsc
