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
/* CVS $Id: l_complex.inl,v 1.16 2014/01/30 17:23:46 cxsc Exp $ */

// -------------------------------------------------------------------------
// ------------------  File l_complex.inl  ---------------------------------
namespace cxsc {

// ---------------- l_complex +/- l_complex --------------------------------

inline l_complex operator +(const l_complex &a, const l_complex &b) throw()
{
   return l_complex(a.re + b.re, a.im + b.im);
}

inline l_complex operator -(const l_complex &a, const l_complex &b) throw()
{
   return l_complex(a.re - b.re, a.im - b.im);
}

inline l_complex l_complex::operator +(const complex &op2) const throw()
{
   return l_complex(this->re + Re(op2), this->im + Im(op2));
}

inline l_complex operator + (const complex& c, const l_complex& lc) throw()
{
   return lc + c;
}

inline l_complex l_complex::operator + (const real& op2) const throw()
{
   return l_complex(this->re + op2, this->im);
}

inline l_complex operator +
                   (const real& r, const l_complex& lc) throw()
{
   return lc + r;
}

inline l_complex l_complex::operator + (const l_real& op2) const throw()
{
   return l_complex(this->re + op2, this->im);
}

inline l_complex operator +
                   (const l_real& r, const l_complex& lc) throw()
{
   return lc + r;
}

inline l_complex l_complex::operator - (const l_real& op2) const throw()
{
   return l_complex(this->re - op2, this->im);
}

inline l_complex operator -
                   (const l_real& r, const l_complex& lc) throw()
{
   return -(lc - r);
}

inline l_complex l_complex::operator -(const complex &op2) const throw()
{
   return l_complex(this->re - Re(op2), this->im - Im(op2));
}

inline l_complex operator - (const complex& c, const l_complex& lc) throw()
{
   return -(lc-c);
}

inline l_complex l_complex::operator - (const real& op2) const throw()
{
   return l_complex(this->re - op2, this->im);
}

inline l_complex operator -
                   (const real& r, const l_complex& lc) throw()
{
   return -(lc-r);
}

inline l_complex l_complex::operator * (const real& op2) const throw()
{
   return l_complex(this->re * op2, this->im * op2);
}

inline l_complex operator *
                   (const real& r, const l_complex& lc) throw()
{
   return lc * r;
}

inline l_complex l_complex::operator * (const l_real& op2) const throw()
{
   return l_complex(this->re * op2, this->im * op2);
}

inline l_complex operator *
                   (const l_real& r, const l_complex& lc) throw()
{
   return lc * r;
}

// ---------------- l_complex +(-)= l_complex|complex|real|l_real ------------
inline l_complex & operator += (l_complex &lc, const l_complex &lc1) throw() 
	{
	    lc = lc + lc1; return lc;
	}
inline l_complex & operator -= (l_complex &lc, const l_complex &lc1) throw() 
	{
	    lc = lc - lc1; return lc;
	}
inline l_complex & operator += (l_complex &lc, const complex &c) throw() 
	{
	    lc = lc + c; return lc;
	}
inline l_complex & operator -= (l_complex &lc, const complex &c) throw() 
	{
	    lc = lc - c; return lc;
	}
inline l_complex & operator += (l_complex &lc, const real &r) throw() 
	{
	    lc = lc + r; return lc;
	}
inline l_complex & operator -= (l_complex &lc, const real &r) throw() 
	{
	    lc = lc - r; return lc;
	}
inline l_complex & operator += (l_complex &lc, const l_real &lr) throw() 
	{
	    lc = lc + lr; return lc;
	}
inline l_complex & operator -= (l_complex &lc, const l_real &lr) throw() 
	{
	    lc = lc - lr; return lc;
	}

// ---------------- l_complex *= l_complex|complex|real|l_real --------------
inline l_complex & operator *= (l_complex &lc, const l_complex &lc1) throw() 
	{
	    lc = lc * lc1; return lc;
	}
inline l_complex & operator *= (l_complex &lc, const complex &c) throw() 
	{
	    lc = lc * c; return lc;
	}
inline l_complex & operator *= (l_complex &lc, const real &r) throw() 
	{
	    lc = lc * r; return lc;
	}
inline l_complex & operator *= (l_complex &lc, const l_real &lr) throw() 
	{
	    lc = lc * lr; return lc;
	}

inline l_complex operator / (const l_complex& a, const complex& b) throw()
{
     return a / l_complex(b);
}

inline l_complex operator / (const l_complex& a, const l_real& b) throw()
{
     return l_complex(a.re/b, a.im/b);
}

inline l_complex operator / (const l_complex& a, const real& b) throw()
{
     return l_complex(a.re/b, a.im/b);
}

inline l_complex operator / (const complex& a, const l_complex& b) throw()
{
    return l_complex(a) / b;
}

inline l_complex operator / (const real& a, const l_complex& b) throw()
{
    return l_complex(a) / b;
}

inline l_complex operator / (const l_real& a, const l_complex& b) throw()
{
    return l_complex(a) / b;
}

inline l_complex& operator /= (l_complex& a, const l_complex& b) throw()
{ return a = a/b; }

inline l_complex& operator /= (l_complex& a, const complex& b) throw()
{ return a = a/b; }

inline l_complex& operator /= (l_complex& a, const real& b) throw()
{ return a = a/b; }

inline l_complex& operator /= (l_complex& a, const l_real& b) throw()
{ return a = a/b; }


// ---------------- Compare Operators ---------------------------------------
inline bool operator! (const l_complex & a) throw() { return !a.re && !a.im; }
inline bool operator== (const l_complex & a, const l_complex & b) throw() 
                                         { return a.re==b.re && a.im==b.im; }
inline bool operator!= (const l_complex & a, const l_complex & b) throw() 
                                         { return a.re!=b.re || a.im!=b.im; }
inline bool operator== (const l_complex & a, const complex & b) throw() 
                                       { return a.re==Re(b) && a.im==Im(b); }
inline bool operator== (const complex & a, const l_complex & b) throw() 
                                       { return Re(a)==b.re && Im(a)==b.im; }
inline bool operator!= (const l_complex & a, const complex & b) throw() 
                                       { return a.re!=Re(b) || a.im!=Im(b); }
inline bool operator!= (const complex & a, const l_complex & b) throw() 
                                       { return Re(a)!=b.re || Im(a)!=b.im; }
inline bool operator== (const l_complex & a, const real & b) throw() 
                                       { return a.re==b && !a.im; }
inline bool operator== (const real & a, const l_complex & b) throw() 
                                       { return a==b.re && !b.im; }
inline bool operator!= (const l_complex & a, const real & b) throw() 
                                       { return a.re!=b || !!a.im; }
inline bool operator!= (const real & a, const l_complex & b) throw() 
                                       { return a!=b.re || !!b.im; }
inline bool operator== (const l_complex & a, const l_real & b) throw() 
                                          { return a.re==b && !a.im; }
inline bool operator== (const l_real & a, const l_complex & b) throw() 
                                          { return a==b.re && !b.im; }
inline bool operator!= (const l_complex & a, const l_real & b) throw() 
                                         { return a.re!=b || !!a.im; }
inline bool operator!= (const l_real & a, const l_complex & b) throw() 
                                         { return a!=b.re || !!b.im; }
inline bool operator== (const l_complex & a, const dotprecision & b) throw() 
                                          { return a.re==b && !a.im; }
inline bool operator== (const dotprecision & a, const l_complex & b) throw() 
                                          { return b.re==a && !b.im; }
inline bool operator!= (const l_complex & a, const dotprecision & b) throw() 
                                          { return a.re!=b || !!a.im; }
inline bool operator!= (const dotprecision & a, const l_complex & b) throw() 
                                          { return b.re!=a || !!b.im; }

inline bool operator ==(const l_complex &c, const cdotprecision &a)
throw() { return(c.re==Re(_l_complex(a)) && c.im==Im(_l_complex(a))); }  
inline bool operator !=(const l_complex &c, const cdotprecision &a)
throw() { return(c.re!=Re(_l_complex(a)) || c.im!=Im(_l_complex(a))); } 
inline bool operator ==(const cdotprecision &a, const l_complex &c)     
throw() { return(c.re==Re(_l_complex(a)) && c.im==Im(_l_complex(a))); }  
inline bool operator !=(const cdotprecision &a, const l_complex &c)     
throw() { return(c.re!=Re(_l_complex(a)) || c.im!=Im(_l_complex(a))); } 

inline l_complex conj(const l_complex& a) throw()
{ return l_complex(a.re,-a.im); }

} // End of namespace cxsc








