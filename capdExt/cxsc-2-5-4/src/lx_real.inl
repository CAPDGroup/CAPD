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

/* CVS $Id: lx_real.inl,v 1.8 2014/01/30 17:23:47 cxsc Exp $ */

/*
**  F. Blomquist, University of Wuppertal, 19.09.2007;
*/

namespace cxsc {

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
// --------- Inline functions and operators related to type lx_real ---------
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------

	inline void times2pown(lx_real &a, const real &n) throw()
{
	a = lx_real(add_real(n,a.ex),a.lr);
}

inline std::string & operator << (std::string &s, const lx_real& a) 
		throw()
// The value of a variable a of type lx_real is copied to a string s.
// s has the form:  {2**(ex)*lr}
{  
	std::stringstream ss;
	string str;

	str = "{2**(";
	s += str;
	ss << SaveOpt << SetPrecision(0,0) << Fixed << a.ex << RestoreOpt;
	ss >> str;
	s += str;
	s+=")*"; 
	s << a.lr; 
	s+='}';
	return s;
}


    inline bool Is_Integer(const real& x)
	// returns 1 if x is an integer value and 
        //           if |x| <= 2^53 - 1 = 9007199254740991.0;
        // otherwise 0 is returned
{
	double dbl1,dbl2;
	dbl1 = _double(x);
	dbl2 = floor(dbl1);
	if (dbl1 == dbl2 && fabs(dbl1) <= Max_Int_R) return 1;
	else return 0;
}

    inline real add_real(const real &a, const real &b) throw()
{
	real res(a+b);
	if (abs(res)>Max_Int_R)
		cxscthrow(REAL_INT_OUT_OF_RANGE(
					 "add_real(const real&, const real&)"));
	return res;
}

    inline real sub_real(const real &a, const real &b) throw()
{
	real res(a-b);
	if (abs(res)>Max_Int_R)
		cxscthrow(REAL_INT_OUT_OF_RANGE("sub_real(const real&, const real&)"));
	return res;
}

    inline int StagPrec(const lx_real &a) throw() 
{ return StagPrec(a.lr); }

	 inline real expo(const lx_real &a) throw()
	 { return (a.ex); }

	 inline int sign(const lx_real &a) throw() 
	 { return sign(a.lr); }

	 inline l_real lr_part(const lx_real &a) throw()
	 { return (a.lr); }

	 inline lx_real abs(const lx_real& a) throw()
	 { return lx_real(a.ex,abs(a.lr)); }
	 
	 inline lx_real adjust(const lx_real &a) throw()
	 { return lx_real(a.ex,adjust(a.lr)); }
	 
	 inline void times2pown_neg(lx_real& a, const real& n) throw()
// Calculating an approximation of  a*2^n, n = 0,-1,-2,...,-9007199254740991.0;
// n MUST be an integer and n MUST not be positive! 
// These conditions are not tested in this function!
// Blomquist, 09.11.2008;
	 {
		 int exal(expo_gr(a.lr));
		 real exa,d,n_d;
		 l_real lia(a.lr);
//		 int k;

		 if (exal>-100000) // a != [0,0]
		 {
			 exa = a.ex;
			 if (exa < -Max_Int_R - n) // exa+n < -9007199254740991;
			 {   // It holds: -Max_Int_R - n in {-Max_Int_R, ...,-1,0},
             // Thus, -Max_Int_R - n is always in integer value.
             // Furthermore it holds: exa in {-Max_Int_R,...,-2,-1}.
				 d = -Max_Int_R - exa; // d in {-Max_Int_R+1,..., -1,0} 
				 n_d = n-d;
             // With exa+n < -Max_Int_R  and with  exa+d = -Max_Int_R
             // it follows:  n-d < 0, and:
             // n-d in {-Max_Int_R,-Max_Int_R+1, ... , -1};
             // Thus, n-d is a negative and integer value.
				 if (n_d < -2100)
					 lia = 0.0;
             else  // -2100 <= n_d <0:
					 Times2pown(lia,n_d);
				 a = lx_real(-Max_Int_R,lia);
			 }
			 else // n+a.ex >= -9007199254740991, so an integer overflow
               // is not possible here!
				 a = lx_real(n+a.ex,lia);
		 }
	 } // times2pown_neg(...)


	 inline lx_real & lx_real::operator = (const lx_real &a) throw()
	 {
		 ex = a.ex;
		 lr = a.lr;
		 return *this;
	 }

	 inline lx_real & lx_real::operator = (const l_real &a) throw()
	 {
		 ex = 0;
		 lr = a;
		 return *this;
	 }

	 inline lx_real & lx_real::operator = (const real &a) throw()
	 {
		 ex = 0;
		 lr = a;
		 return *this;
	 }

// -----------------------------------------------------

	 inline bool eq_zero(const lx_real &a) throw()
	 {  return (a.lr == 0 );  }

	 inline bool gr_zero(const lx_real &a) throw()
	 {  return (a.lr > 0 );  }

	 inline bool ge_zero(const lx_real &a) throw()
	 {  return (a.lr >= 0);  }

	 inline bool sm_zero(const lx_real &a) throw()
	 {  return (a.lr < 0 );  }

	 inline bool se_zero(const lx_real &a) throw()
	 {  return (a.lr <= 0);  }

// -----------------------------------------------------

	 inline lx_real operator -(const lx_real &a) throw()
	 { return lx_real(a.ex,-a.lr); }

	 inline lx_real operator +(const lx_real &a) throw()
	 { return a; }

	 inline lx_real operator + (const lx_real& a, const l_real& b) throw()
	 { return a + lx_real(b); }
	 inline lx_real operator + (const l_real& a, const lx_real& b) throw()
	 { return lx_real(a) + b; }
	 inline lx_real operator + (const lx_real& a, const real& b) throw()
	 { return a + lx_real(l_real(b)); }
	 inline lx_real operator + (const real& a, const lx_real& b) throw()
	 { return lx_real(l_real(a)) + b; }

	 inline lx_real & operator +=(lx_real& a, const lx_real& b) throw()
	 { return a = a+b; }
	 inline lx_real & operator +=(lx_real& a, const l_real& b) throw()
	 { return a = a+b; }
	 inline lx_real & operator +=(lx_real& a, const real& b) throw()
	 { return a = a+b; }

	 inline lx_real operator - (const lx_real& a, const lx_real& b) throw()
	 { return a + lx_real(-b); }
	 inline lx_real operator - (const lx_real& a, const l_real& b) throw()
	 { return a + lx_real(-b); }
	 inline lx_real operator - (const l_real& a, const lx_real& b) throw()
	 { return lx_real(a) + lx_real(-b); }
	 inline lx_real operator - (const lx_real& a, const real& b) throw()
	 { return a + lx_real(-b); }
	 inline lx_real operator - (const real& a, const lx_real& b) throw()
	 { return lx_real(a) + lx_real(-b); }

	 inline lx_real & operator -=(lx_real& a, const lx_real& b) throw()
	 { return a = a-b; }
	 inline lx_real & operator -=(lx_real& a, const l_real& b) throw()
	 { return a = a-b; }
	 inline lx_real & operator -=(lx_real& a, const real& b) throw()
	 { return a = a-b; }

	 inline lx_real operator * (const lx_real& a, const l_real& b) throw()
	 { return a * lx_real(b); }
	 inline lx_real operator * (const l_real& a, const lx_real& b) throw()
	 { return lx_real(a) * b; }
	 inline lx_real operator * (const lx_real& a, const real& b) throw()
	 { return a * lx_real(b); }
	 inline lx_real operator * (const real& a, const lx_real& b) throw()
	 { return lx_real(a) * b; }

	 inline lx_real & operator *=(lx_real& a, const lx_real& b) throw()
	 { return a = a*b; }
	 inline lx_real & operator *=(lx_real& a, const l_real& b) throw()
	 { return a = a*b; }
	 inline lx_real & operator *=(lx_real& a, const real& b) throw()
	 { return a = a*b; }

	 inline lx_real operator / (const lx_real& a, const l_real& b) throw()
	 { return a / lx_real(b); }
	 inline lx_real operator / (const l_real& a, const lx_real& b) throw()
	 { return lx_real(a) / b; }
	 inline lx_real operator / (const lx_real& a, const real& b) throw()
	 { return a / lx_real(b); }
	 inline lx_real operator / (const real& a, const lx_real& b) throw()
	 { return lx_real(a) / b; }

	 inline lx_real & operator /=(lx_real& a, const lx_real& b) throw()
	 { return a = a/b; }
	 inline lx_real & operator /=(lx_real& a, const l_real& b) throw()
	 { return a = a/b; }
	 inline lx_real & operator /=(lx_real& a, const real& b) throw()
	 { return a = a/b; }

	 inline bool operator ! (const lx_real& a) throw()
	 { return !a.lr; }

	 inline bool operator == (const lx_real &a, const l_real &b) throw()
	 {  return (a==lx_real(b));  }

	 inline bool operator == (const l_real &a, const lx_real &b) throw()
	 {  return (lx_real(a)==b);  }

	 inline bool operator == (const lx_real &a, const real &b) throw()
	 {  return (a==lx_real(b));  }

	 inline bool operator == (const real &a, const lx_real &b) throw()
	 {  return (lx_real(a)==b);  }


	 inline bool operator != (const lx_real &a, const lx_real &b) throw()
	 {  return !(a==b); }

	 inline bool operator != (const lx_real &a, const l_real &b) throw()
	 {  return !(a==lx_real(b));  }

	 inline bool operator != (const l_real &a, const lx_real &b) throw()
	 {  return !(lx_real(a)==b);  }

	 inline bool operator != (const lx_real &a, const real &b) throw()
	 {  return !(a==lx_real(b));  }

	 inline bool operator != (const real &a, const lx_real &b) throw()
	 {  return !(lx_real(a)==b);  }

	 inline bool operator <= (const lx_real &a, const lx_real &b) throw()
	 { return !(a>b); }

	 inline bool operator < (const lx_real &a, const lx_real &b) throw()
	 {  return (b>a); }

	 inline bool operator >= (const lx_real &a, const lx_real &b) throw()
	 { return !(a<b); }

// ---------------------------------------------------------

	 inline bool operator > (const real &a, const lx_real &b) throw()
	 { return lx_real(a)>b; }

	 inline bool operator <= (const real &a, const lx_real &b) throw()
	 { return !(a>b); }

	 inline bool operator <  (const real &a, const lx_real &b) throw()
	 {  return b>lx_real(a); }

	 inline bool operator >= (const real &a, const lx_real &b) throw()
	 { return !(a<b); }

// ---------------------------------------------------------

	 inline bool operator >  (const lx_real &a, const real &b) throw()
	 { return a>lx_real(b); }

	 inline bool operator <= (const lx_real &a, const real &b) throw()
	 { return !(a>b); }

	 inline bool operator <  (const lx_real &a, const real &b) throw()
	 {  return b>a; }

	 inline bool operator >= (const lx_real &a, const real &b) throw()
	 { return !(a<b); }

// ---------------------------------------------------------

	 inline bool operator >  (const l_real &a, const lx_real &b) throw()
	 { return lx_real(a)>b; }

	 inline bool operator <= (const l_real &a, const lx_real &b) throw()
	 { return !(a>b); }

	 inline bool operator <  (const l_real &a, const lx_real &b) throw()
	 {  return b>lx_real(a); }

	 inline bool operator >= (const l_real &a, const lx_real &b) throw()
	 { return !(a<b); }

// ---------------------------------------------------------

	 inline bool operator >  (const lx_real &a, const l_real &b) throw()
	 { return a>lx_real(b); }

	 inline bool operator <= (const lx_real &a, const l_real &b) throw()
	 { return !(a>b); }

	 inline bool operator <  (const lx_real &a, const l_real &b) throw()
	 {  return b>a; }

	 inline bool operator >= (const lx_real &a, const l_real &b) throw()
	 { return !(a<b); }

// -----------------------------------------------------


	 inline lx_real max(const lx_real& a, const lx_real& b)
	 { return (b>a)? b : a; }

	 inline lx_real min(const lx_real& a, const lx_real& b)
	 { return (b>a)? a : b; }

	 inline real cutint(const real& x) throw()
// y = cutint(x) delivers the truncated part y of x.
// If y is not an integer value then 9007199254740992.0
// is returned;      
// For all integer values y it holds:
//               |y| <= Max_Int_R := 9007199254740991.0
// Examples:
// y = cutint(-0.1);     --->  y = 0;
// y = cutint(-123.5);   --->  y = -123.0;
// y = cutint(+123.5);   --->  y = +123.0;
// y = cutint(9007199254740991.8); --->  y = 9007199254740992.0;
// y = cutint(9007199254740992.0); --->  y = 9007199254740992.0;
// y = cutint(-1e20);              --->  y = 9007199254740992.0;
// Blomquist, 26.05.2008;
	 {
		 real res(x);
		 double dbl;
		 bool neg;
		 neg = x<0;
		 if (neg) res = -res; // res = |x|
		 if (res>Max_Int_R) res = 9007199254740992.0;
		 else
		 {
			 dbl = _double(res);
			 dbl = floor(dbl);
			 res = real(dbl);
			 if (neg) res = -res;
		 }
		 return res;
	 }

} // end namespace cxsc
