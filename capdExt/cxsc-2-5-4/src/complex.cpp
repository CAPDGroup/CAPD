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

/* CVS $Id: complex.cpp,v 1.29 2014/01/30 17:23:44 cxsc Exp $ */
/* New versions of product(...), quotient(...); Blomquist, 26.10.02; */

#include "complex.hpp"
#include "dot.hpp"
#include "cdot.hpp"
#include "real.hpp"
#include "rmath.hpp"
#include "idot.hpp"
#include "interval.hpp"

#include "cinterval.hpp"

namespace cxsc {

complex::complex(const cdotprecision & a) throw()
{
   *this=rnd(a);
}

complex & complex::operator =(const cdotprecision & a) throw()
{
   return *this=rnd(a);
}

bool operator== (const complex & a, const dotprecision & b)    throw() { return !a.im && a.re==b; }
bool operator== (const dotprecision & a, const complex & b)    throw() { return !b.im && a==b.re; }
bool operator!= (const complex & a, const dotprecision & b)    throw() { return !!a.im || a.re!=b; }
bool operator!= (const dotprecision & a, const complex & b)    throw() { return !!b.im || a!=b.re; }



// ---- Ausgabefunkt. ---------------------------------------

std::ostream & operator << (std::ostream &s, const complex& a) throw()
{
   s << '('          
     << a.re << ','  
     << a.im       
     << ')';
   return s;
}
std::string & operator << (std::string &s, const complex &a) throw()
{
   s+='(';
   s << a.re;
   s+=',';
   s << a.im; 
   s+=')';
   return s;
}

std::istream & operator >> (std::istream &s, complex &a) throw()
{
   char c;

   skipeolnflag = inpdotflag = true;
   c = skipwhitespacessinglechar (s, '(');
   if (inpdotflag) 
      s.putback(c);

   s >> a.re;

   skipeolnflag = inpdotflag = true;
   c = skipwhitespacessinglechar (s, ',');
   if (inpdotflag) s.putback(c);

   s >> a.im >> RestoreOpt;

   if (!waseolnflag) 
   {
      skipeolnflag = false, inpdotflag = true;
      c = skipwhitespaces (s);
      if (inpdotflag && c != ')') 
         s.putback(c);
   }

   /*if (a.re > a.im) {
      errmon (ERR_INTERVAL(EMPTY));
   } */         
   return s;
}

std::string & operator >> (std::string &s, complex &a) throw()
{
   s = skipwhitespacessinglechar (s, '(');
   s >> SaveOpt >> RndDown >> a.re;
   s = skipwhitespacessinglechar (s, ',');
   s >> RndUp >> a.im >> RestoreOpt;
   s = skipwhitespaces (s);

   if (s[0] == ')') 
      s.erase(0,1);

    /*if (a.re > a.im) {
      errmon (ERR_INTERVAL(EMPTY));
    } */                         
   return s;
}

void operator >>(const std::string &s,complex &a) throw()
{
   std::string r(s);
   r>>a;
}
void operator >>(const char *s,complex &a) throw()
{
   std::string r(s);
   r>>a;
}

complex exp(const complex& x) throw() { return mid(exp(cinterval(x))); }
complex expm1(const complex& x) throw() { return mid(expm1(cinterval(x))); }
complex exp2(const complex& x) throw() { return mid(exp2(cinterval(x))); }
complex exp10(const complex& x) throw() { return mid(exp10(cinterval(x))); }
complex cos(const complex& x) throw() { return mid(cos(cinterval(x))); }
complex sin(const complex& x) throw() { return mid(sin(cinterval(x))); }
complex cosh(const complex& x) throw() { return mid(cosh(cinterval(x))); }
complex sinh(const complex& x) throw() { return mid(sinh(cinterval(x))); }

complex tan(const complex& x) throw() { return mid(tan(cinterval(x))); }
complex cot(const complex& x) throw() { return mid(cot(cinterval(x))); }
complex tanh(const complex& x) throw() { return mid(tanh(cinterval(x))); }
complex coth(const complex& x) throw() { return mid(coth(cinterval(x))); }

real arg(const complex& x) throw() { return mid(arg(cinterval(x))); }
real Arg(const complex& x) throw() { return mid(Arg(cinterval(x))); }

complex ln(const complex& x) throw() { return mid(ln(cinterval(x))); }
complex lnp1(const complex& x) throw() { return mid(lnp1(cinterval(x))); }
complex log2(const complex& x) throw() { return mid(log2(cinterval(x))); }
complex log10(const complex& x) throw() { return mid(log10(cinterval(x))); }

complex sqr(const complex& x) throw() { return mid(sqr(cinterval(x))); }

complex sqrt(const complex& x) throw() { return mid(sqrt(cinterval(x))); }
complex sqrtp1m1(const complex& x) throw() { return mid(sqrtp1m1(cinterval(x))); }
complex sqrt1px2(const complex& x) throw() { return mid(sqrt1px2(cinterval(x))); }
complex sqrtx2m1(const complex& x) throw() { return mid(sqrtx2m1(cinterval(x))); }
complex sqrt1mx2(const complex& x) throw() { return mid(sqrt1mx2(cinterval(x))); }
complex sqrt(const complex& x, int d) throw() 
{ return mid(sqrt(cinterval(x),d)); }

std::list<complex> sqrt_all( const complex& c )
{ 
	complex z;
	z = sqrt(c);
		
	std::list<complex> res;
	res.push_back(  z );
	res.push_back( -z );

	return res;
} 	// end sqrt_all

std::list<complex> sqrt_all( const complex& z, int n )
			//
//  sqrt_all(z,n) computes a list of n values approximating all n-th roots of z
			//
//  For n >=3, computing the optimal approximations is very expensive
//  and thus not considered cost-effective.
			//
//  Hence, the polar form is used to calculate sqrt_all(z,n)
			//
{
	std::list<complex> res;

	if( n == 0 )
	{
		res.push_back( complex(1,0) );
		return res;
	}
	else if( n == 1 )
	{
		res.push_back(z);
		return res;
	}
	else if( n == 2 ) return sqrt_all( z );
	else
	{
		real
				arg_z = arg( z ), root_abs_z = sqrt( abs( z ), n );

		for(int k = 0; k < n; k++)
		{
			real arg_k = ( arg_z + 2 * k * mid(Pi_interval) ) / n;

			res.push_back( complex( root_abs_z * cos( arg_k ),
								         root_abs_z * sin( arg_k ) ) );
		}
		return res;
	}
}
//-- end sqrt_all -------------------------------------------------------------
	

complex power_fast(const complex& x,int d) throw() 
{ return mid(power_fast(cinterval(x),d)); }
complex power(const complex& x,int d) throw() 
{ return mid(power(cinterval(x),d)); }
complex pow(const complex& x, const real& r) throw() 
{ return mid(pow(cinterval(x),interval(r))); }
complex pow(const complex& x, const complex& y) throw() 
{ return mid(pow(cinterval(x),cinterval(y))); }

complex asin(const complex& x) throw() { return mid(asin(cinterval(x))); }
complex acos(const complex& x) throw() { return mid(acos(cinterval(x))); }
complex asinh(const complex& x) throw() { return mid(asinh(cinterval(x))); }
complex acosh(const complex& x) throw() { return mid(acosh(cinterval(x))); }
complex atan(const complex& x) throw() { return mid(atan(cinterval(x))); }
complex acot(const complex& x) throw() { return mid(acot(cinterval(x))); }
complex atanh(const complex& x) throw() { return mid(atanh(cinterval(x))); }
complex acoth(const complex& x) throw() { return mid(acoth(cinterval(x))); }

} // namespace cxsc

