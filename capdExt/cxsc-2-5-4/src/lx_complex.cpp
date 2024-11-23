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

/* CVS $Id: lx_complex.cpp,v 1.9 2014/01/30 17:23:47 cxsc Exp $ */

/*
**  F. Blomquist, University of Wuppertal, 19.09.2007;
*/

#include "lx_complex.hpp"
#include "lx_cinterval.hpp"

namespace cxsc {
// ----------------------------------------------------------------------------
// ----- Functions and operators related to type lx_complex -------------------
// ----------------------------------------------------------------------------

	l_complex & l_complex::operator = (const lx_complex& a) throw()
{
	l_real u,v;
	u = Re(a);  v = Im(a);
	return *this = l_complex(u,v);
}
	
	complex & complex::operator = (const lx_complex& a) throw()
{
	l_complex x;
	complex y;
		
	x = a;
	y = x;
	return *this = y;
}

	std::string & operator >> (std::string& s, lx_complex& a) throw()
// Writes string s to variable a of type lx_complex;
// and returns an empty string s;
// Example:  s = "({2,3} , {4,5})" delivers a value a
// with:    a ~ 10^2 * 3.0 + i*10^4 * 5.0;
{
	string su;
	int i;
	std::cout << "Halo 1" << std::endl;
	s = skipwhitespacessinglechar (s, '(');
	std::cout << "s = " << s << std::endl;
	i = s.find("}");
	std::cout << "i = " << i << std::endl;
	su = s.substr(0,i+1);
	std::cout << "su = " << su << std::endl;
	su >> a.re;
		
	s.erase(0,i+1);
	s = skipwhitespacessinglechar (s, ',');
	std::cout << "s = " << s << std::endl;
	s >> a.im;

	s = "";
	return s;
}
	
	void operator >> (const std::string &s, lx_complex &a) throw()
{
   // Writes string s to variable a of type lx_real;
	std::string r(s);
	r >> a;
}
	
	void operator >> (const char *s, lx_complex& a) throw()
{
	std::string r(s);
	r >> a;
}	
	
	lx_real abs  (const lx_complex& a) throw()
{ return sqrtx2y2(a.re,a.im); }
	lx_real abs2 (const lx_complex& a) throw()
	{ return a.re*a.re + a.im*a.im; }

// -----------------------------------------------------------------------------
// --------- Elementary functions related to type lx_complex -------------------
// -----------------------------------------------------------------------------

lx_complex sqr(const lx_complex& z) throw() { return (z*z); }
lx_complex sqrt(const lx_complex& z) throw()
{ return mid(sqrt(lx_cinterval(z))); }
	
lx_complex sqrt(const lx_complex& z,int n) throw()
{ return mid( sqrt(lx_cinterval(z),n) ); }
	
lx_complex exp(const lx_complex& z) throw()
{ return mid(exp(lx_cinterval(z))); }
	
lx_complex exp2(const lx_complex& z) throw()
{ return mid(exp2(lx_cinterval(z))); }

lx_complex exp10(const lx_complex& z) throw()
{ return mid(exp10(lx_cinterval(z))); }

lx_complex sin(const lx_complex& z) throw()
{ return mid(sin(lx_cinterval(z))); }

lx_complex cos(const lx_complex& z) throw()
{ return mid(cos(lx_cinterval(z))); }

lx_complex tan(const lx_complex& z) throw()
{ return mid(tan(lx_cinterval(z))); }

lx_complex cot(const lx_complex& z) throw()
{ return mid(cot(lx_cinterval(z))); }	

lx_complex asin(const lx_complex& z) throw()
{ return mid(asin(lx_cinterval(z))); }

lx_complex acos(const lx_complex& z) throw()
{ return mid(acos(lx_cinterval(z))); }

lx_complex atan(const lx_complex& z) throw()
{ return mid(atan(lx_cinterval(z))); }

lx_complex acot(const lx_complex& z) throw()
{ return mid(acot(lx_cinterval(z))); }

lx_complex sinh(const lx_complex& z) throw()
{ return mid(sinh(lx_cinterval(z))); }

lx_complex cosh(const lx_complex& z) throw()
{ return mid(cosh(lx_cinterval(z))); }

lx_complex tanh(const lx_complex& z) throw()
{ return mid(tanh(lx_cinterval(z))); }

lx_complex coth(const lx_complex& z) throw()
{ return mid(coth(lx_cinterval(z))); }

lx_complex asinh(const lx_complex& z) throw()
{ return mid(asinh(lx_cinterval(z))); }

lx_complex acosh(const lx_complex& z) throw()
{ return mid(acosh(lx_cinterval(z))); }

lx_complex atanh(const lx_complex& z) throw()
{ return mid(atanh(lx_cinterval(z))); }

lx_complex acoth(const lx_complex& z) throw()
{ return mid(acoth(lx_cinterval(z))); }

// sqrt_all(c) computes a list of 2 values for all square roots of c
std::list<lx_complex> sqrt_all( const lx_complex& c )
{ 
	lx_complex lc;
	lc = sqrt(c);

	std::list<lx_complex> res;
	res.push_back(  lc );
	res.push_back( -lc );

	return res;
} 	// end sqrt_all
	
lx_real arg(const lx_complex& z) throw()
{ return mid(arg(lx_cinterval(z))); }

lx_real Arg(const lx_complex& z) throw()
{ return mid(Arg(lx_cinterval(z))); }

std::list<lx_complex> sqrt_all( const lx_complex& z, int n )
														//
//  sqrt_all(z,n) computes a list of n values approximating all n-th roots of z
														//
//  For n >=3, computing the optimal approximations is very expensive
//  and thus not considered cost-effective.
														//
//  Hence, the polar form is used to calculate sqrt_all(z,n)
														//
{
	std::list<lx_complex> res;

	if( n == 0 )
	{
		res.push_back( lx_complex(l_real(1),l_real(0)) );
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
		lx_real
				arg_z = arg( z ), root_abs_z = sqrt( abs( z ), n );

		for(int k = 0; k < n; k++)
		{
			lx_real arg_k = ( arg_z + 2 * k * Pi_lx_real() ) / n;
			res.push_back( lx_complex( root_abs_z * cos( arg_k ),
					root_abs_z * sin( arg_k ) ) );
		}
		return res;
	}
}
												//
//-- end sqrt_all -------------------------------------------------------------

lx_complex ln(const lx_complex& z) throw()
{ return mid(ln(lx_cinterval(z))); }

lx_complex log2(const lx_complex& z) throw()
{ return mid(log2(lx_cinterval(z))); }
lx_complex log10(const lx_complex& z) throw()
{ return mid(log10(lx_cinterval(z))); }

lx_complex power_fast(const lx_complex& z, const real& n) throw()
{
	if( n == 0 )
		return lx_complex(l_real(1),l_real(0));
	else 
		if( n == 1 ) return z;
		else 
			if( n == -1 ) return 1 / z;
			else 
				if( n == 2 ) return sqr(z);
				else
				{
					lx_real abs_z = abs(z);
					if( ((n < 0) && (abs_z == 0.0)) || !(Is_Integer(n)))
	   			//  z contains 0
						cxscthrow (STD_FKT_OUT_OF_DEF(
				"lx_complex power_fast(const lx_complex& z, const real& n ); z = 0 or n is not integer."));
					if( abs_z == 0.0 )
						return lx_complex(l_real(0),l_real(0));
					else
					{
						lx_real arg_z = arg(z);
						lx_real abs_z_n = exp( n * ln( abs_z ) );

						return lx_complex( abs_z_n * cos( n * arg_z ),
												 abs_z_n * sin( n * arg_z ) );
					}
				}
} // End: power_fast

lx_complex power(const lx_complex& x, const real& n) throw()
{
	if( !(Is_Integer(n)) )
	//  n is not an integer
	cxscthrow(STD_FKT_OUT_OF_DEF(
		"lx_complex power(const lx_complex& z, const real& n); n is not integer."));
	
	real zhi(2.0), N(n), r;
	lx_complex one = lx_complex(l_real(1),l_real(0));
	double dbl;
	lx_complex y, neu, X(x);
	
	if (x == one) y = x;
	else 
		if (N == 0.0) y = one;
		else 
		{
			if (N == 1) y = x;
			else 
				if (N == 2) y = sqr(x);
				else 
				{
					if (N < 0) 
					{
						X = 1.0/X;
						N = -N;
					}
		    		// Initialisierung
					if ( !Is_Integer(N/2) )
						y = X;
					else y = one;
					neu = sqr(X);   // neu = X*X;
					do 
					{
						dbl = _double(N/zhi);
						dbl = floor(dbl);
						r = (real) dbl;
						if ( !Is_Integer( r/2 ) )
							y *= neu;
						zhi += zhi;
						if (zhi <= N) 
							neu = sqr(neu); // neu = neu * neu;
					} while (zhi <= N);
				}
		}
	return y;
} // end power(z,n)

lx_complex pow(const lx_complex& z, const lx_real& p) throw()
{ return mid( pow( lx_cinterval(z) , lx_interval(p) ) ); }

lx_complex pow(const lx_complex& z, const lx_complex& p) throw()	
{ return mid( pow( lx_cinterval(z) , lx_cinterval(p) ) ); }

lx_complex sqrt1px2(const lx_complex& z) throw()
{ return mid(sqrt1px2(lx_cinterval(z))); }

lx_complex sqrt1mx2(const lx_complex& z) throw()
{ return mid(sqrt1mx2(lx_cinterval(z))); }

lx_complex sqrtx2m1(const lx_complex& z) throw()
{ return mid(sqrtx2m1(lx_cinterval(z))); }

lx_complex sqrtp1m1(const lx_complex& z) throw()
{ return mid(sqrtp1m1(lx_cinterval(z))); }

lx_complex expm1(const lx_complex& z) throw()
{ return mid(expm1(lx_cinterval(z))); }

lx_complex lnp1(const lx_complex& z) throw()
{ return mid(lnp1(lx_cinterval(z))); }

} // namespace cxsc
