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

/* CVS $Id: l_cmath.cpp,v 1.11 2014/01/30 17:23:46 cxsc Exp $ */

#include "l_cmath.hpp"

namespace cxsc {
	
	l_complex sqrt(const l_complex& z) throw()
// Computation of sqrt(z); stagprec <= stagmax=19 determines the maximum 
// accuracy of about 16*19 = 304 decimal digits.
// The branch cut of this sqrt-function is the negative imaginary axis.
// stagmax=19 is predetermined by the internal used function 
// l_real sqrt(const l_real&), which has its own maximum precision of 
// stagprec=19;
	{
		l_real x=Re(z), y=Im(z), w;
		if ( zero_(x) && zero_(y) ) return l_complex(0,0);
    // Now z != 0
		int stagsave = stagprec, stagmax = 19;
		if (stagprec > stagmax) stagprec = stagmax;
    // stagprec>19 makes no sense, because stagmax(sqrt(...))=19;
		int ex = expo(x[1]), exy = expo(y[1]);
		if (exy > ex) ex = exy; // ex: maximum of the exponents;
		ex = 400-ex; // ex: optimal scaling factor
		if (ex%2) ex--;  // ex is now even;
		times2pown(x,ex);  times2pown(y,ex); // scaling with 2^ex
		bool neg = sign(x[1]) < 0;
		if (neg) x = -x;
		w = abs( l_complex(x,y) ) + x;
		times2pown(w,1);
		w = sqrt(w);
		if (neg)
		{
			x = abs(y)/w;
			y = sign(y[1]) < 0 ? -w : w;
			times2pown(y,-1); 
		} else
		{
			x = w;
			times2pown(x,-1); 
			y /= w;
		}
    ex /= 2; // Backscaling of the current result with 2^(-ex):
	 times2pown(x,-ex);  times2pown(y,-ex);
	 stagprec = stagsave; // restore old value of the stagprec variable
	 return l_complex(x,y);
	} // sqrt
	
	l_complex sqrtp1m1(const l_complex& z) throw()
	{ return mid(sqrtp1m1(l_cinterval(z))); }
	
	l_complex sqrt1px2(const l_complex& z) throw()
	{ return mid(sqrt1px2(l_cinterval(z))); }
	
	l_complex sqrtx2m1(const l_complex& z) throw()
	{ return mid(sqrtx2m1(l_cinterval(z))); }
	
	l_complex sqrt1mx2(const l_complex& z) throw()
	{ return mid(sqrt1mx2(l_cinterval(z))); }

	l_complex exp(const l_complex& z) throw()
	{ return mid(exp(l_cinterval(z))); }
	
	l_complex expm1(const l_complex& z) throw()
	{ return mid(expm1(l_cinterval(z))); }
	
	l_complex exp2(const l_complex& z) throw()
	{ return mid(exp2(l_cinterval(z))); }
	
	l_complex exp10(const l_complex& z) throw()
	{ return mid(exp10(l_cinterval(z))); }

	l_complex sin(const l_complex& z) throw()
	{ return mid(sin(l_cinterval(z))); }

	l_complex cos(const l_complex& z) throw()
	{ return mid(cos(l_cinterval(z))); }

	l_complex tan(const l_complex& z) throw()
	{ return mid(tan(l_cinterval(z))); }

	l_complex cot(const l_complex& z) throw()
	{ return mid(cot(l_cinterval(z))); }	
	
	l_complex asin(const l_complex& z) throw()
	{ return mid(asin(l_cinterval(z))); }

	l_complex acos(const l_complex& z) throw()
	{ return mid(acos(l_cinterval(z))); }

	l_complex atan(const l_complex& z) throw()
	{ return mid(atan(l_cinterval(z))); }

	l_complex acot(const l_complex& z) throw()
	{ return mid(acot(l_cinterval(z))); }	
	
	l_complex sinh(const l_complex& z) throw()
	{ return mid(sinh(l_cinterval(z))); }

	l_complex cosh(const l_complex& z) throw()
	{ return mid(cosh(l_cinterval(z))); }

	l_complex tanh(const l_complex& z) throw()
	{ return mid(tanh(l_cinterval(z))); }

	l_complex coth(const l_complex& z) throw()
	{ return mid(coth(l_cinterval(z))); }
	
	l_complex asinh(const l_complex& z) throw()
	{ return mid(asinh(l_cinterval(z))); }

	l_complex acosh(const l_complex& z) throw()
	{ return mid(acosh(l_cinterval(z))); }

	l_complex atanh(const l_complex& z) throw()
	{ return mid(atanh(l_cinterval(z))); }

	l_complex acoth(const l_complex& z) throw()
	{ return mid(acoth(l_cinterval(z))); }

   // sqrt_all(c) computes a list of 2 values for all square roots of c
	std::list<l_complex> sqrt_all( const l_complex& c )
	{ 
		l_complex lc;
		lc = sqrt(c);
		
		std::list<l_complex> res;
		res.push_back(  lc );
		res.push_back( -lc );

		return res;
	} 	// end sqrt_all
	
	l_complex sqrt(const l_complex& z, int n) throw()
	{ return mid(sqrt(l_cinterval(z),n)); }		
			
	l_real arg(const l_complex& z) throw()
	{ return mid(arg(l_cinterval(z))); }
	
	l_real Arg(const l_complex& z) throw()
	{ return mid(Arg(l_cinterval(z))); }
	
	std::list<l_complex> sqrt_all( const l_complex& z, int n )
//
//  sqrt_all(z,n) computes a list of n values approximating all n-th roots of z
//
//  For n >=3, computing the optimal approximations is very expensive
//  and thus not considered cost-effective.
//
//  Hence, the polar form is used to calculate sqrt_all(z,n)
//
	{
		std::list<l_complex> res;

		if( n == 0 )
		{
			res.push_back( l_complex(1,0) );
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
			l_real
				arg_z = arg( z ), root_abs_z = sqrt( abs( z ), n );

			for(int k = 0; k < n; k++)
			{
				l_real arg_k = ( arg_z + 2 * k * mid(Pi_l_interval()) ) / n;

				res.push_back( l_complex( root_abs_z * cos( arg_k ),
									           root_abs_z * sin( arg_k ) ) );
			}
			return res;
		}
	}
	//
//-- end sqrt_all -------------------------------------------------------------
	
	l_complex ln(const l_complex& z) throw()
	{ return mid(ln(l_cinterval(z))); }
	
	l_complex lnp1(const l_complex& z) throw()
	{ return mid(lnp1(l_cinterval(z))); }
	
	l_complex power_fast(const l_complex& z, int n) throw()
	{
		if( n == 0 )
			return l_complex(1,0);
	   else 
			if( n == 1 ) return z;
	      else 
				if( n == -1 ) return 1 / z;
	         else 
					if( n == 2 ) return sqr(z); 
	else
	{
		l_real abs_z = abs(z);

		if( n < 0 && abs_z == 0.0 )
	   //  z contains 0
			cxscthrow (STD_FKT_OUT_OF_DEF
				("l_complex power_fast(const l_complex& z, int n ); z == 0."));
		if( abs_z == 0.0 )
			return l_complex(0,0);
		else
		{
			l_real arg_z = arg(z);
			l_real abs_z_n = exp( n * ln( abs_z ) );

			return l_complex( abs_z_n * cos( n * arg_z ),
									abs_z_n * sin( n * arg_z ) );
		}
	}
}

l_complex power(const l_complex& z, int n) throw()
{ return mid( power(l_cinterval(z),n) ); }
	
l_complex log2(const l_complex& z) throw()
{ return mid(log2(l_cinterval(z))); }

l_complex log10(const l_complex& z) throw()
{ return mid(log10(l_cinterval(z))); }
	
l_complex pow(const l_complex& z, const l_real& p) throw()
{ return mid( pow( l_cinterval(z) , l_interval(p) ) ); }
	
l_complex pow(const l_complex& z, const l_complex& p) throw()	
{ return mid( pow( l_cinterval(z) , l_cinterval(p) ) ); }



} // namespace cxsc
