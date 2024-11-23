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

/* CVS $Id: dotk.inl,v 1.24 2014/01/30 17:23:45 cxsc Exp $ */

#ifndef _CXSC_DOTK_INL_INCLUDED
#define _CXSC_DOTK_INL_INCLUDED

#if !ROUND_ASM && !ROUND_C99_SAVE && !ROUND_C96_SAVE && !ROUND_C99_QUICK && !ROUND_C96_QUICK && !_MSC_VER
#define _CXSC_DOTK_ROUND_SOFT
#endif

#include "real.hpp"
#include "rmatrix.hpp"

#ifndef _CXSC_DOTK_ROUND_SOFT
#ifdef _MSC_VER
#include <xmmintrin.h>
#else
#include <fenv.h>
#endif
#endif

#ifdef CXSC_USE_FMA
#include <cmath>
#endif


namespace cxsc {

#if !defined(_CXSC_DOTK_ROUND_SOFT) && !defined(CXSC_USE_BLAS)
	/*
	   Sets the rounding mode according to the parameter r:
	   r=-1: Round downwards
           r=0 : Round to nearest
           r=1 : Round upwards
	   r=2 : Round toward zero
	*/
#ifndef _MSC_VER
	static inline void setround(int r) {
		switch(r) {
			case -1:
				fesetround(FE_DOWNWARD);
				break;
			case 0:
				fesetround(FE_TONEAREST);
				break;
			case 1:
				fesetround(FE_UPWARD);
				break;
			case 2: 
				fesetround(FE_TOWARDZERO);
				break;
			default:
				fesetround(FE_TONEAREST);
		}
	}

#else

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
	

	/*
	   Return an int value corresponding to the current rounding mode:
	   -1: Round downwards
           0 : Round to nearest
           1 : Round upwards
	   2 : Round toward zero
	*/
	
#ifndef _MSC_VER
	static inline int getround() {
		switch(fegetround()) {
			case FE_DOWNWARD:
				return -1;
				break;
			case FE_TONEAREST:
				return 0;
				break;
			case FE_UPWARD:
				return 1;
				break;
			case FE_TOWARDZERO:
				return 2;
				break;
			default:
				return 0;
		}	
	}

#else

#include <xmmintrin.h>
	static inline int getround() {
		switch(_MM_GET_ROUNDING_MODE()) {
			case _MM_ROUND_DOWN:
				return -1;
				break;
			case _MM_ROUND_NEAREST:
				return 0;
				break;
			case _MM_ROUND_UP:
				return 1;
				break;
			case _MM_ROUND_TOWARD_ZERO:
				return 2;
				break;
			default:
				return 0;
		}	
	}
#endif
	
#endif //ifndef CXSC_ROUND_SOFT

	
	/*
	Converts the floating point sum of two floating point variables x and y into 
	a new sum of two floating point values a and b (x+y=a+b in the set of real 
	numbers), where a is the result of the floating point summation and b is a
	correction term ("the error of the floating point result")
	
	\param x First addend
	\param y Second addend
	\param a The result of the floating point summation
	\param b The error of the floating point summation
	*/ 
	static inline void TwoSum(const real x, const real y, real &a, real &b)   {
		a = x + y;
		b = a - x;
		b = ((x - (a - b)) + (y - b));
	}
	
	/*
	Splits a floating point number x into two floating point number a and b
	with x=a+b, so that the mantissa of a and b are not overlapping and |y|<=|x|
	
	\param x Floating point number to be split
	\param x_h First part of x
	\param x_t Second part of x
	*/
	static inline void Split(const real &x, real &x_h, real &x_t)   {
		x_t = Factor * x;
		x_h = x_t - (x_t - x);
		x_t = x - x_h;
	}
	
	/*
	Computes the product of two floating point numbers as well as the error
	occuring during multiplication. A product x*y of two floating point numbers
	is transformed into the sum a+b of two other floating point numbers, where
	a is the result of the ordinary floating point multiplication and b
	is a correction value depicting the error.
	
	\param x First factor
	\param y Second factor
	\param a Result of the floating point multilpication
	\param b Error of the floating point multiplication
	*/
	static inline void TwoProduct(const real &x, const real &y, real &a, real &b) {
#ifdef CXSC_USE_FMA
#ifdef CXSC_PPC64
		a = x * y;
		double r;
		asm volatile ("fnmsub %0, %1, %2, %3\n" : "=f" (r) : "f"(_double(x)), "f" (_double(y)), "f"(_double(a)) );
		b = -r;
#else
		a = x * y;
		b = fma(_double(x), _double(y), _double(-a));
#endif
#else
		real x1,x2,y1,y2;
		a = x * y;
		Split(x,x1,x2);
		Split(y,y1,y2);
		b = x2 * y2 - (((a - x1 * y1) - x2 * y1) - x1 * y2);
#endif
	} 

	/*
	Processes the error values according to the desired precision and 
	computes the final result of the DotK algorithm for K>=3 as well as an error
	bound.
	
	\param p Array consisting of correction values and floating point result of
		the dot product
	\param n Dimension of vector
	\param k Runs of SumK algorithm
	\param err Error bound
	\param val The final result of the Dotk-Algorithm stored in an accumulator
	*/
	static inline void SumK(real* p, int n, int k, real &err, dotprecision& val) {
		real corr = 0.0, tmperr = 0.0, res = 0.0;
		val += p[n-1];
		res += p[n-1];
		
		for(int j=1 ; j<k ; j++) {
			for(int i=1 ; i<n-1 ; i++) 
				TwoSum(p[i],p[i-1],p[i],p[i-1]);    

			val += p[n-2];
			res += p[n-2];
			p[n-2] = 0.0;
		}
		
		for(int i=0 ; i<n-2 ; i++) {
			corr += p[i];
			tmperr += abs(p[i]);
		}
	

		real alpha, delta, error;
	
		delta = (n*Epsilon) / (1.0-2*n*Epsilon);
		alpha = (Epsilon*abs(res)) + (delta*tmperr+3*MinReal/Epsilon);
		error = alpha / (1.0 - 2*Epsilon);

		err = addu(err, error);
		
		val += corr;
	}


	/*
	Processes the error values according to the desired precision and 
	computes the final result of the DotK algorithm for K>=3. This version does
	NOT compute an error bound.
	
	\param p Array consisting of correction values and floating point result of
		the dot product
	\param n Dimension of vector
	\param k Runs of SumK algorithm
	\param val The final result of the Dotk-Algorithm stored in an accumulator
	*/
	template<typename T>
	static inline void SumK_NoErr(real* p, int n, int k, T& val) {
		real corr = 0.0;
		val += p[n-1];
		
		for(int j=1 ; j<k ; j++) {
			for(int i=1 ; i<n-1 ; i++) 
				TwoSum(p[i],p[i-1],p[i],p[i-1]);    

			val += p[n-2];
			p[n-2] = 0.0;
		}
		
		for(int i=0 ; i<n-2 ; i++) {
			corr += p[i];
		}
		
		val += corr;
	}


	/*
	Computes the dot product of two real vectors in k-fold double precision.
        Algorithms used for values of K (k is member of dotpecision object):
        k=0: Use accumulator (maximum accuracy)
        k=1: Use floating point computations (error bound computed via rounding mode switch)
        k>=2: Use DotK-algorithm
        The final result with an error bound is stored in an accumlator.
	
	\param val Accumulator storing the final result
	\param x First vector of dot product
	\param y Second vector of dot product
	*/
	template<typename S, typename T>
	inline void addDot(dotprecision &val, const S &x, const T &y) {
		int n = Ub(x)-Lb(x)+1;
		int lb1 = Lb(x);
		int lb2 = Lb(y);
		real res,err=0.0;

#ifndef _CXSC_DOTK_ROUND_SOFT
		int rnd;
		//Check rounding mode
		if((rnd=getround()) != 0) {
			setround(0);
		}
#endif
	
		if(val.k == 0) { //use accumulator
	
			for(int i=1 ; i<=n ; i++)
				accumulate(val, x[i+lb1-1], y[i+lb2-1]);
	
		} else if(val.k == 1) { //use floating point
	
			real resd = 0.0, resu = 0.0;
	
#ifndef _CXSC_DOTK_ROUND_SOFT
			setround(-1);
			for(int i=1 ; i<=n ; i++)
				resd += x[i+lb1-1] * y[i+lb2-1];
	
			setround(1);
			for(int i=1 ; i<=n ; i++)
				resu += x[i+lb1-1] * y[i+lb2-1];
	
			setround(0);
			res = resd+(resu-resd)*0.5;

                        setround(1);
			val.err += (resu-res);
#else
			for(int i=1 ; i<=n ; i++)
				resd = addd(resd, muld(x[i+lb1-1], y[i+lb2-1]) );
	
			for(int i=1 ; i<=n ; i++)
				resu = addu(resu, mulu(x[i+lb1-1], y[i+lb2-1]) );
	
			res = resd+(resu-resd)*0.5;

			val.err = addu(val.err, subu(resu,res) );			
#endif
			val += res;
	
		} else if(val.k == 2) { //use DotK optimized for K=2
	
			real p, s, h, r, q, t;
		
			TwoProduct(x[lb1],y[lb2],p,s);
	
			err += abs(s);
	
			for(int i=2 ; i<=n ; i++) {
				TwoProduct(x[lb1+i-1],y[lb2+i-1],h,r);
				TwoSum(p,h,p,q);
				t = q + r;
				s += t;
				err += abs(t);
			}
	
			val += p;
			val += s;

			res = p+s;
			real alpha, delta, error;
	
			delta = (n*Epsilon) / (1.0-2*n*Epsilon);
			alpha = (Epsilon*abs(res)) + (delta*err+3*MinReal/Epsilon);
			error = alpha / (1.0 - 2*Epsilon);

			val.err = addu(val.err, error);
	
	
		} else { //use DotK
	
			real r = 0, h;
			real* t = new real[2*n];
		
			for(int i=1 ; i<=n ; i++) {
				TwoProduct(x[lb1+i-1], y[lb2+i-1], h, t[i-1]);
				TwoSum(r, h, r, t[n+i-2]);
			}
	
			t[2*n-1] = r;
			SumK(t, 2*n, val.k-1, err, val);
	
			val.err = addu(val.err, err);
	
			delete[] t;
	
		}

#ifndef _CXSC_DOTK_ROUND_SOFT
		//Reset rounding mode to former value
		setround(rnd);
#endif
	
	}


	/*
	Computes the dot product of two real vectors in k-fold double precision.
        Algorithms used for values of K (k is member of dotpecision object):
        k=0: Use accumulator (maximum accuracy)
        k=1: Use floating point computations 
        k>=2: Use DotK-algorithm
        This version does not compute an error bound and is used for computations
	in the operators.
	
	\param val Accumulator storing the final result
	\param x First vector of dot product
	\param y Second vector of dot product
	*/
	template<typename S, typename T>
	inline void addDot_op(dotprecision &val, const S &x, const T &y) {
		int n = Ub(x)-Lb(x)+1;
		int lb1 = Lb(x);
		int lb2 = Lb(y);
		real res = 0.0;

#ifndef _CXSC_DOTK_ROUND_SOFT
		int rnd;
		//Check rounding mode
		if((rnd=getround()) != 0) {
			setround(0);
		}
#endif

		if(val.k == 0) { //Use accumulator
	
			for(int i=1 ; i<=n ; i++)
				accumulate(val, x[i+lb1-1], y[i+lb2-1]);
	
		} else if(val.k == 1) { //Use floating point
	
			for(int i=1 ; i<=n ; i++)
				res += x[i+lb1-1] * y[i+lb2-1];
	
			val += res;
	
		} else if(val.k == 2) { //Use optimized DotK
	
			real p, s, h, r, q, t;
	
			TwoProduct(x[lb1],y[lb2],p,s);
		
			for(int i=2 ; i<=n ; i++) {
				TwoProduct(x[lb1+i-1],y[lb2+i-1],h,r);
				TwoSum(p,h,p,q);
				t = q + r;
				s += t;
			}
	
			val += p;
			val += s;
	
		} else { //Use DotK
	
			real r = 0, h;
			real* t = new real[2*n];
	
			for(int i=1 ; i<=n ; i++) {
				TwoProduct(x[lb1+i-1], y[lb2+i-1], h, t[i-1]);
				TwoSum(r, h, r, t[n+i-2]);
			}
	
			t[2*n-1] = r;
			SumK_NoErr(t, 2*n, val.k-1, val);
	
			delete[] t;
	
		}

#ifndef _CXSC_DOTK_ROUND_SOFT
		//Reset rounding mode to former value
		setround(rnd);
#endif
	
	}


	/*
	Computes the sum of the elements of a real vector in k-fold double precision.
        Algorithms used for values of K (k is member of dotpecision object):
        k=0: Use accumulator (maximum accuracy)
        k=1: Use floating point computations (error bound computed via rounding mode switch)
        k>=2: Use DotK-algorithm
        The final result with an error bound is stored in an accumlator.
	
	\param val Accumulator storing the final result
	\param x First vector of dot product
	\param y Second vector of dot product
	*/
	template<typename S>
	inline void addSum(dotprecision &val, const S &x) {
		int n = Ub(x)-Lb(x)+1;
		int lb = Lb(x);
		real res,err=0.0;

#ifndef _CXSC_DOTK_ROUND_SOFT
		int rnd;
		//Check rounding mode
		if((rnd=getround()) != 0) {
			setround(0);
		}
#endif
	
		if(val.k == 0) { //use accumulator
	
			for(int i=0 ; i<n ; i++)
				val += x[i+lb];
	
		} else if(val.k == 1) { //use floating point
	
			real resd = x[lb], resu = x[lb];
	
#ifndef _CXSC_DOTK_ROUND_SOFT
			setround(-1);
			for(int i=1 ; i<n ; i++)
				resd += x[i+lb];
	
			setround(1);
			for(int i=1 ; i<n ; i++)
				resu += x[i+lb];
	
			setround(0);
			res = resd+(resu-resd)*0.5;

                        setround(1);
			val.err += (resu-res);
#else
			for(int i=1 ; i<n ; i++)
				resd = addd(resd,x[i+lb]);
	
			for(int i=1 ; i<n ; i++)
				resu = addu(resu,x[i+lb]);
	
			res = resd+(resu-resd)*0.5;

			val.err = addu(val.err,subu(resu,res));
#endif

			val += res;
	
		} else if(val.k == 2) { //use DotK optimized for K=2
	
			real p=x[lb], s=0.0, q;
		
			for(int i=1 ; i<n ; i++) {
				TwoSum(p,x[i+lb],p,q);
				s += q;
				err += abs(q);
			}
	
			val += p;
			val += s;

			res = p+s;
			real alpha, delta, error;
	
			delta = (n*Epsilon) / (1.0-2*n*Epsilon);
			alpha = (Epsilon*abs(res)) + (delta*err+3*MinReal/Epsilon);
			error = alpha / (1.0 - 2*Epsilon);

			val.err = addu(val.err, error);
	
	
		} else { //use DotK
	
			real r = x[lb];
			real* t = new real[n];
		
			for(int i=1 ; i<n ; i++) {
				TwoSum(r, x[i+lb], r, t[i-1]);
			}
	
			t[n-1] = r;
			SumK(t, n, val.k-1, err, val);
	
			val.err = addu(val.err, err);
	
			delete[] t;
	
		}

#ifndef _CXSC_DOTK_ROUND_SOFT
		//Reset rounding mode to former value
		setround(rnd);
#endif
	
	}


	

}

#endif // _CXSC_DOTK_INL_INCLUDED
