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

/* CVS $Id: cdotk.inl,v 1.16 2014/01/30 17:23:43 cxsc Exp $ */

#ifndef _CXSC_CDOTK_INL_INCLUDED
#define _CXSC_CDOTK_INL_INCLUDED

#include "dotk.inl"

namespace cxsc {

	/*
	Computes the dot product of two complex vectors in k-fold double precision.
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
	static inline void addDot(cdotprecision &val, const S &x, const T &y) {
		int n = VecLen(x);
		int lb1 = Lb(x);
		int lb2 = Lb(y);
		int k = val.get_k();

#ifndef _CXSC_DOTK_ROUND_SOFT
		int rnd;
		//Check rounding mode
		if((rnd=getround()) != 0) {
			setround(0);
		}
#endif
	
		if(k == 0) { //Use accumulator

			for(int i=0 ; i<n ; i++)
				accumulate(val, x[i+lb1], y[i+lb2]);

		} else if(k == 1) { //Use floating point

			real res_re_d = 0.0, res_re_u = 0.0;
			real res_im_d = 0.0, res_im_u = 0.0;
	
#ifndef _CXSC_DOTK_ROUND_SOFT
			setround(-1);
			for(int i=0 ; i<n ; i++) {
				res_re_d += Re(x[i+lb1]) * Re(y[i+lb2]) - Im(x[i+lb1]) * Im(y[i+lb2]);
				res_im_d += Re(x[i+lb1]) * Im(y[i+lb2]) + Im(x[i+lb1]) * Re(y[i+lb2]);
			}
	
			setround(1);
			for(int i=0 ; i<n ; i++) {
				res_re_u += Re(x[i+lb1]) * Re(y[i+lb2]) - Im(x[i+lb1]) * Im(y[i+lb2]);
				res_im_u += Re(x[i+lb1]) * Im(y[i+lb2]) + Im(x[i+lb1]) * Re(y[i+lb2]);
			}
	
			setround(0);
			real res_re = res_re_d+(res_re_u-res_re_d)*0.5;
			real res_im = res_im_d+(res_im_u-res_im_d)*0.5;

                        setround(1);
			Re(val).err += res_re_u - res_re;
			Im(val).err += res_im_u - res_im;

#else

			for(int i=0 ; i<n ; i++) {
				res_re_d = addd( res_re_d, subd( muld(Re(x[i+lb1]),Re(y[i+lb2])) , muld(Im(x[i+lb1]),Im(y[i+lb2])) ) );
				res_im_d = addd( res_im_d, addd( muld(Re(x[i+lb1]),Im(y[i+lb2])) , muld(Im(x[i+lb1]),Re(y[i+lb2])) ) );
			}
	
			for(int i=0 ; i<n ; i++) {
				res_re_u = addu( res_re_u, subu( mulu(Re(x[i+lb1]),Re(y[i+lb2])) , mulu(Im(x[i+lb1]),Im(y[i+lb2])) ) );
				res_im_u = addu( res_im_u, addu( mulu(Re(x[i+lb1]),Im(y[i+lb2])) , mulu(Im(x[i+lb1]),Re(y[i+lb2])) ) );
			}
	
			real res_re = res_re_d+(res_re_u-res_re_d)*0.5;
			real res_im = res_im_d+(res_im_u-res_im_d)*0.5;

			Re(val).err = addu(Re(val).err, subu(res_re_u,res_re) );
			Im(val).err = addu(Im(val).err, subu(res_im_u,res_im) );

#endif

			val += complex(res_re,res_im);
		} else if(k == 2) { //Use optimized DotK

			real p_re = 0.0, p_im = 0.0;
			real s_re = 0.0, s_im = 0.0, h, r, q, t;
			real err_re = 0.0, err_im = 0.0;
	
			for(int i=0 ; i<n ; i++) {
				TwoProduct(Re(x[lb1+i]), Re(y[lb2+i]), h, r);
				TwoSum(p_re, h, p_re, q);
				t = q + r;
				s_re += t;
				err_re += abs(t);
	
				TwoProduct(-Im(x[lb1+i]), Im(y[lb2+i]), h, r);
				TwoSum(p_re, h, p_re, q);
				t = q + r;
				s_re += t;
				err_re += abs(t);
	
				TwoProduct(Re(x[lb1+i]), Im(y[lb2+i]), h, r);
				TwoSum(p_im, h, p_im, q);
				t = q + r;
				s_im += t;
				err_im += abs(t);
	
				TwoProduct(Im(x[lb1+i]), Re(y[lb2+i]), h, r);
				TwoSum(p_im, h, p_im, q);
				t = q + r;
				s_im += t;
				err_im += abs(t);
			}

			val += complex(p_re, p_im);
			val += complex(s_re, s_im);

			//Compute error bound
 			complex res = complex(p_re+s_re, p_im+s_im);
			real alpha_re, alpha_im, delta, error_re, error_im;
			delta = (2*n*Epsilon) / (1.0-4*n*Epsilon);

			alpha_re = (Epsilon*abs(Re(res))) + (delta*err_re+3*MinReal/Epsilon);
			error_re = alpha_re / (1.0 - 2*Epsilon);      
			Re(val).err = addu(Re(val).err, error_re);

			alpha_im = (Epsilon*abs(Im(res))) + (delta*err_im+3*MinReal/Epsilon);      
			error_im = alpha_im / (1.0 - 2*Epsilon);            
			Im(val).err = addu(Im(val).err, error_im);  

			
		} else { //Use DotK

			real r_re = 0, r_im=0, h;
			real* t_re = new real[4*n];
			real* t_im = new real[4*n];
		
			for(int i=1 ; i<=n ; i++) {
				int ind2 = 2*i;
				int ind1 = ind2 - 1;
		
				TwoProduct(Re(x[lb1+i-1]), Re(y[lb2+i-1]), h, t_re[ind1-1]);
				TwoSum(r_re, h, r_re, t_re[2*n+ind1-2]);
		
				TwoProduct(-Im(x[lb1+i-1]), Im(y[lb2+i-1]), h, t_re[ind2-1]);
				TwoSum(r_re, h, r_re, t_re[2*n+ind2-2]);
	
				TwoProduct(Re(x[lb1+i-1]), Im(y[lb2+i-1]), h, t_im[ind1-1]);
				TwoSum(r_im, h, r_im, t_im[2*n+ind1-2]);
	
				TwoProduct(Im(x[lb1+i-1]), Re(y[lb2+i-1]), h, t_im[ind2-1]);
				TwoSum(r_im, h, r_im, t_im[2*n+ind2-2]);
			}
	
			t_re[4*n-1] = r_re;
			t_im[4*n-1] = r_im;

			SumK(t_re, 4*n, k-1, Re(val).err, Re(val));
			SumK(t_im, 4*n, k-1, Im(val).err, Im(val));
	
			delete [] t_re;
			delete [] t_im;
		}

		//Reset rounding mode
#ifndef _CXSC_DOTK_ROUND_SOFT
		setround(rnd);
#endif
	}


	/*
	Computes the dot product of two complex vectors in k-fold double precision.
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
	static inline void addDot_op(cdotprecision &val, const S &x, const T &y) {
		int n = VecLen(x);
		int lb1 = Lb(x);
		int lb2 = Lb(y);
		int k = val.get_k();

#ifndef _CXSC_DOTK_ROUND_SOFT
		int rnd;
		//Check rounding mode
		if((rnd=getround()) != 0) {
			setround(0);
		}
#endif
	
		if(k == 0) { //Use accumulator

			for(int i=0 ; i<n ; i++)
				accumulate(val, x[i+lb1], y[i+lb2]);

		} else if(k == 1) { //Use floating point

			real res_re = 0.0;
			real res_im = 0.0;
	
			for(int i=0 ; i<n ; i++) {
				res_re += Re(x[i+lb1]) * Re(y[i+lb2]) - Im(x[i+lb1]) * Im(y[i+lb2]);
				res_im += Re(x[i+lb1]) * Im(y[i+lb2]) + Im(x[i+lb1]) * Re(y[i+lb2]);
			}

			val += complex(res_re,res_im);

		} else if(k == 2) { //Use optimized DotK

			real p_re = 0, p_im = 0;
			real s_re = 0, s_im = 0, h, r, q, t;
	
			for(int i=0 ; i<n ; i++) {
				TwoProduct(Re(x[lb1+i]), Re(y[lb2+i]), h, r);
				TwoSum(p_re, h, p_re, q);
				t = q + r;
				s_re += t;
	
				TwoProduct(-Im(x[lb1+i]), Im(y[lb2+i]), h, r);
				TwoSum(p_re, h, p_re, q);
				t = q + r;
				s_re += t;
	
				TwoProduct(Re(x[lb1+i]), Im(y[lb2+i]), h, r);
				TwoSum(p_im, h, p_im, q);
				t = q + r;
				s_im += t;
	
				TwoProduct(Im(x[lb1+i]), Re(y[lb2+i]), h, r);
				TwoSum(p_im, h, p_im, q);
				t = q + r;
				s_im += t;
			}

			val += complex(p_re, p_im);
			val += complex(s_re, s_im);
			
		} else { //Use DotK

			real r_re = 0, r_im=0, h;
			real* t_re = new real[4*n];
			real* t_im = new real[4*n];
		
			for(int i=1 ; i<=n ; i++) {
				int ind2 = 2*i;
				int ind1 = ind2 - 1;
		
				TwoProduct(Re(x[lb1+i-1]), Re(y[lb2+i-1]), h, t_re[ind1-1]);
				TwoSum(r_re, h, r_re, t_re[2*n+ind1-2]);
		
				TwoProduct(-Im(x[lb1+i-1]), Im(y[lb2+i-1]), h, t_re[ind2-1]);
				TwoSum(r_re, h, r_re, t_re[2*n+ind2-2]);
	
				TwoProduct(Re(x[lb1+i-1]), Im(y[lb2+i-1]), h, t_im[ind1-1]);
				TwoSum(r_im, h, r_im, t_im[2*n+ind1-2]);
	
				TwoProduct(Im(x[lb1+i-1]), Re(y[lb2+i-1]), h, t_im[ind2-1]);
				TwoSum(r_im, h, r_im, t_im[2*n+ind2-2]);
			}
	
			t_re[4*n-1] = r_re;
			t_im[4*n-1] = r_im;

			SumK_NoErr(t_re, 4*n, k-1, Re(val));
			SumK_NoErr(t_im, 4*n, k-1, Im(val));
	
			delete [] t_re;
			delete [] t_im;
		}

#ifndef _CXSC_DOTK_ROUND_SOFT
		//Reset rounding mode
		setround(rnd);
#endif

	}

}

#endif // _CXSC_CDOTK_INL_INCLUDED
