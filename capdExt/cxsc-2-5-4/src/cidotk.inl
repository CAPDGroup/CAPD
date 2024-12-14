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

/* CVS $Id: cidotk.inl,v 1.18 2014/01/30 17:23:43 cxsc Exp $ */

#ifndef _CXSC_CIDOTK_INL_INCLUDED
#define _CXSC_CIDOTK_INL_INCLUDED

#include "idotk.inl"

namespace cxsc {

#ifndef _CXSC_DOTK_ROUND_SOFT
	/*
	Compute complex interval dot product in floating point using midpoint radius representation
	for the real and imaginary part
	*/
	static inline cinterval floatcidot(const civector& x, const civector& y, int n, int lb1, int lb2) {
		cvector xmid(n),xrad(n),ymid(n),yrad(n);
		real crad_re = 0.0, crad_im = 0.0;
		complex res_inf(0.0,0.0),res_sup(0.0,0.0);

		setround(1);

		for(int i=0 ; i<n ; i++) {
			xmid[i+1] = Inf(x[lb1+i]) + 0.5*(Sup(x[lb1+i])-Inf(x[lb1+i]));
			ymid[i+1] = Inf(y[lb2+i]) + 0.5*(Sup(y[lb2+i])-Inf(y[lb2+i]));				

		}

		for(int i=0 ; i<n ; i++) {
			xrad[i+1] = xmid[i+1] - Inf(x[i+lb1]);
			yrad[i+1] = ymid[i+1] - Inf(y[i+lb2]);
		}

		for(int i=0 ; i<n ; i++) {
			crad_re += Re(xrad[i+1]) * (abs(Re(ymid[i+1])) + Re(yrad[i+1])) + abs(Re(xmid[i+1])) * Re(yrad[i+1]);
			crad_re += Im(xrad[i+1]) * (abs(Im(ymid[i+1])) + Im(yrad[i+1])) + abs(Im(xmid[i+1])) * Im(yrad[i+1]);

			crad_im += Re(xrad[i+1]) * (abs(Im(ymid[i+1])) + Im(yrad[i+1])) + abs(Re(xmid[i+1])) * Im(yrad[i+1]);
			crad_im += Im(xrad[i+1]) * (abs(Re(ymid[i+1])) + Re(yrad[i+1])) + abs(Im(xmid[i+1])) * Re(yrad[i+1]);

			Re(res_sup) += Re(xmid[i+1]) * Re(ymid[i+1]) - Im(xmid[i+1]) * Im(ymid[i+1]);
			Im(res_sup) += Re(xmid[i+1]) * Im(ymid[i+1]) + Im(xmid[i+1]) * Re(ymid[i+1]);
		}

		res_sup += complex(crad_re,crad_im);

		setround(-1);
		for(int i=0 ; i<n ; i++) {
			Re(res_inf) += Re(xmid[i+1]) * Re(ymid[i+1]) - Im(xmid[i+1]) * Im(ymid[i+1]);
			Im(res_inf) += Re(xmid[i+1]) * Im(ymid[i+1]) + Im(xmid[i+1]) * Re(ymid[i+1]);
		}

		res_inf -= complex(crad_re,crad_im);

		return cinterval(res_inf,res_sup);
	}


	/*
	Compute complex interval dot product in floating point using midpoint radius representation
	for the real and imaginary part
	*/
	static inline cinterval floatcidot(const civector& x, const cvector& y, int n, int lb1, int lb2) {
		cvector xmid(n),xrad(n);
		real crad_re = 0.0, crad_im = 0.0;
		complex res_inf(0.0,0.0),res_sup(0.0,0.0);

		setround(1);

		for(int i=0 ; i<n ; i++) {
			xmid[i+1] = Inf(x[lb1+i]) + 0.5*(Sup(x[lb1+i])-Inf(x[lb1+i]));

		}

		for(int i=0 ; i<n ; i++) {
			xrad[i+1] = xmid[i+1] - Inf(x[i+lb1]);
		}

		for(int i=0 ; i<n ; i++) {
			crad_re += Re(xrad[i+1]) * abs(Re(y[i+lb2])) + Im(xrad[i+1]) * abs(Im(y[i+lb2]));
			crad_im += Im(xrad[i+1]) * abs(Re(y[i+lb2])) + Re(xrad[i+1]) * abs(Im(y[i+lb2]));
			Re(res_sup) += Re(xmid[i+1]) * Re(y[i+lb2]) - Im(xmid[i+1]) * Im(y[i+lb2]);
			Im(res_sup) += Re(xmid[i+1]) * Im(y[i+lb2]) + Im(xmid[i+1]) * Re(y[i+lb2]);
		}

		res_sup += complex(crad_re,crad_im);

		setround(-1);
		for(int i=0 ; i<n ; i++) {
			Re(res_inf) += Re(xmid[i+1]) * Re(y[i+lb2]) - Im(xmid[i+1]) * Im(y[i+lb2]);
			Im(res_inf) += Re(xmid[i+1]) * Im(y[i+lb2]) + Im(xmid[i+1]) * Re(y[i+lb2]);
		}

		res_inf -= complex(crad_re,crad_im);

		return cinterval(res_inf,res_sup);
	}

	/*
	Compute complex interval dot product in floating point using midpoint radius representation
	for the real and imaginary part
	*/
	static inline cinterval floatcidot(const cvector& x, const civector& y, int n, int lb1, int lb2) {		
		return floatcidot(y,x,n,lb2,lb1);
	}
#endif

	/*
	Computes the dot product of two complex interval vectors in k-fold double precision.
        Algorithms used for values of K (k is member of dotpecision object):
        k=0: Use accumulator (maximum accuracy)
        k=1: Use floating point computations (mid-rad)
        k>=2: Use DotK-algorithm
        The final result with an error bound is stored in an accumlator.
	
	\param val Accumulator storing the final result
	\param x First vector of dot product
	\param y Second vector of dot product
	*/
	template<typename S, typename T>
	static inline void addDot( cidotprecision &val, const S &x, const T &y ) {
		int n = VecLen(x);
		int lb1 = Lb(x);
		int lb2 = Lb(y);
		int k = val.get_k();
		real err_inf_re=0.0, err_sup_re=0.0, err_inf_im=0.0, err_sup_im=0.0;

#ifndef _CXSC_DOTK_ROUND_SOFT
		int rnd;
		//Check rounding mode
		if((rnd=getround()) != 0) {
			setround(0);
		}
#endif
	
		if(k == 0) { //Use accumulator

			for(int i=0 ; i<n ; i++) 
				accumulate(val, x[lb1+i], y[lb2+i]);

		} else if(k == 1) { //Use floating point
#ifndef _CXSC_DOTK_ROUND_SOFT
			val += floatcidot(x,y,n,lb1,lb2);
#else
			cinterval res(0.0);
			for(int i=0 ; i<n ; i++) 
				res += x[lb1+i] * y[lb2+i];
			val += res;
#endif
		} else if(k == 2) { //Use optimized DotK

			interval p_re(0,0), p_im(0,0), h;
			real s_inf_re=0, s_inf_im=0, s_sup_re=0, s_sup_im=0;
			real q_inf, q_sup, r_inf, r_sup;
			real t;
	
			for(int i=0 ; i<n ; i++) {
				TwoProduct(Re(x[lb1+i]), Re(y[lb2+i]), h, r_inf, r_sup);
				TwoSum(p_re, h, p_re, q_inf, q_sup);
	
				t = q_inf + r_inf;
				s_inf_re += t;
				err_inf_re += abs(t);
	
				t = q_sup + r_sup;
				s_sup_re += t;
				err_sup_re += abs(t);
	
	
				TwoProduct(-Im(x[lb1+i]), Im(y[lb2+i]), h, r_inf, r_sup);
				TwoSum(p_re, h, p_re, q_inf, q_sup);
	
				t = q_inf + r_inf;
				s_inf_re += t;
				err_inf_re += abs(t);
	
				t = q_sup + r_sup;
				s_sup_re += t;
				err_sup_re += abs(t);
	
	
				TwoProduct(Re(x[lb1+i]), Im(y[lb2+i]), h, r_inf, r_sup);
				TwoSum(p_im, h, p_im, q_inf, q_sup);
	
				t = q_inf + r_inf;
				s_inf_im += t;
				err_inf_im += abs(t);
	
				t = q_sup + r_sup;
				s_sup_im += t;
				err_sup_im += abs(t);
	
	
				TwoProduct(Im(x[lb1+i]), Re(y[lb2+i]), h, r_inf, r_sup);
				TwoSum(p_im, h, p_im, q_inf, q_sup);
	
				t = q_inf + r_inf;
				s_inf_im += t;
				err_inf_im += abs(t);
	
				t = q_sup + r_sup;
				s_sup_im += t;
				err_sup_im += abs(t);    
			}
		
			dotprecision& val_inf_re = InfRe(val);
			dotprecision& val_inf_im = InfIm(val);
			dotprecision& val_sup_re = SupRe(val);
			dotprecision& val_sup_im = SupIm(val);

			val_inf_re += Inf(p_re);
			val_inf_re += s_inf_re;
			val_sup_re += Sup(p_re);
			val_sup_re += s_sup_re;

			val_inf_im += Inf(p_im);
			val_inf_im += s_inf_im;
			val_sup_im += Sup(p_im);
			val_sup_im += s_sup_im;

			//Compute error bound
 			real alpha, delta, error;
			delta = (2*n*Epsilon) / (1.0-4*n*Epsilon);

			alpha = (Epsilon*abs(Inf(p_re)+s_inf_re)) + (delta*err_inf_re+3*MinReal/Epsilon);
			error = alpha / (1.0 - 2*Epsilon);      
			val_inf_re -= error;

			alpha = (Epsilon*abs(Sup(p_re)+s_sup_re)) + (delta*err_sup_re+3*MinReal/Epsilon);
			error = alpha / (1.0 - 2*Epsilon);      
			val_sup_re += error;

			alpha = (Epsilon*abs(Inf(p_im)+s_inf_im)) + (delta*err_inf_im+3*MinReal/Epsilon);
			error = alpha / (1.0 - 2*Epsilon);      
			val_inf_im -= error;

			alpha = (Epsilon*abs(Sup(p_im)+s_sup_im)) + (delta*err_sup_im+3*MinReal/Epsilon);
			error = alpha / (1.0 - 2*Epsilon);      
			val_sup_im += error;
		
		} else { //Use DotK
		
			interval r_re(0,0), r_im(0,0), h;
			real* t_inf_re = new real[4*n];
			real* t_sup_re = new real[4*n];
			real* t_inf_im = new real[4*n];
			real* t_sup_im = new real[4*n];
	
			for(int i=1 ; i<=n ; i++) {
				int ind2 = 2*i;
				int ind1 = ind2 - 1;                     
			
				TwoProduct(Re(x[lb1+i-1]), Re(y[lb2+i-1]), h, t_inf_re[ind1-1], t_sup_re[ind1-1]);
				TwoSum(r_re, h, r_re, t_inf_re[2*n+ind1-2], t_sup_re[2*n+ind1-2]);
	
				TwoProduct(-Im(x[lb1+i-1]), Im(y[lb2+i-1]), h, t_inf_re[ind2-1], t_sup_re[ind2-1]);
				TwoSum(r_re, h, r_re, t_inf_re[2*n+ind2-2], t_sup_re[2*n+ind2-2]);
	
				TwoProduct(Re(x[lb1+i-1]), Im(y[lb2+i-1]), h, t_inf_im[ind1-1], t_sup_im[ind1-1]);
				TwoSum(r_im, h, r_im, t_inf_im[2*n+ind1-2], t_sup_im[2*n+ind1-2]);
	
				TwoProduct(Im(x[lb1+i-1]), Re(y[lb2+i-1]), h, t_inf_im[ind2-1], t_sup_im[ind2-1]);
				TwoSum(r_im, h, r_im, t_inf_im[2*n+ind2-2], t_sup_im[2*n+ind2-2]);
			}              
	
			t_inf_re[4*n-1] = Inf(r_re);
			t_sup_re[4*n-1] = Sup(r_re);
			t_inf_im[4*n-1] = Inf(r_im);
			t_sup_im[4*n-1] = Sup(r_im);
		
			SumK(t_inf_re, 4*n, k-1, err_inf_re, InfRe(val));
			SumK(t_sup_re, 4*n, k-1, err_sup_re, SupRe(val));
			SumK(t_inf_im, 4*n, k-1, err_inf_im, InfIm(val));
			SumK(t_sup_im, 4*n, k-1, err_sup_im, SupIm(val));

			InfRe(val) -= err_inf_re;
			SupRe(val) += err_sup_re;
			InfIm(val) -= err_inf_im;
			SupIm(val) += err_sup_im;		

			delete [] t_inf_re;
			delete [] t_sup_re;
			delete [] t_inf_im;
			delete [] t_sup_im;
	
		}         

#ifndef _CXSC_DOTK_ROUND_SOFT
		//Reset rounding mode
		setround(rnd);
#endif
	
	}

}

#endif // _CXSC_CIDOTK_INL_INCLUDED
 
