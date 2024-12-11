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

/* CVS $Id: idotk.inl,v 1.18 2014/01/30 17:23:45 cxsc Exp $ */

#ifndef _CXSC_IDOTK_INL_INCLUDED
#define _CXSC_IDOTK_INL_INCLUDED

#include "real.hpp"
#include "dotk.inl"
#include "ivector.hpp"


namespace cxsc {

	/*
	Computes the floating point summation of two intervals as well as the error 
	occuring during this summations for the infimum and the supremum
	
	\param x First addend
	\param y Second addend
	\param a The result of the floating point summation
	\param b_i The error of the floating point summation for the infimum
	\param b_s The error of the floating point summation for the supremum
	
	\sa static inline void TwoSum(const real x, const real y, real &a, real &b)
	*/ 
	static inline void TwoSum(const interval &x, const interval &y, interval &a, 
					real &b_i, real &b_s) {
		real a_i,a_s;
		TwoSum(Inf(x),Inf(y),a_i,b_i);
		TwoSum(Sup(x),Sup(y),a_s,b_s);
		a = interval(a_i,a_s);
	}  
	


/*	static INLINE void TwoSum(const interval &x, const real &y, interval &a, 
					real &b_i, real &b_s) {
		real a_i,a_s;
		TwoSum(Inf(x),y,a_i,b_i);
		TwoSum(Sup(x),y,a_s,b_s);
		a = interval(a_i,a_s);
	}  */

	/*
	Computes the product of two intervals as well as the error
	occuring during multiplication. The multiplication is performed using the
	TwoProduct function for two real variables an utilizing the following
	multiplication table for intervals:
	
	
			bi,bs >= 0       bi < 0, bs >=0       bi,bs < 0
	----------------+------------------+------------------+----------------
	ai,as >= 0      I   ai*bi, as*bs   I   as*bi, as*bs   I   as*bi, ai*bs
	----------------+------------------+------------------+----------------
	ai < 0, as >= 0 I   ai*bs, as*bs   I      ....        I   as*bi, ai*bi
	----------------+------------------+------------------+----------------
	ai,as < 0       I   ai*bs, as*bi   I   ai*bs, ai*bi   I   as*bs, ai*bi
	----------------+------------------+------------------+----------------
	
	.... :  min(ai*bs, as*bi), max(ai*ai, as*as)
	
	\param x First factor
	\param y Second factor
	\param a Result of the floating point multiplication
	\param b_inf Error of the floating point multiplication of the infimum
	\param b_sup Error of the floating point multiplication of the supremum
	\sa static inline void TwoProduct(const real &x, const real &y, real &a, real &b) 
	*/
	static inline void TwoProduct(const interval &x, const interval &y, interval &a, real &b_inf, real &b_sup) {
		real a_inf, a_sup;
	
	
		if(Inf(x) >= 0  &&  Sup(x) >= 0) {
	
			if(Inf(y) >= 0  &&  Sup(y) >= 0) {
				TwoProduct(Inf(x), Inf(y), a_inf, b_inf);
				TwoProduct(Sup(x), Sup(y), a_sup, b_sup);
			} else if(Inf(y) < 0  &&  Sup(y) >= 0) {
				TwoProduct(Sup(x), Inf(y), a_inf, b_inf);
				TwoProduct(Sup(x), Sup(y), a_sup, b_sup);
			} else {
				TwoProduct(Sup(x), Inf(y), a_inf, b_inf);
				TwoProduct(Inf(x), Sup(y), a_sup, b_sup);
			}
	
		} else if(Inf(x) < 0  &&  Sup(x) >= 0) {
	
			if(Inf(y) >= 0  &&  Sup(y) >= 0) {
				TwoProduct(Inf(x), Sup(y), a_inf, b_inf);
				TwoProduct(Sup(x), Sup(y), a_sup, b_sup);
			} else if(Inf(y) < 0  &&  Sup(y) >= 0) {

				real ta1,ta2,tb1,tb2;

				TwoProduct(Inf(x), Sup(y), ta1, tb1);
				TwoProduct(Sup(x), Inf(y), ta2, tb2);
				if((ta1 < ta2) || (ta1 == ta2 && tb1 < tb2)) {
					a_inf = ta1;
					b_inf = tb1;
				} else {
					a_inf = ta2;
					b_inf = tb2;
				}

				TwoProduct(Inf(x), Inf(y), ta1, tb1);
				TwoProduct(Sup(x), Sup(y), ta2, tb2);
				if((ta1 > ta2) || (ta1 == ta2 && tb1 > tb2)) {
					a_sup = ta1;
					b_sup = tb1;
				} else {
					a_sup = ta2;
					b_sup = tb2;
				}
				
			} else {
				TwoProduct(Sup(x), Inf(y), a_inf, b_inf);
				TwoProduct(Inf(x), Inf(y), a_sup, b_sup);
			}
	
		} else {
	
			if(Inf(y) >= 0  &&  Sup(y) >= 0) {
				TwoProduct(Inf(x), Sup(y), a_inf, b_inf);
				TwoProduct(Sup(x), Inf(y), a_sup, b_sup);
			} else if(Inf(y) < 0  &&  Sup(y) >= 0) {
				TwoProduct(Inf(x), Sup(y), a_inf, b_inf);
				TwoProduct(Inf(x), Inf(y), a_sup, b_sup);
			} else {
				TwoProduct(Sup(x), Sup(y), a_inf, b_inf);
				TwoProduct(Inf(x), Inf(y), a_sup, b_sup);
			}
	
		}
	
		SetInf(a,a_inf); SetSup(a,a_sup);
	}

	
	/*
	Computes the product interval*real as well as the error
	occuring during multiplication. The multiplication is performed using the
	TwoProduct function for two real variables an utilizing the following
	multiplication table for intervals:
	
	
			bi,bs >= 0       bi < 0, bs >=0       bi,bs < 0
	----------------+------------------+------------------+----------------
	ai,as >= 0      I   ai*bi, as*bs   I   as*bi, as*bs   I   as*bi, ai*bs
	----------------+------------------+------------------+----------------
	ai < 0, as >= 0 I   ai*bs, as*bs   I      ....        I   as*bi, ai*bi
	----------------+------------------+------------------+----------------
	ai,as < 0       I   ai*bs, as*bi   I   ai*bs, ai*bi   I   as*bs, ai*bi
	----------------+------------------+------------------+----------------
	
	.... :  min(ai*bs, as*bi), max(ai*ai, as*as)
	
	\param x First factor
	\param y Second factor
	\param a Result of the floating point multiplication
	\param b_inf Error of the floating point multiplication of the infimum
	\param b_sup Error of the floating point multiplication of the supremum
	\sa static inline void TwoProduct(const real &x, const real &y, real &a, real &b) 
	*/
	static inline void TwoProduct(const interval &x, const real &y, interval &a, real &b_inf, real &b_sup) {
		real a_inf, a_sup;
		a = interval(0,0);
	
		if(Inf(x) >= 0  &&  Sup(x) >= 0) {
	
			if(y >= 0) {
				TwoProduct(Inf(x), y, a_inf, b_inf);
				TwoProduct(Sup(x), y, a_sup, b_sup);
			} else {
				TwoProduct(Sup(x), y, a_inf, b_inf);
				TwoProduct(Inf(x), y, a_sup, b_sup);
			}
	
		} else if(Inf(x) < 0  &&  Sup(x) >= 0) {
	
			if(y >= 0) {
				TwoProduct(Inf(x), y, a_inf, b_inf);
				TwoProduct(Sup(x), y, a_sup, b_sup);
			} else {
				TwoProduct(Sup(x), y, a_inf, b_inf);
				TwoProduct(Inf(x), y, a_sup, b_sup);
			}
	
		} else {
	
			if(y >= 0) {
				TwoProduct(Inf(x), y, a_inf, b_inf);
				TwoProduct(Sup(x), y, a_sup, b_sup);
			} else {
				TwoProduct(Sup(x), y, a_inf, b_inf);
				TwoProduct(Inf(x), y, a_sup, b_sup);
			}
	
		}
	
		SetInf(a,a_inf); SetSup(a,a_sup);
	}

	static inline void TwoProduct(const real &y, const interval &x, interval &a, real &b_inf, real &b_sup) {
		real a_inf, a_sup;
		a = interval(0,0);
	
		if(Inf(x) >= 0  &&  Sup(x) >= 0) {
	
			if(y >= 0) {
				TwoProduct(Inf(x), y, a_inf, b_inf);
				TwoProduct(Sup(x), y, a_sup, b_sup);
			} else {
				TwoProduct(Sup(x), y, a_inf, b_inf);
				TwoProduct(Inf(x), y, a_sup, b_sup);
			}
	
		} else if(Inf(x) < 0  &&  Sup(x) >= 0) {
	
			if(y >= 0) {
				TwoProduct(Inf(x), y, a_inf, b_inf);
				TwoProduct(Sup(x), y, a_sup, b_sup);
			} else {
				TwoProduct(Sup(x), y, a_inf, b_inf);
				TwoProduct(Inf(x), y, a_sup, b_sup);
			}
	
		} else {
	
			if(y >= 0) {
				TwoProduct(Inf(x), y, a_inf, b_inf);
				TwoProduct(Sup(x), y, a_sup, b_sup);
			} else {
				TwoProduct(Sup(x), y, a_inf, b_inf);
				TwoProduct(Inf(x), y, a_sup, b_sup);
			}
	
		}
	
		SetInf(a,a_inf); SetSup(a,a_sup);
	}

	//Sorts inf and sup of two interval vectors into two real vector for the computation of the infimum 
	//and two real vectors for the computation of the supremum. Rounding mode must be set to upwards!
	static inline void sort(const ivector &x, const ivector &y, rvector& x_inf, rvector& y_inf, rvector &x_sup, rvector &y_sup, int n, int lb1, int lb2) {
		
		for(int i=0 ; i<n ; i++) {

			if(Inf(x[i+lb1]) >= 0  &&  Sup(x[i+lb1]) >= 0) {
		
				if(Inf(y[i+lb2]) >= 0  &&  Sup(y[i+lb2]) >= 0) {
					x_inf[i+1] = Inf(x[i+lb1]);
					y_inf[i+1] = Inf(y[i+lb2]);
					x_sup[i+1] = Sup(x[i+lb1]);
					y_sup[i+1] = Sup(y[i+lb2]);
				} else if(Inf(y[i+lb2]) < 0  &&  Sup(y[i+lb2]) >= 0) {
					x_inf[i+1] = Sup(x[i+lb1]);
					y_inf[i+1] = Inf(y[i+lb2]);
					x_sup[i+1] = Sup(x[i+lb1]);
					y_sup[i+1] = Sup(y[i+lb2]);
				} else {
					x_inf[i+1] = Sup(x[i+lb1]);
					y_inf[i+1] = Inf(y[i+lb2]);
					x_sup[i+1] = Inf(x[i+lb1]);
					y_sup[i+1] = Sup(y[i+lb2]);
				}
		
			} else if(Inf(x[i+lb1]) < 0  &&  Sup(x[i+lb1]) >= 0) {
		
				if(Inf(y[i+lb2]) >= 0  &&  Sup(y[i+lb2]) >= 0) {
					x_inf[i+1] = Inf(x[i+lb1]);
					y_inf[i+1] = Sup(y[i+lb2]);
					x_sup[i+1] = Sup(x[i+lb1]);
					y_sup[i+1] = Sup(y[i+lb2]);
				} else if(Inf(y[i+lb2]) < 0  &&  Sup(y[i+lb2]) >= 0) {

#ifndef _CXSC_DOTK_ROUND_SOFT
					if((-Inf(x[i+lb1]))*Sup(y[i+lb2]) >= Sup(x[i+lb1])*(-Inf(y[i+lb2]))) {
#else
					if( mulu(-Inf(x[i+lb1]),Sup(y[i+lb2])) >= mulu(Sup(x[i+lb1]),-Inf(y[i+lb2])) ) {
#endif
						x_inf[i+1] = Inf(x[i+lb1]);
						y_inf[i+1] = Sup(y[i+lb2]);
					} else {
						x_inf[i+1] = Sup(x[i+lb1]);
						y_inf[i+1] = Inf(y[i+lb2]);
					}

#ifndef _CXSC_DOTK_ROUND_SOFT
					if(Inf(x[i+lb1])*Inf(y[i+lb2]) >= Sup(x[i+lb1])*Sup(y[i+lb2])) {
#else
					if(mulu(Inf(x[i+lb1]),Inf(y[i+lb2])) >= mulu(Sup(x[i+lb1]),Sup(y[i+lb2]))) {
#endif
						x_sup[i+1] = Inf(x[i+lb1]);
						y_sup[i+1] = Inf(y[i+lb2]);
					} else {
						x_sup[i+1] = Sup(x[i+lb1]);
						y_sup[i+1] = Sup(y[i+lb2]);
					}
					
				} else {
					x_inf[i+1] = Sup(x[i+lb1]);
					y_inf[i+1] = Inf(y[i+lb2]);
					x_sup[i+1] = Inf(x[i+lb1]);
					y_sup[i+1] = Inf(y[i+lb2]);
				}
		
			} else {
		
				if(Inf(y[i+lb2]) >= 0  &&  Sup(y[i+lb2]) >= 0) {
					x_inf[i+1] = Inf(x[i+lb1]);
					y_inf[i+1] = Sup(y[i+lb2]);
					x_sup[i+1] = Sup(x[i+lb1]);
					y_sup[i+1] = Inf(y[i+lb2]);
				} else if(Inf(y[i+lb2]) < 0  &&  Sup(y[i+lb2]) >= 0) {
					x_inf[i+1] = Inf(x[i+lb1]);
					y_inf[i+1] = Sup(y[i+lb2]);
					x_sup[i+1] = Inf(x[i+lb1]);
					y_sup[i+1] = Inf(y[i+lb2]);
				} else {
					x_inf[i+1] = Sup(x[i+lb1]);
					y_inf[i+1] = Sup(y[i+lb2]);
					x_sup[i+1] = Inf(x[i+lb1]);
					y_sup[i+1] = Inf(y[i+lb2]);
				}
		
			}
		
		}
	}

	//Sorts inf and sup of an interval vector and a real vector into two real vector for the computation of the infimum 
	//and two real vectors for the computation of the supremum. Rounding mode must be set to upwards!
	static inline void sort(const ivector &x, const rvector &y, rvector& x_inf, rvector& y_inf, rvector &x_sup, rvector &y_sup, int n, int lb1, int lb2) {
		
		for(int i=0 ; i<n ; i++) {

			if(Inf(x[i+lb1]) >= 0  &&  Sup(x[i+lb1]) >= 0) {
		
				if(y[i+lb2] >= 0) {
					x_inf[i+1] = Inf(x[i+lb1]);
					y_inf[i+1] = y[i+lb2];
					x_sup[i+1] = Sup(x[i+lb1]);
					y_sup[i+1] = y[i+lb2];
				} else {
					x_inf[i+1] = Sup(x[i+lb1]);
					y_inf[i+1] = y[i+lb2];
					x_sup[i+1] = Inf(x[i+lb1]);
					y_sup[i+1] = y[i+lb2];
				}
		
			} else if(Inf(x[i+lb1]) < 0  &&  Sup(x[i+lb1]) >= 0) {
		
				if(y[i+lb2] >= 0) {
					x_inf[i+1] = Inf(x[i+lb1]);
					y_inf[i+1] = y[i+lb2];
					x_sup[i+1] = Sup(x[i+lb1]);
					y_sup[i+1] = y[i+lb2];
				} else {
					x_inf[i+1] = Sup(x[i+lb1]);
					y_inf[i+1] = y[i+lb2];
					x_sup[i+1] = Inf(x[i+lb1]);
					y_sup[i+1] = y[i+lb2];
				}
		
			} else {
		
				if(y[i+lb2] >= 0) {
					x_inf[i+1] = Inf(x[i+lb1]);
					y_inf[i+1] = y[i+lb2];
					x_sup[i+1] = Sup(x[i+lb1]);
					y_sup[i+1] = y[i+lb2];
				} else {
					x_inf[i+1] = Sup(x[i+lb1]);
					y_inf[i+1] = y[i+lb2];
					x_sup[i+1] = Inf(x[i+lb1]);
					y_sup[i+1] = y[i+lb2];
				}
		
			}
		
		}
	}

	//Sorts inf and sup of an interval vector and a real vector into two real vector for the computation of the infimum 
	//and two real vectors for the computation of the supremum. Rounding mode must be set to upwards!
	static inline void sort(const rvector &y, const ivector &x, rvector& x_inf, rvector& y_inf, rvector &x_sup, rvector &y_sup, int n, int lb2, int lb1) {
		for(int i=0 ; i<n ; i++) {

			if(Inf(x[i+lb1]) >= 0  &&  Sup(x[i+lb1]) >= 0) {
		
				if(y[i+lb2] >= 0) {
					x_inf[i+1] = Inf(x[i+lb1]);
					y_inf[i+1] = y[i+lb2];
					x_sup[i+1] = Sup(x[i+lb1]);
					y_sup[i+1] = y[i+lb2];
				} else {
					x_inf[i+1] = Sup(x[i+lb1]);
					y_inf[i+1] = y[i+lb2];
					x_sup[i+1] = Inf(x[i+lb1]);
					y_sup[i+1] = y[i+lb2];
				}
		
			} else if(Inf(x[i+lb1]) < 0  &&  Sup(x[i+lb1]) >= 0) {
		
				if(y[i+lb2] >= 0) {
					x_inf[i+1] = Inf(x[i+lb1]);
					y_inf[i+1] = y[i+lb2];
					x_sup[i+1] = Sup(x[i+lb1]);
					y_sup[i+1] = y[i+lb2];
				} else {
					x_inf[i+1] = Sup(x[i+lb1]);
					y_inf[i+1] = y[i+lb2];
					x_sup[i+1] = Inf(x[i+lb1]);
					y_sup[i+1] = y[i+lb2];
				}
		
			} else {
		
				if(y[i+lb2] >= 0) {
					x_inf[i+1] = Inf(x[i+lb1]);
					y_inf[i+1] = y[i+lb2];
					x_sup[i+1] = Sup(x[i+lb1]);
					y_sup[i+1] = y[i+lb2];
				} else {
					x_inf[i+1] = Sup(x[i+lb1]);
					y_inf[i+1] = y[i+lb2];
					x_sup[i+1] = Inf(x[i+lb1]);
					y_sup[i+1] = y[i+lb2];
				}
		
			}
		
		}
	}


	/*
	Computes the interval dot product in floating point using midpoint radius representation
	*/
/*	static inline interval floatidot(const ivector& x, const ivector& y, int n, int lb1, int lb2) {
		rvector xmid(n),xrad(n),ymid(n),yrad(n);
		real crad = 0.0;
		real res_inf=0.0,res_sup=0.0;

		setround(1);

		for(int i=0 ; i<n ; i++) {
			xmid[i+1] = Inf(x[lb1+i]) + 0.5*(Sup(x[lb1+i])-Inf(x[lb1+i]));
			ymid[i+1] = Inf(y[lb2+i]) + 0.5*(Sup(y[lb2+i])-Inf(y[lb2+i]));				

		}

		for(int i=0 ; i<n ; i++) {
			xrad[i+1] = xmid[i+1] - Inf(x[i+lb1]);
			yrad[i+1] = ymid[i+1] - Inf(y[i+lb1]);
		}

		for(int i=0 ; i<n ; i++) {
			crad += xrad[i+1] * (abs(ymid[i+1]) + yrad[i+1]) + abs(xmid[i+1]) * yrad[i+1];
			res_sup += xmid[i+1] * ymid[i+1];
		}

		res_sup += crad;

		setround(-1);
		for(int i=0 ; i<n ; i++) {
			res_inf += xmid[i+1] * ymid[i+1];
		}

		res_inf -= crad;
		setround(0);
		return interval(res_inf,res_sup);
	}*/

	/*
	Computes the interval dot product in floating point using midpoint radius representation
	*/
/*	static inline interval floatidot(const ivector& x, const rvector& y, int n, int lb1, int lb2) {
		rvector xmid(n),xrad(n);
		real crad = 0.0;
		real res_inf=0.0,res_sup=0.0;

		setround(1);

		for(int i=0 ; i<n ; i++) {
			xmid[i+1] = Inf(x[lb1+i]) + 0.5*(Sup(x[lb1+i])-Inf(x[lb1+i]));

		}

		for(int i=0 ; i<n ; i++) {
			xrad[i+1] = xmid[i+1] - Inf(x[i+lb1]);
		}

		for(int i=0 ; i<n ; i++) {
			crad += xrad[i+1] * abs(y[i+lb2]);
			res_sup += xmid[i+1] * y[i+lb2];
		}

		res_sup += crad;

		setround(-1);
		for(int i=0 ; i<n ; i++) {
			res_inf += xmid[i+1] * y[i+lb2];
		}

		res_inf -= crad;
		setround(0);
		return interval(res_inf,res_sup);
	}*/

	/*
	Computes the interval dot product in floating point using midpoint radius representation
	*/
/*	static inline interval floatidot(const rvector& x, const ivector& y, int n, int lb1, int lb2) {		
		return floatidot(y,x,n,lb2,lb1);
	}*/

	/*
	Computes the dot product of two interval vectors in k-fold double precision.
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
	static inline void addDot(idotprecision &val, const S &x, const T &y) {
		int n = VecLen(x);
		int lb1 = Lb(x);
		int lb2 = Lb(y);
		real err_inf = 0.0;
		real err_sup = 0.0;
		int k = val.get_k();

#ifndef _CXSC_DOTK_ROUND_SOFT
		int rnd;
		//Check rounding mode
		if((rnd=getround()) != 0) {
			setround(0);
		}
#endif
	
		if(k == 0) { //Use accumulator
			
			for(int i=1 ; i<=n ; i++)
				accumulate(val, x[i+lb1-1], y[i+lb2-1]);

		} else if(k == 1) { //Use floating point
			
#ifndef _CXSC_DOTK_ROUND_SOFT
			rvector x_inf(n),y_inf(n),x_sup(n),y_sup(n);
			interval r(0.0,0.0);

			setround(1);
                        sort(x,y,x_inf,y_inf,x_sup,y_sup,n,lb1,lb2);
			for(int i=1 ; i<=n ; i++)
				Sup(r) += x_sup[i] * y_sup[i];

			setround(-1);
			for(int i=1 ; i<=n ; i++)
				Inf(r) += x_inf[i] * y_inf[i];

			val += r;
#else
			interval r(0.0,0.0);
			for(int i=1 ; i<=n ; i++)
				r += x[i+lb1-1] * y[i+lb2-1];
			val += r;
#endif
			/*val += floatidot(x,y,n,lb1,lb2);*/

		} else if(k == 2) { //Use optimized DotK

			interval p,h;
			real s_inf, s_sup, r_inf, r_sup;
			real p_inf, p_sup, q_inf, q_sup;
			real t;
	
			TwoProduct(x[lb1], y[lb2], p, s_inf, s_sup);
	
			err_inf += abs(s_inf);
			err_sup += abs(s_sup);
	
			for(int i=2 ; i<=n ; i++) {
				TwoProduct(x[lb1+i-1], y[lb2+i-1], h, r_inf, r_sup);
				TwoSum(p, h, p, q_inf, q_sup);
	
				t = q_inf + r_inf;
				s_inf += t;
				err_inf += abs(t);
	
				t = q_sup + r_sup;
				s_sup += t;
				err_sup += abs(t);
			}
		
			dotprecision& val_inf = Inf(val);
			dotprecision& val_sup = Sup(val);

			val_inf += Inf(p);
			val_inf += s_inf;

			val_sup += Sup(p);
			val_sup += s_sup;

			real alpha, delta, error;
			delta = (n*Epsilon) / (1.0-2*n*Epsilon);

			//Error infimum	
			alpha = (Epsilon*abs(Inf(p))) + (delta*err_inf+3*MinReal/Epsilon);
			error = alpha / (1.0 - 2*Epsilon);
			val_inf -= error;

			//Error supremum
			alpha = (Epsilon*abs(Sup(p))) + (delta*err_sup+3*MinReal/Epsilon);
			error = alpha / (1.0 - 2*Epsilon);
			val_sup += error;

			//UncheckedSetInf(val,val_inf);
			//UncheckedSetSup(val,val_sup);                          
		
		} else { //Use DotK
		
			interval r(0,0), h;
			real* t_inf = new real[2*n];
			real* t_sup = new real[2*n];
	
			for(int i=1 ; i<=n ; i++) {
				TwoProduct(x[lb1+i-1], y[lb2+i-1], h, t_inf[i-1], t_sup[i-1]);
				TwoSum(r, h, r, t_inf[n+i-2], t_sup[n+i-2]);
			}              
	
			t_inf[2*n-1] = Inf(r);
			t_sup[2*n-1] = Sup(r);

			dotprecision& val_inf = Inf(val);
			dotprecision& val_sup = Sup(val);

			SumK(t_inf, 2*n, k-1, err_inf, val_inf);
			SumK(t_sup, 2*n, k-1, err_sup, val_sup);

			val_inf -= err_inf;
			val_sup += err_sup;
		
			//UncheckedSetInf(val, val_inf);
			//UncheckedSetSup(val, val_sup);

			delete [] t_inf;
			delete [] t_sup;	

		}

#ifndef _CXSC_DOTK_ROUND_SOFT
		//Reset rounding mode
		setround(rnd);
#endif

	}


}

#endif // _CXSC_IDOTK_INL_INCLUDED
