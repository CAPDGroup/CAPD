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

/* CVS $Id: cxsc_blas.inl,v 1.20 2014/01/30 17:23:44 cxsc Exp $ */

/*
**  FastPLSS: A library of (parallel) verified linear (interval) system 
**  solvers using C-XSC (V 0.2)
*/

#include <fenv.h>
#include "cxsc_blas.hpp"
 
namespace cxsc {

#ifndef _CXSC_BLAS_SETROUND
#define _CXSC_BLAS_SETROUND
	/*
	   Sets the rounding mode according to the parameter r:
	   r=-1: Round downwards
           r=0 : Round to nearest
           r=1 : Round upwards
	   r=2 : Round toward zero
	*/
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

	/*
	   Sets the rounding mode according to the parameter r:
	   r=-1: Round downwards
           r=0 : Round to nearest
           r=1 : Round upwards
	   r=2 : Round toward zero
	*/
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
#endif

#ifdef _CXSC_BLAS_RVECTOR
#ifndef _CXSC_BLAS_RVECTOR_INC
#define _CXSC_BLAS_RVECTOR_INC
inline void blasdot(const rvector& x, const rvector& y, real& res) {
   int n = VecLen(x);
   int inc=1;
   
   double* xd = (double*) x.to_blas_array();
   double* yd = (double*) y.to_blas_array();
   
   res = cblas_ddot(n, xd, inc, yd, inc);
}

inline void blasdot(const rvector& x, const rvector& y, interval& res) {
   int n = VecLen(x);
   int inc=1;
   int rnd = getround();
   
   double* xd = (double*) x.to_blas_array();
   double* yd = (double*) y.to_blas_array();
   
   setround(-1);
   SetInf(res, cblas_ddot(n, xd, inc, yd, inc));
   setround(1);
   SetSup(res, cblas_ddot(n, xd, inc, yd, inc));
   setround(rnd);
}
#endif
#endif

#if defined(_CXSC_BLAS_CVECTOR) && defined(_CXSC_BLAS_IVECTOR)
#if !defined(_CXSC_BLAS_CVECTOR_INC) || !defined(_CXSC_BLAS_IVECTOR_INC)
inline void blasdot(const cvector& x, const ivector& y, cinterval& res) {
   const int n = VecLen(x);
   const int inc=1;
   const int lb1 = Lb(x);
   const int lb2 = Lb(y);
   int rnd = getround();

   rvector x_inf(n),x_sup(n),y_inf(n),y_sup(n);

   double* dxi = x_inf.to_blas_array();
   double* dxs = x_sup.to_blas_array();
   double* dyi = y_inf.to_blas_array();
   double* dys = y_sup.to_blas_array();
   double res_inf,res_sup;

   setround(1);

   bsort(Re(x),y,x_inf,y_inf,x_sup,y_sup,n,lb1,lb2);

   res_sup = cblas_ddot(n, dxs, inc, dys, inc);

   //setround(-1);

   res_inf = -cblas_ddot(n, dxi, inc, dyi, inc);

   SetRe(res, interval(res_inf,res_sup));

   bsort(Im(x),y,x_inf,y_inf,x_sup,y_sup,n,lb1,lb2);

   res_sup = cblas_ddot(n, dxs, inc, dys, inc);

   //setround(-1);

   res_inf = -cblas_ddot(n, dxi, inc, dyi, inc);

   SetIm(res, interval(res_inf,res_sup));

   setround(rnd);
}

inline void blasdot(const ivector& x, const cvector& y, cinterval& res) {
   const int n = VecLen(x);
   const int inc=1;
   const int lb1 = Lb(x);
   const int lb2 = Lb(y);
   int rnd = getround();

   rvector x_inf(n),x_sup(n),y_inf(n),y_sup(n);

   double* dxi = x_inf.to_blas_array();
   double* dxs = x_sup.to_blas_array();
   double* dyi = y_inf.to_blas_array();
   double* dys = y_sup.to_blas_array();
   double res_inf,res_sup;

   setround(1);

   bsort(x,Re(y),x_inf,y_inf,x_sup,y_sup,n,lb1,lb2);

   res_sup = cblas_ddot(n, dxs, inc, dys, inc);

   //setround(-1);

   res_inf = -cblas_ddot(n, dxi, inc, dyi, inc);

   SetRe(res, interval(res_inf,res_sup));

   bsort(x,Im(y),x_inf,y_inf,x_sup,y_sup,n,lb1,lb2);

   res_sup = cblas_ddot(n, dxs, inc, dys, inc);

   //setround(-1);

   res_inf = -cblas_ddot(n, dxi, inc, dyi, inc);

   SetIm(res, interval(res_inf,res_sup));

   setround(rnd);
}
#endif
#endif

#ifdef _CXSC_BLAS_CVECTOR
#ifndef _CXSC_BLAS_CVECTOR_INC
#define _CXSC_BLAS_CVECTOR_INC
inline void blasdot(const rvector& x, const cvector& y, complex& res) {
   int n = VecLen(x);
   int inc=1, inc2=2;
   double res_re, res_im;
   
   double* xd = (double*) x.to_blas_array();
   double* yd = (double*) y.to_blas_array();
   
   res_re = cblas_ddot(n, xd, inc, yd, inc2);
   res_im = cblas_ddot(n, xd, inc, yd+1, inc2);

   res = complex(res_re,res_im);
}

inline void blasdot(const cvector& x, const rvector& y, complex& res) {
   int n = VecLen(x);
   int inc=1, inc2=2;
   double res_re, res_im;
   
   double* xd = (double*) x.to_blas_array();
   double* yd = (double*) y.to_blas_array();
   
   res_re = cblas_ddot(n, xd, inc2, yd, inc);
   res_im = cblas_ddot(n, xd+1, inc2, yd, inc);

   res = complex(res_re,res_im);
}

inline void blasdot(const cvector& x, const cvector& y, complex& res) {
   int n = VecLen(x);
   int inc=1;
   
   double* xd = (double*) x.to_blas_array();
   double* yd = (double*) y.to_blas_array();
   
   cblas_zdotu_sub(n, xd, inc, yd, inc, (void*)&res);
}

inline void blasdot(const rvector& x, const cvector& y, cinterval& res) {
   interval re,im;

   blasdot(x,Re(y),re);
   blasdot(x,Im(y),im);

   SetRe(res,re);
   SetIm(res,im);
}

inline void blasdot(const cvector& x, const rvector& y, cinterval& res) {
   interval re,im;

   blasdot(Re(x),y,re);
   blasdot(Im(x),y,im);

   SetRe(res,re);
   SetIm(res,im);
}


inline void blasdot(const cvector& x, const cvector& y, cinterval& res) {
   int rnd = getround();
   int n = VecLen(x);
   int inc=2;
   double res_re, res_im;
   
   double* xd = (double*) x.to_blas_array();
   double* yd = (double*) y.to_blas_array();
   
   setround(1);
   res_re = - cblas_ddot(n, xd+1, inc, yd+1, inc);

   setround(-1);
   res_re += cblas_ddot(n, xd, inc, yd, inc);
   res_im = cblas_ddot(n, xd, inc, yd+1, inc);
   res_im += cblas_ddot(n, xd+1, inc, yd, inc);
   UncheckedSetInf(res,complex(res_re,res_im));

   res_re = - cblas_ddot(n, xd+1, inc, yd+1, inc);
   setround(1);
   res_re += cblas_ddot(n, xd, inc, yd, inc);
   res_im = cblas_ddot(n, xd, inc, yd+1, inc);
   res_im += cblas_ddot(n, xd+1, inc, yd, inc);
   UncheckedSetSup(res,complex(res_re,res_im));
   
   setround(rnd);
}
#endif
#endif

#ifdef _CXSC_BLAS_IVECTOR
#ifndef _CXSC_BLAS_IVECTOR_INC
#define _CXSC_BLAS_IVECTOR_INC
inline void bsort(const ivector &x, const ivector &y, rvector& x_inf, rvector& y_inf, rvector &x_sup, rvector &y_sup, int n, int lb1, int lb2) {
	
	for(int i=0 ; i<n ; i++) {

		if(Inf(x[i+lb1]) >= 0  &&  Sup(x[i+lb1]) >= 0) {
	
			if(Inf(y[i+lb2]) >= 0  &&  Sup(y[i+lb2]) >= 0) {
				x_inf[i+1] = -Inf(x[i+lb1]);
				y_inf[i+1] = Inf(y[i+lb2]);
				x_sup[i+1] = Sup(x[i+lb1]);
				y_sup[i+1] = Sup(y[i+lb2]);
			} else if(Inf(y[i+lb2]) < 0  &&  Sup(y[i+lb2]) >= 0) {
				x_inf[i+1] = -Sup(x[i+lb1]);
				y_inf[i+1] = Inf(y[i+lb2]);
				x_sup[i+1] = Sup(x[i+lb1]);
				y_sup[i+1] = Sup(y[i+lb2]);
			} else {
				x_inf[i+1] = -Sup(x[i+lb1]);
				y_inf[i+1] = Inf(y[i+lb2]);
				x_sup[i+1] = Inf(x[i+lb1]);
				y_sup[i+1] = Sup(y[i+lb2]);
			}
	
		} else if(Inf(x[i+lb1]) < 0  &&  Sup(x[i+lb1]) >= 0) {
	
			if(Inf(y[i+lb2]) >= 0  &&  Sup(y[i+lb2]) >= 0) {
				x_inf[i+1] = -Inf(x[i+lb1]);
				y_inf[i+1] = Sup(y[i+lb2]);
				x_sup[i+1] = Sup(x[i+lb1]);
				y_sup[i+1] = Sup(y[i+lb2]);
			} else if(Inf(y[i+lb2]) < 0  &&  Sup(y[i+lb2]) >= 0) {

				if((-Inf(x[i+lb1]))*Sup(y[i+lb2]) >= Sup(x[i+lb1])*(-Inf(y[i+lb2]))) {
					x_inf[i+1] = -Inf(x[i+lb1]);
					y_inf[i+1] = Sup(y[i+lb2]);
				} else {
					x_inf[i+1] = -Sup(x[i+lb1]);
					y_inf[i+1] = Inf(y[i+lb2]);
				}

				if(Inf(x[i+lb1])*Inf(y[i+lb2]) >= Sup(x[i+lb1])*Sup(y[i+lb2])) {
					x_sup[i+1] = Inf(x[i+lb1]);
					y_sup[i+1] = Inf(y[i+lb2]);
				} else {
					x_sup[i+1] = Sup(x[i+lb1]);
					y_sup[i+1] = Sup(y[i+lb2]);
				}
				
			} else {
				x_inf[i+1] = Sup(x[i+lb1]);
				y_inf[i+1] = -Inf(y[i+lb2]);
				x_sup[i+1] = Inf(x[i+lb1]);
				y_sup[i+1] = Inf(y[i+lb2]);
			}
	
		} else {
	
			if(Inf(y[i+lb2]) >= 0  &&  Sup(y[i+lb2]) >= 0) {
				x_inf[i+1] = -Inf(x[i+lb1]);
				y_inf[i+1] = Sup(y[i+lb2]);
				x_sup[i+1] = Sup(x[i+lb1]);
				y_sup[i+1] = Inf(y[i+lb2]);
			} else if(Inf(y[i+lb2]) < 0  &&  Sup(y[i+lb2]) >= 0) {
				x_inf[i+1] = -Inf(x[i+lb1]);
				y_inf[i+1] = Sup(y[i+lb2]);
				x_sup[i+1] = Inf(x[i+lb1]);
				y_sup[i+1] = Inf(y[i+lb2]);
			} else {
				x_inf[i+1] = -Sup(x[i+lb1]);
				y_inf[i+1] = Sup(y[i+lb2]);
				x_sup[i+1] = Inf(x[i+lb1]);
				y_sup[i+1] = Inf(y[i+lb2]);
			}
	
		}
	
	}
}

//Sorts inf and sup of an interval vector and a real vector into two real vector for the computation of the infimum 
//and two real vectors for the computation of the supremum. Rounding mode must be set to upwards!
inline void bsort(const ivector &x, const rvector &y, rvector& x_inf, rvector& y_inf, rvector &x_sup, rvector &y_sup, int n, int lb1, int lb2) {
	
	for(int i=0 ; i<n ; i++) {

		if(Inf(x[i+lb1]) >= 0  &&  Sup(x[i+lb1]) >= 0) {
	
			if(y[i+lb2] >= 0) {
				x_inf[i+1] = -Inf(x[i+lb1]);
				y_inf[i+1] = y[i+lb2];
				x_sup[i+1] = Sup(x[i+lb1]);
				y_sup[i+1] = y[i+lb2];
			} else {
				x_inf[i+1] = -Sup(x[i+lb1]);
				y_inf[i+1] = y[i+lb2];
				x_sup[i+1] = Inf(x[i+lb1]);
				y_sup[i+1] = y[i+lb2];
			}
	
		} else if(Inf(x[i+lb1]) < 0  &&  Sup(x[i+lb1]) >= 0) {
	
			if(y[i+lb2] >= 0) {
				x_inf[i+1] = -Inf(x[i+lb1]);
				y_inf[i+1] = y[i+lb2];
				x_sup[i+1] = Sup(x[i+lb1]);
				y_sup[i+1] = y[i+lb2];
			} else {
				x_inf[i+1] = -Sup(x[i+lb1]);
				y_inf[i+1] = y[i+lb2];
				x_sup[i+1] = Inf(x[i+lb1]);
				y_sup[i+1] = y[i+lb2];
			}
	
		} else {
	
			if(y[i+lb2] >= 0) {
				x_inf[i+1] = -Inf(x[i+lb1]);
				y_inf[i+1] = y[i+lb2];
				x_sup[i+1] = Sup(x[i+lb1]);
				y_sup[i+1] = y[i+lb2];
			} else {
				x_inf[i+1] = -Sup(x[i+lb1]);
				y_inf[i+1] = y[i+lb2];
				x_sup[i+1] = Inf(x[i+lb1]);
				y_sup[i+1] = y[i+lb2];
			}
	
		}
	
	}
}

//Sorts inf and sup of an interval vector and a real vector into two real vector for the computation of the infimum 
//and two real vectors for the computation of the supremum. Rounding mode must be set to upwards!
inline void bsort(const rvector &y, const ivector &x, rvector& x_inf, rvector& y_inf, rvector &x_sup, rvector &y_sup, int n, int lb2, int lb1) {
	for(int i=0 ; i<n ; i++) {

		if(Inf(x[i+lb1]) >= 0  &&  Sup(x[i+lb1]) >= 0) {
	
			if(y[i+lb2] >= 0) {
				x_inf[i+1] = -Inf(x[i+lb1]);
				y_inf[i+1] = y[i+lb2];
				x_sup[i+1] = Sup(x[i+lb1]);
				y_sup[i+1] = y[i+lb2];
			} else {
				x_inf[i+1] = -Sup(x[i+lb1]);
				y_inf[i+1] = y[i+lb2];
				x_sup[i+1] = Inf(x[i+lb1]);
				y_sup[i+1] = y[i+lb2];
			}
	
		} else if(Inf(x[i+lb1]) < 0  &&  Sup(x[i+lb1]) >= 0) {
	
			if(y[i+lb2] >= 0) {
				x_inf[i+1] = -Inf(x[i+lb1]);
				y_inf[i+1] = y[i+lb2];
				x_sup[i+1] = Sup(x[i+lb1]);
				y_sup[i+1] = y[i+lb2];
			} else {
				x_inf[i+1] = -Sup(x[i+lb1]);
				y_inf[i+1] = y[i+lb2];
				x_sup[i+1] = Inf(x[i+lb1]);
				y_sup[i+1] = y[i+lb2];
			}
	
		} else {
	
			if(y[i+lb2] >= 0) {
				x_inf[i+1] = -Inf(x[i+lb1]);
				y_inf[i+1] = y[i+lb2];
				x_sup[i+1] = Sup(x[i+lb1]);
				y_sup[i+1] = y[i+lb2];
			} else {
				x_inf[i+1] = -Sup(x[i+lb1]);
				y_inf[i+1] = y[i+lb2];
				x_sup[i+1] = Inf(x[i+lb1]);
				y_sup[i+1] = y[i+lb2];
			}
	
		}
	
	}
}

inline void blasdot(const ivector& x, const ivector& y, interval& res) {
   const int n = VecLen(x);
   const int inc=1;
   const int lb1 = Lb(x);
   const int lb2 = Lb(y);
   int rnd = getround();

   rvector x_inf(n),x_sup(n),y_inf(n),y_sup(n);

   double* dxi = x_inf.to_blas_array();
   double* dxs = x_sup.to_blas_array();
   double* dyi = y_inf.to_blas_array();
   double* dys = y_sup.to_blas_array();
   double res_inf,res_sup;

   setround(1);

   bsort(x,y,x_inf,y_inf,x_sup,y_sup,n,lb1,lb2);

   res_sup = cblas_ddot(n, dxs, inc, dys, inc);

   //setround(-1);

   res_inf = -cblas_ddot(n, dxi, inc, dyi, inc);

   setround(rnd);

   res = interval(res_inf,res_sup);
}

inline void blasdot(const ivector& x, const rvector& y, interval& res) {
   const int n = VecLen(x);
   const int inc=1;
   const int lb1 = Lb(x);
   const int lb2 = Lb(y);
   int rnd = getround();

   rvector x_inf(n),x_sup(n),y_inf(n),y_sup(n);

   double* dxi = x_inf.to_blas_array();
   double* dxs = x_sup.to_blas_array();
   double* dyi = y_inf.to_blas_array();
   double* dys = y_sup.to_blas_array();
   double res_inf,res_sup;

   setround(1);

   bsort(x,y,x_inf,y_inf,x_sup,y_sup,n,lb1,lb2);

   res_sup = cblas_ddot(n, dxs, inc, dys, inc);

   //setround(-1);

   res_inf = -cblas_ddot(n, dxi, inc, dyi, inc);

   setround(rnd);

   res = interval(res_inf,res_sup);
}

inline void blasdot(const rvector& x, const ivector& y, interval& res) {
   const int n = VecLen(x);
   const int inc=1;
   const int lb1 = Lb(x);
   const int lb2 = Lb(y);
   int rnd = getround();

   rvector x_inf(n),x_sup(n),y_inf(n),y_sup(n);

   double* dxi = x_inf.to_blas_array();
   double* dxs = x_sup.to_blas_array();
   double* dyi = y_inf.to_blas_array();
   double* dys = y_sup.to_blas_array();
   double res_inf,res_sup;

   setround(1);

   bsort(x,y,x_inf,y_inf,x_sup,y_sup,n,lb1,lb2);

   res_sup = cblas_ddot(n, dxs, inc, dys, inc);

   //setround(-1);

   res_inf = -cblas_ddot(n, dxi, inc, dyi, inc);

   setround(rnd);

   res = interval(res_inf,res_sup);
}
#endif
#endif


#ifdef _CXSC_BLAS_CIVECTOR
#ifndef _CXSC_BLAS_CIVECTOR_INC
#define _CXSC_BLAS_CIVECTOR_INC
inline void blasdot(const rvector& x, const civector& y, cinterval& res) {
   const int n = VecLen(x);
   const int inc=1;
   const int lb1 = Lb(x);
   const int lb2 = Lb(y);
   int rnd = getround();

   rvector x_inf(n),x_sup(n),y_inf(n),y_sup(n);

   double* dxi = x_inf.to_blas_array();
   double* dxs = x_sup.to_blas_array();
   double* dyi = y_inf.to_blas_array();
   double* dys = y_sup.to_blas_array();
   double res_inf,res_sup;

   setround(1);

   bsort(x,Re(y),x_inf,y_inf,x_sup,y_sup,n,lb1,lb2);

   res_sup = cblas_ddot(n, dxs, inc, dys, inc);

   //setround(-1);

   res_inf = -cblas_ddot(n, dxi, inc, dyi, inc);

   SetRe(res, interval(res_inf,res_sup));

   bsort(x,Im(y),x_inf,y_inf,x_sup,y_sup,n,lb1,lb2);

   res_sup = cblas_ddot(n, dxs, inc, dys, inc);

   //setround(-1);

   res_inf = -cblas_ddot(n, dxi, inc, dyi, inc);

   SetIm(res, interval(res_inf,res_sup));

   setround(rnd);
}


inline void blasdot(const civector& x, const rvector& y, cinterval& res) {
   const int n = VecLen(x);
   const int inc=1;
   const int lb1 = Lb(x);
   const int lb2 = Lb(y);
   int rnd = getround();

   rvector x_inf(n),x_sup(n),y_inf(n),y_sup(n);

   double* dxi = x_inf.to_blas_array();
   double* dxs = x_sup.to_blas_array();
   double* dyi = y_inf.to_blas_array();
   double* dys = y_sup.to_blas_array();
   double res_inf,res_sup;

   setround(1);

   bsort(Re(x),y,x_inf,y_inf,x_sup,y_sup,n,lb1,lb2);

   res_sup = cblas_ddot(n, dxs, inc, dys, inc);

   //setround(-1);

   res_inf = -cblas_ddot(n, dxi, inc, dyi, inc);

   SetRe(res, interval(res_inf,res_sup));

   bsort(Im(x),y,x_inf,y_inf,x_sup,y_sup,n,lb1,lb2);

   res_sup = cblas_ddot(n, dxs, inc, dys, inc);

   //setround(-1);

   res_inf = -cblas_ddot(n, dxi, inc, dyi, inc);

   SetIm(res, interval(res_inf,res_sup));

   setround(rnd);
}

inline void blasdot(const civector& x, const civector& y, cinterval& res) {
   const int n = VecLen(x);
   const int inc=1;
   const int lb1 = Lb(x);
   const int lb2 = Lb(y);
   int rnd = getround();

   cvector xmid(n),xrad(n),ymid(n),yrad(n);
   real crad_re = 0.0, crad_im = 0.0;
   complex res_inf(0.0,0.0),res_sup(0.0,0.0);
   real tmp1,tmp2,tmp3,tmp4;

   double* dxmid = xmid.to_blas_array();
   double* dymid = ymid.to_blas_array();

   setround(1);

   for(int i=0 ; i<n ; i++) {
      xmid[i+1] = Inf(x[lb1+i]) + 0.5*(Sup(x[lb1+i])-Inf(x[lb1+i]));
      ymid[i+1] = Inf(y[lb2+i]) + 0.5*(Sup(y[lb2+i])-Inf(y[lb2+i]));				
      xrad[i+1] = xmid[i+1] - Inf(x[i+lb1]);
      yrad[i+1] = ymid[i+1] - Inf(y[i+lb2]);

      tmp1 = abs(Re(ymid[i+1])) + Re(yrad[i+1]);
      tmp2 = abs(Im(ymid[i+1])) + Im(yrad[i+1]);
      tmp3 = abs(Re(xmid[i+1]));
      tmp4 = abs(Im(xmid[i+1]));

      crad_re += Re(xrad[i+1]) * tmp1 + tmp3 * Re(yrad[i+1]);
      crad_re += Im(xrad[i+1]) * tmp2 + tmp4 * Im(yrad[i+1]);

      crad_im += Re(xrad[i+1]) * tmp2 + tmp3 * Im(yrad[i+1]);
      crad_im += Im(xrad[i+1]) * tmp1 + tmp4 * Re(yrad[i+1]);
   }

   cblas_zdotu_sub(n, dxmid, inc, dymid, inc, &res_sup);

   res_sup += complex(crad_re,crad_im);

   setround(-1);

   cblas_zdotu_sub(n, dxmid, inc, dymid, inc, &res_inf);

   res_inf -= complex(crad_re,crad_im);

   setround(rnd);
   
   UncheckedSetInf(res,res_inf);
   UncheckedSetSup(res,res_sup);

}


inline void blasdot(const civector& x, const cvector& y, cinterval& res) {
   const int n = VecLen(x);
   const int inc=1;
   const int lb1 = Lb(x);
   int rnd = getround();

   cvector xmid(n),xrad(n);
   real crad_re = 0.0, crad_im = 0.0;
   complex res_inf(0.0,0.0),res_sup(0.0,0.0);
   real tmp1,tmp2,tmp3,tmp4;

   double* dxmid = xmid.to_blas_array();
   double* dy    = y.to_blas_array();

   setround(1);

   for(int i=0 ; i<n ; i++) {
      xmid[i+1] = Inf(x[lb1+i]) + 0.5*(Sup(x[lb1+i])-Inf(x[lb1+i]));
      xrad[i+1] = xmid[i+1] - Inf(x[i+lb1]);

      tmp1 = abs(Re(y[i+1]));
      tmp2 = abs(Im(y[i+1]));

      crad_re += Re(xrad[i+1]) * tmp1;
      crad_re += Im(xrad[i+1]) * tmp2;

      crad_im += Re(xrad[i+1]) * tmp2;
      crad_im += Im(xrad[i+1]) * tmp1;
   }

   cblas_zdotu_sub(n, dxmid, inc, dy, inc, &res_sup);

   res_sup += complex(crad_re,crad_im);

   setround(-1);

   cblas_zdotu_sub(n, dxmid, inc, dy, inc, &res_inf);

   res_inf -= complex(crad_re,crad_im);

   setround(rnd);
   
   UncheckedSetInf(res,res_inf);
   UncheckedSetSup(res,res_sup);

}

inline void blasdot(const cvector& x, const civector& y, cinterval& res) {
   const int n = VecLen(x);
   const int inc=1;
   const int lb1 = Lb(x);
   int rnd = getround();

   cvector ymid(n),yrad(n);
   real crad_re = 0.0, crad_im = 0.0;
   complex res_inf(0.0,0.0),res_sup(0.0,0.0);
   real tmp1,tmp2,tmp3,tmp4;

   double* dx    = x.to_blas_array();
   double* dymid = ymid.to_blas_array();

   setround(1);

   for(int i=0 ; i<n ; i++) {
      ymid[i+1] = Inf(y[lb1+i]) + 0.5*(Sup(y[lb1+i])-Inf(y[lb1+i]));
      yrad[i+1] = ymid[i+1] - Inf(y[i+lb1]);

      tmp1 = abs(Re(x[i+1]));
      tmp2 = abs(Im(x[i+1]));

      crad_re += Re(yrad[i+1]) * tmp1;
      crad_re += Im(yrad[i+1]) * tmp2;

      crad_im += Re(yrad[i+1]) * tmp2;
      crad_im += Im(yrad[i+1]) * tmp1;
   }

   cblas_zdotu_sub(n, dx, inc, dymid, inc, &res_sup);

   res_sup += complex(crad_re,crad_im);

   setround(-1);

   cblas_zdotu_sub(n, dx, inc, dymid, inc, &res_inf);

   res_inf -= complex(crad_re,crad_im);

   setround(rnd);
   
   UncheckedSetInf(res,res_inf);
   UncheckedSetSup(res,res_sup);

}


inline void blasdot(const civector& x, const ivector& y, cinterval& res) {
   const int n = VecLen(x);
   const int inc=1, inc2=2;
   const int lb1 = Lb(x);
   const int lb2 = Lb(y);
   int rnd = getround();

   cvector xmid(n),xrad(n);
   rvector ymid(n),yrad(n);
   real crad_re = 0.0, crad_im = 0.0;
   complex res_inf(0.0,0.0),res_sup(0.0,0.0);
   real tmp1,tmp2,tmp3,tmp4;

   double* dxmid = xmid.to_blas_array();
   double* dymid = ymid.to_blas_array();

   setround(1);

   for(int i=0 ; i<n ; i++) {
      xmid[i+1] = Inf(x[lb1+i]) + 0.5*(Sup(x[lb1+i])-Inf(x[lb1+i]));
      ymid[i+1] = Inf(y[lb2+i]) + 0.5*(Sup(y[lb2+i])-Inf(y[lb2+i]));				
      xrad[i+1] = xmid[i+1] - Inf(x[i+lb1]);
      yrad[i+1] = ymid[i+1] - Inf(y[i+lb1]);

      tmp1 = abs(ymid[i+1]) + yrad[i+1];
      tmp2 = abs(ymid[i+1]) + yrad[i+1];
      tmp3 = abs(Re(xmid[i+1]));
      tmp4 = abs(Im(xmid[i+1]));

      crad_re += Re(xrad[i+1]) * tmp1 + tmp3 * yrad[i+1];
      crad_re += Im(xrad[i+1]) * tmp2 + tmp4 * yrad[i+1];

      crad_im += Re(xrad[i+1]) * tmp2 + tmp3 * yrad[i+1];
      crad_im += Im(xrad[i+1]) * tmp1 + tmp4 * yrad[i+1];
   }

   SetRe(res_sup, cblas_ddot(n, dxmid, inc2, dymid, inc));
   SetIm(res_sup, cblas_ddot(n, dxmid+1, inc2, dymid, inc));

   res_sup += complex(crad_re,crad_im);

   setround(-1);

   SetRe(res_sup, cblas_ddot(n, dxmid, inc2, dymid, inc));
   SetIm(res_sup, cblas_ddot(n, dxmid+1, inc2, dymid, inc));

   res_inf -= complex(crad_re,crad_im);

   setround(rnd);

   UncheckedSetInf(res,res_inf);
   UncheckedSetSup(res,res_sup);
}


inline void blasdot(const ivector& x, const civector& y, cinterval& res) {
   const int n = VecLen(x);
   const int inc=1, inc2=2;
   const int lb1 = Lb(x);
   const int lb2 = Lb(y);
   int rnd = getround();

   rvector xmid(n),xrad(n);
   cvector ymid(n),yrad(n);
   real crad_re = 0.0, crad_im = 0.0;
   complex res_inf(0.0,0.0),res_sup(0.0,0.0);
   real tmp1,tmp2,tmp3,tmp4;

   double* dxmid = xmid.to_blas_array();
   double* dymid = ymid.to_blas_array();

   setround(1);

   for(int i=0 ; i<n ; i++) {
      xmid[i+1] = Inf(x[lb1+i]) + 0.5*(Sup(x[lb1+i])-Inf(x[lb1+i]));
      ymid[i+1] = Inf(y[lb2+i]) + 0.5*(Sup(y[lb2+i])-Inf(y[lb2+i]));				
      xrad[i+1] = xmid[i+1] - Inf(x[i+lb1]);
      yrad[i+1] = ymid[i+1] - Inf(y[i+lb1]);

      tmp1 = abs(Re(ymid[i+1])) + Re(yrad[i+1]);
      tmp2 = abs(Im(ymid[i+1])) + Im(yrad[i+1]);
      tmp3 = abs(xmid[i+1]);
      tmp4 = abs(xmid[i+1]);

      crad_re += xrad[i+1] * tmp1 + tmp3 * Re(yrad[i+1]);
      crad_re += xrad[i+1] * tmp2 + tmp4 * Im(yrad[i+1]);

      crad_im += xrad[i+1] * tmp2 + tmp3 * Im(yrad[i+1]);
      crad_im += xrad[i+1] * tmp1 + tmp4 * Re(yrad[i+1]);
   }

   SetRe(res_sup, cblas_ddot(n, dxmid, inc, dymid, inc2));
   SetIm(res_sup, cblas_ddot(n, dxmid, inc, dymid+1, inc2));

   res_sup += complex(crad_re,crad_im);

   setround(-1);

   SetRe(res_sup, cblas_ddot(n, dxmid, inc, dymid, inc2));
   SetIm(res_sup, cblas_ddot(n, dxmid, inc, dymid+1, inc2));

   res_inf -= complex(crad_re,crad_im);

   setround(rnd);

   UncheckedSetInf(res,res_inf);
   UncheckedSetSup(res,res_sup);
}
#endif
#endif
/***************************************************************************/


#ifdef _CXSC_BLAS_RMATRIX
#ifndef _CXSC_BLAS_RMATVEC_INC
#define _CXSC_BLAS_RMATVEC_INC
inline void blasmvmul(const rmatrix& A, const rvector& x, rvector& r) {
   
   double* DA = (double*) A.to_blas_array();  
   double* dx = (double*) x.to_blas_array();
   double* dr = (double*) r.to_blas_array();    
   
   const int m = Ub(A,1) - Lb(A,1) + 1;
   const int n = Ub(A,2) - Lb(A,2) + 1;
   const double alpha = 1.0, beta = 0.0;
   const int inc = 1;

   cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n,
               alpha, DA, n, dx, inc, beta, dr, inc);
}
#endif
#endif

#ifdef _CXSC_BLAS_RMATRIX
#ifdef _CXSC_BLAS_IVECTOR
#ifndef _CXSC_BLAS_RIMATVEC_INC
#define _CXSC_BLAS_RIMATVEC_INC
inline void blasmvmul(const rmatrix& A, const rvector& x, ivector& r) {
   int rnd = getround();
   double* DA = (double*) A.to_blas_array();  
   double* dx = (double*) x.to_blas_array();
   
   const int m = Ub(A,1) - Lb(A,1) + 1;
   const int n = Ub(A,2) - Lb(A,2) + 1;
   const double alpha = 1.0, beta = 0.0;
   const int inc = 1;
   double* dr = new double[n];    

   setround(1);
   cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n,
               alpha, DA, n, dx, inc, beta, dr, inc);
   for(int i=0 ; i<n ; i++)
     UncheckedSetSup(r[i+Lb(r)], dr[i]);

   setround(-1);
   cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n,
               alpha, DA, n, dx, inc, beta, dr, inc);
   for(int i=0 ; i<n ; i++)
     UncheckedSetInf(r[i+Lb(r)], dr[i]);

   delete[] dr;
   setround(rnd);
}

inline void blasmvmul(const rmatrix& A, const ivector& x, ivector& r) {
   const int n = VecLen(x);
   const int m = Ub(A,1)-Lb(A,1)+1;
   const int inc=1;
   const int lb1 = Lb(A,2);
   const int lb2 = Lb(x);
   int rnd = getround();

   rvector x_inf(n),x_sup(n),y_inf(n),y_sup(n);

   double* dxi = x_inf.to_blas_array();
   double* dxs = x_sup.to_blas_array();
   double* dyi = y_inf.to_blas_array();
   double* dys = y_sup.to_blas_array();
   double res_inf,res_sup;

   setround(1);
 
   for(int i=1 ; i<=m ; i++) {
     bsort(A[i+Lb(A,1)-1],x,x_inf,y_inf,x_sup,y_sup,n,lb1,lb2);
     res_sup = cblas_ddot(n, dxs, inc, dys, inc);
     UncheckedSetSup(r[i], res_sup);
    // setround(-1);
     res_inf = cblas_ddot(n, dxi, inc, dyi, inc);
     UncheckedSetInf(r[i], -res_inf);
   }

   setround(rnd);
}
#endif
#endif
#endif

#ifdef _CXSC_BLAS_IMATRIX
#ifdef _CXSC_BLAS_IVECTOR
#ifndef _CXSC_BLAS_IRMATVEC_INC
#define _CXSC_BLAS_IRMATVEC_INC
inline void blasmvmul(const imatrix& A, const rvector& x, ivector& r) {
   const int n = VecLen(x);
   const int m = Ub(A,1)-Lb(A,1)+1;
   const int inc=1;
   const int lb1 = Lb(A,2);
   const int lb2 = Lb(x);
   int rnd = getround();

   rvector x_inf(n),x_sup(n),y_inf(n),y_sup(n);

   double* dxi = x_inf.to_blas_array();
   double* dxs = x_sup.to_blas_array();
   double* dyi = y_inf.to_blas_array();
   double* dys = y_sup.to_blas_array();
   double res_inf,res_sup;

   setround(1);
 
   for(int i=1 ; i<=m ; i++) {
     bsort(A[i+Lb(A,1)-1],x,x_inf,y_inf,x_sup,y_sup,n,lb1,lb2);
     res_sup = cblas_ddot(n, dxs, inc, dys, inc);
     UncheckedSetSup(r[i], res_sup);
    // setround(-1);
     res_inf = cblas_ddot(n, dxi, inc, dyi, inc);
     UncheckedSetInf(r[i], -res_inf);
   }

   setround(rnd);

}
#endif
#endif
#endif

#ifdef _CXSC_BLAS_IMATRIX
#ifdef _CXSC_BLAS_IVECTOR
#ifndef _CXSC_BLAS_IIMATVEC_INC
#define _CXSC_BLAS_IIMATVEC_INC
inline void blasmvmul(const imatrix& A, const ivector& x, ivector& r) {
   const int n = VecLen(x);
   const int m = Ub(A,1)-Lb(A,1)+1;
   const int inc=1;
   const int lb1 = Lb(A,2);
   const int lb2 = Lb(x);
   int rnd = getround();

   rvector x_inf(n),x_sup(n),y_inf(n),y_sup(n);

   double* dxi = x_inf.to_blas_array();
   double* dxs = x_sup.to_blas_array();
   double* dyi = y_inf.to_blas_array();
   double* dys = y_sup.to_blas_array();
   double res_inf,res_sup;

   setround(1);
 
   for(int i=1 ; i<=m ; i++) {
     bsort(A[i+Lb(A,1)-1],x,x_inf,y_inf,x_sup,y_sup,n,lb1,lb2);
     res_sup = cblas_ddot(n, dxs, inc, dys, inc);
     UncheckedSetSup(r[i], res_sup);
    // setround(-1);
     res_inf = cblas_ddot(n, dxi, inc, dyi, inc);
     UncheckedSetInf(r[i], -res_inf);
   }

   setround(rnd);

}
#endif
#endif
#endif

#ifdef _CXSC_BLAS_CMATRIX
#ifdef _CXSC_BLAS_CIVECTOR
#ifndef _CXSC_BLAS_CIMATVEC_INC
#define _CXSC_BLAS_CIMATVEC_INC
inline void blasmvmul(const cmatrix& A, const ivector& x, civector& r) {
   const int m = Ub(A,1)-Lb(A,1)+1;

   for(int i=0 ; i<m ; i++)
     blasdot(A[i+Lb(A,1)], x, r[i+Lb(r)]);
}
#endif
#endif
#endif

#ifdef _CXSC_BLAS_IMATRIX
#ifdef _CXSC_BLAS_CIVECTOR
#ifndef _CXSC_BLAS_ICMATVEC_INC
#define _CXSC_BLAS_ICMATVEC_INC
inline void blasmvmul(const imatrix& A, const cvector& x, civector& r) {
   const int m = Ub(A,1)-Lb(A,1)+1;

   for(int i=0 ; i<m ; i++)
     blasdot(A[i+Lb(A,1)], x, r[i+Lb(r)]);
}
#endif
#endif
#endif

#ifdef _CXSC_BLAS_RMATRIX
#ifdef _CXSC_BLAS_CIVECTOR
#ifndef _CXSC_BLAS_RCIMATVEC_INC
#define _CXSC_BLAS_RCIMATVEC_INC
inline void blasmvmul(const rmatrix& A, const civector& x, civector& r) {
   const int m = Ub(A,1)-Lb(A,1)+1;

   for(int i=0 ; i<m ; i++)
     blasdot(A[i+Lb(A,1)], x, r[i+Lb(r)]);
}
#endif
#endif
#endif

#ifdef _CXSC_BLAS_CIMATRIX
#ifdef _CXSC_BLAS_CIVECTOR
#ifndef _CXSC_BLAS_CIRMATVEC_INC
#define _CXSC_BLAS_CIRMATVEC_INC
inline void blasmvmul(const cimatrix& A, const rvector& x, civector& r) {
   const int m = Ub(A,1)-Lb(A,1)+1;

   for(int i=0 ; i<m ; i++)
     blasdot(A[i+Lb(A,1)], x, r[i+Lb(r)]);
}
#endif
#endif
#endif

#ifdef _CXSC_BLAS_CMATRIX
#ifdef _CXSC_BLAS_CIVECTOR
#ifndef _CXSC_BLAS_CCIMATVEC_INC
#define _CXSC_BLAS_CCIMATVEC_INC
inline void blasmvmul(const cmatrix& A, const civector& x, civector& r) {
   const int m = Ub(A,1)-Lb(A,1)+1;

   for(int i=0 ; i<m ; i++)
     blasdot(A[i+Lb(A,1)], x, r[i+Lb(r)]);
}

inline void blasmvmul(const cmatrix& A, const cvector& x, civector& r) {
   ivector tmp1, tmp2;
   blasmvmul(Re(A),Re(x),tmp1);
   blasmvmul(Im(A),Im(x),tmp2);
   SetRe(r,tmp1-tmp2);
   blasmvmul(Re(A),Im(x),tmp1);
   blasmvmul(Im(A),Re(x),tmp2);
   SetIm(r,tmp1+tmp2);
}

inline void blasmvmul(const cmatrix& A, const rvector& x, civector& r) {
   int n = VecLen(r);
   ivector re(n),im(n);
   
   blasmvmul(Re(A),x,re);
   blasmvmul(Im(A),x,im);

   SetRe(r,re);
   SetIm(r,im);
}

inline void blasmvmul(const rmatrix& A, const cvector& x, civector& r) {
   int n = VecLen(r);
   ivector re(n),im(n);
   
   blasmvmul(A,Re(x),re);
   blasmvmul(A,Im(x),im);

   SetRe(r,re);
   SetIm(r,im);
}
#endif
#endif
#endif

#ifdef _CXSC_BLAS_CIMATRIX
#ifdef _CXSC_BLAS_CIVECTOR
#ifndef _CXSC_BLAS_CICMATVEC_INC
#define _CXSC_BLAS_CICMATVEC_INC
inline void blasmvmul(const cimatrix& A, const cvector& x, civector& r) {
   const int m = Ub(A,1)-Lb(A,1)+1;

   for(int i=0 ; i<m ; i++)
     blasdot(A[i+Lb(A,1)], x, r[i+Lb(r)]);
}
#endif
#endif
#endif

#ifdef _CXSC_BLAS_IMATRIX
#ifdef _CXSC_BLAS_CIVECTOR
#ifndef _CXSC_BLAS_ICIMATVEC_INC
#define _CXSC_BLAS_ICIMATVEC_INC
inline void blasmvmul(const imatrix& A, const civector& x, civector& r) {
   const int m = Ub(A,1)-Lb(A,1)+1;

   for(int i=0 ; i<m ; i++)
     blasdot(A[i+Lb(A,1)], x, r[i+Lb(r)]);
}
#endif
#endif
#endif

#ifdef _CXSC_BLAS_CIMATRIX
#ifdef _CXSC_BLAS_CIVECTOR
#ifndef _CXSC_BLAS_CIIMATVEC_INC
#define _CXSC_BLAS_CIIMATVEC_INC
inline void blasmvmul(const cimatrix& A, const ivector& x, civector& r) {
   const int m = Ub(A,1)-Lb(A,1)+1;

   for(int i=0 ; i<m ; i++)
     blasdot(A[i+Lb(A,1)], x, r[i+Lb(r)]);
}
#endif
#endif
#endif

#ifdef _CXSC_BLAS_CIMATRIX
#ifdef _CXSC_BLAS_CIVECTOR
#ifndef _CXSC_BLAS_CICIMATVEC_INC
#define _CXSC_BLAS_CICIMATVEC_INC
inline void blasmvmul(const cimatrix& A, const civector& x, civector& r) {
   const int m = Ub(A,1)-Lb(A,1)+1;

   for(int i=0 ; i<m ; i++)
     blasdot(A[i+Lb(A,1)], x, r[i+Lb(r)]);
}
#endif
#endif
#endif

#ifdef _CXSC_BLAS_CMATRIX
#ifdef _CXSC_BLAS_CVECTOR
#ifndef _CXSC_BLAS_CCMATVEC_INC
#define _CXSC_BLAS_CCMATVEC_INC
inline void blasmvmul(const cmatrix& A, const cvector& x, cvector& r) {
   
   double* DA = (double*) A.to_blas_array();  
   double* dx = (double*) x.to_blas_array();
   double* dr = (double*) r.to_blas_array();    
   
   const int m = Ub(A,1) - Lb(A,1) + 1;
   const int n = Ub(A,2) - Lb(A,2) + 1;
   const complex alpha(1.0,0.0);
   const complex beta(0.0,0.0);
   const int inc = 1;

   cblas_zgemv(CblasRowMajor, CblasNoTrans, m, n,
               &alpha, DA, n, dx, inc, &beta, dr, inc);
}
#endif
#endif
#endif

#ifdef _CXSC_BLAS_RMATRIX
#ifdef _CXSC_BLAS_CVECTOR
#ifndef _CXSC_BLAS_RCMATVEC_INC
#define _CXSC_BLAS_RCMATVEC_INC
inline void blasmvmul(const rmatrix& A, const cvector& x, cvector& r) {
   
   double* DA = (double*) A.to_blas_array();  
   double* dx = (double*) x.to_blas_array();
   double* dr = (double*) r.to_blas_array();    
   
   const int m = Ub(A,1) - Lb(A,1) + 1;
   const int n = Ub(A,2) - Lb(A,2) + 1;
   const double alpha = 1.0;
   const double beta = 0.0;
   const int inc = 2;

   cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n,
               alpha, DA, m, dx, inc, beta, dr, inc);

   cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n,
               alpha, DA, n, dx+1, inc, beta, dr+1, inc);
}
#endif
#endif
#endif

#ifdef _CXSC_BLAS_CMATRIX
#ifdef _CXSC_BLAS_CVECTOR
#ifndef _CXSC_BLAS_CRMATVEC_INC
#define _CXSC_BLAS_CRMATVEC_INC
inline void blasmvmul(const cmatrix& A, const rvector& x, cvector& r) {
   cvector tmp(x);

   double* DA = (double*) A.to_blas_array();  
   double* dx = (double*) tmp.to_blas_array();
   double* dr = (double*) r.to_blas_array();    
   
   const int m = Ub(A,1) - Lb(A,1) + 1;
   const int n = Ub(A,2) - Lb(A,2) + 1;
   const complex alpha(1.0,0.0);
   const complex beta(0.0,0.0);
   const int inc = 1;

   cblas_zgemv(CblasRowMajor, CblasNoTrans, m, n,
               &alpha, DA, n, dx, inc, &beta, dr, inc);
}
#endif
#endif
#endif

/***************************************************************************/

#ifdef _CXSC_BLAS_RMATRIX
#ifndef _CXSC_BLAS_RMATRIX_INC
#define _CXSC_BLAS_RMATRIX_INC
inline void blasmatmul(const rmatrix &A, const rmatrix &B, rmatrix &C) {
   
   double* DA = (double*) A.to_blas_array();
   double* DB = (double*) B.to_blas_array();
   double* DC = (double*) C.to_blas_array();       

   const int m = Ub(A,1) - Lb(A,1) + 1;
   const int n = Ub(A,2) - Lb(A,2) + 1;
   const int o = Ub(C,2) - Lb(C,2) + 1;
   const double alpha = 1.0, beta = 0.0;
   
   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, o,
               n, alpha, DA, n, DB, o, beta, DC, o);    
}
#endif
#endif

#ifdef _CXSC_BLAS_IMATRIX
#ifndef _CXSC_BLAS_IMATRIX_INC
#define _CXSC_BLAS_IMATRIX_INC
inline void blasmatmul(const rmatrix &A, const rmatrix &B, imatrix &C) {
   int rnd = getround();   
   double* DA = (double*) A.to_blas_array();
   double* DB = (double*) B.to_blas_array();
   int ind;

   const int m = Ub(A,1) - Lb(A,1) + 1;
   const int n = Ub(A,2) - Lb(A,2) + 1;
   const int o = Ub(C,2) - Lb(C,2) + 1;
   double* DC = new double[m*o];
   const double alpha = 1.0, beta = 0.0;
   
   setround(1);
   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, o,
               n, alpha, DA, n, DB, o, beta, DC, o);    

   for(int i=0 ; i<m ; i++) {
      for(int j=0 ; j<o ; j++) {
         ind = i*o+j;
         UncheckedSetSup(C[i+Lb(C,1)][j+Lb(C,2)], DC[ind]);
      }
   }

   setround(-1);
   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, o,
               n, alpha, DA, n, DB, o, beta, DC, o);    

   for(int i=0 ; i<m ; i++) {
      for(int j=0 ; j<o ; j++) {
         ind = i*o+j;
         UncheckedSetInf(C[i+Lb(C,1)][j+Lb(C,2)], DC[ind]);
      }
   }

   delete[] DC;

   setround(rnd);
}


inline void blasmatmul(const imatrix &A, const imatrix &B, imatrix &C) {
   const int m = Ub(A,1) - Lb(A,1) + 1;
   const int n = Ub(A,2) - Lb(A,2) + 1;
   const int o = Ub(C,2) - Lb(C,2) + 1;

   int lbC_1 = Lb(C,1); int lbC_2 = Lb(C,2);
   Resize(C);

   double* Amid = new double[m*n];
   double* Aabs = new double[m*n];
   double* Arad = new double[m*n];
   double* Bmid = new double[n*o];
   double* Babs = new double[n*o];
   double* Brad = new double[n*o];
   double* C1   = new double[m*o];
   double* C2   = new double[m*o];
   int rnd = getround();

   const double alpha = 1.0, beta = 0.0;

   setround(1);

   int ind;

   //Compute mid and rad of A
   for(int i=Lb(A,1) ; i<=Ub(A,1) ; i++) {
      for(int j=Lb(A,2) ; j<=Ub(A,2) ; j++) {
         ind = (i-Lb(A,1))*n+(j-Lb(A,2));
         Amid[ind] = _double(Inf(A[i][j]) + 0.5*(Sup(A[i][j]) - Inf(A[i][j])));
         Arad[ind] = _double(Amid[ind] - Inf(A[i][j]));
         Aabs[ind] = fabs(Amid[ind]);
      }
   }

   //Compute mid and rad of B
   for(int i=Lb(B,1) ; i<=Ub(B,1) ; i++) {
      for(int j=Lb(B,2) ; j<=Ub(B,2) ; j++) {
         ind = (i-Lb(B,1))*o+(j-Lb(B,2));
         Bmid[ind] = _double(Inf(B[i][j]) + 0.5*(Sup(B[i][j]) - Inf(B[i][j])));
         Brad[ind] = _double(Bmid[ind] - Inf(B[i][j]));
         Babs[ind] = fabs(Bmid[ind]);
      }
   }

   setround(-1);

   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, o,
               n, alpha, Amid, n, Bmid, o, beta, C1, o);    

   setround(1);

   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, o,
               n, alpha, Amid, n, Bmid, o, beta, C2, o);    

   delete[] Amid;
   delete[] Bmid;

   double* Cmid = new double[m*o];
   double* Crad = new double[m*o];

   for(int i=0 ; i<m ; i++) {
      for(int j=0 ; j<o ; j++) {
         ind = i*o+j;
         Cmid[ind] = C1[ind] + 0.5*(C2[ind] - C1[ind]);
         Crad[ind] = Cmid[ind] - C1[ind];
      }
   }

   for(int i=0 ; i<n ; i++) {
      for(int j=0 ; j<o ; j++) {
         ind = i*o+j;
         Babs[ind] += Brad[ind];
      }
   }
   
   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, o,
               n, alpha, Arad, n, Babs, o, beta, C1, o);    

   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, o,
               n, alpha, Aabs, n, Brad, o, beta, C2, o);    

   delete[] Aabs;
   delete[] Arad;
   delete[] Babs;
   delete[] Brad;

   Resize(C, lbC_1, lbC_1+m-1, lbC_2, lbC_2+o-1);
   //C = imatrix(lbC_1, lbC_1+m-1, lbC_2, lbC_2+o-1);

   for(int i=0 ; i<m ; i++) {
      for(int j=0 ; j<o ; j++) {
         ind = i*o+j;
         Crad[ind] += C1[ind] + C2[ind];
         UncheckedSetSup(C[i+Lb(C,1)][j+Lb(C,2)], Cmid[ind] + Crad[ind]);
      }
   }

   setround(-1);

   for(int i=0 ; i<m ; i++) {
      for(int j=0 ; j<o ; j++) {
         ind = i*o+j;
         UncheckedSetInf(C[i+Lb(C,1)][j+Lb(C,2)], Cmid[ind] - Crad[ind]);
      }
   }

   setround(rnd);

   delete[] C1;
   delete[] C2;
   delete[] Crad;
   delete[] Cmid;
   
}


inline void blasmatmul(const rmatrix &A, const imatrix &B, imatrix &C) {

   const int m = Ub(A,1) - Lb(A,1) + 1;
   const int n = Ub(A,2) - Lb(A,2) + 1;
   const int o = Ub(B,2) - Lb(B,2) + 1;

   int lb1_C = Lb(C,1);
   int lb2_C = Lb(C,2);
   Resize(C);

   double* DA = (double*)A.to_blas_array();
   double* Bmid = new double[n*o];
   double* Brad = new double[n*o];
   double* C1   = new double[m*o];
   double* C2   = new double[m*o];
   int rnd = getround();

   const double alpha = 1.0;
   double beta = 0.0;

   setround(1);

   int ind;

   //Compute mid and rad of B
   for(int i=Lb(B,1) ; i<=Ub(B,1) ; i++) {
      for(int j=Lb(B,2) ; j<=Ub(B,2) ; j++) {
         ind = (i-Lb(B,1))*o+(j-Lb(B,2));
         Bmid[ind] = _double(Inf(B[i][j]) + 0.5*(Sup(B[i][j]) - Inf(B[i][j])));
         Brad[ind] = _double(Bmid[ind] - Inf(B[i][j]));
      }
   }

   setround(-1);

   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, o,
               n, alpha, DA, n, Bmid, o, beta, C1, o);    

   setround(1);

   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, o,
               n, alpha, DA, n, Bmid, o, beta, C2, o);    

   delete[] Bmid;

   double* Cmid = new double[m*o];
   for(int i=0 ; i<m ; i++) {
      for(int j=0 ; j<o ; j++) {
         ind = i*o+j;
         Cmid[ind] = C1[ind] + 0.5*(C2[ind] - C1[ind]);
      }
   }

   delete[] C2;

   double* Crad = new double[m*o];
   for(int i=0 ; i<m ; i++) {
      for(int j=0 ; j<o ; j++) {
         ind = i*o+j;
         Crad[ind] = Cmid[ind] - C1[ind];
      }
   }

   delete[] C1;

   double* Aabs = new double[m*n];
   //Compute abs(A)
   for(int i=Lb(A,1) ; i<=Ub(A,1) ; i++) {
      for(int j=Lb(A,2) ; j<=Ub(A,2) ; j++) {
         ind = (i-Lb(A,1))*n+(j-Lb(A,2));
         Aabs[ind] = fabs(DA[ind]);
      }
   }

   beta = 1.0;

   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, o,
               n, alpha, Aabs, n, Brad, o, beta, Crad, o);    

   delete[] Aabs;
   delete[] Brad;

   Resize(C, lb1_C, lb1_C+m-1, lb2_C, lb2_C+o-1);
   //C = imatrix(lb1_C, lb1_C+m-1, lb2_C, lb2_C+o-1);

   for(int i=0 ; i<m ; i++) {
      for(int j=0 ; j<o ; j++) {
         ind = i*o+j;
         UncheckedSetSup(C[i+Lb(C,1)][j+Lb(C,2)], Cmid[ind] + Crad[ind]);
      }
   }

   setround(-1);

   for(int i=0 ; i<m ; i++) {
      for(int j=0 ; j<o ; j++) {
         ind = i*o+j;
         UncheckedSetInf(C[i+Lb(C,1)][j+Lb(C,2)], Cmid[ind] - Crad[ind]);
      }
   }

   setround(rnd);

   delete[] Crad;
   delete[] Cmid;
}


inline void blasmatmul(const imatrix &A, const rmatrix &B, imatrix &C) {
   
   const int m = Ub(A,1) - Lb(A,1) + 1;
   const int n = Ub(A,2) - Lb(A,2) + 1;
   const int o = Ub(C,2) - Lb(C,2) + 1;

   int lb1_C = Lb(C,1);
   int lb2_C = Lb(C,2);
   Resize(C);

   double* DB   = (double*) B.to_blas_array();
   double* Amid = new double[m*n];
   double* Arad = new double[m*n];
   double* C1   = new double[m*o];
   double* C2   = new double[m*o];
   int rnd = getround();

   const double alpha = 1.0;
   double beta = 0.0;

   setround(1);

   int ind;

   //Compute mid and rad of A
   for(int i=Lb(A,1) ; i<=Ub(A,1) ; i++) {
      for(int j=Lb(A,2) ; j<=Ub(A,2) ; j++) {
         ind = (i-Lb(A,1))*n+(j-Lb(A,2));
         Amid[ind] = _double(Inf(A[i][j]) + 0.5*(Sup(A[i][j]) - Inf(A[i][j])));
         Arad[ind] = _double(Amid[ind] - Inf(A[i][j]));
      }
   }

   setround(-1);

   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, o,
               n, alpha, Amid, n, DB, o, beta, C1, o);    

   setround(1);

   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, o,
               n, alpha, Amid, n, DB, o, beta, C2, o);    

   delete[] Amid;

   double* Cmid = new double[m*o];
   for(int i=0 ; i<m ; i++) {
      for(int j=0 ; j<o ; j++) {
         ind = i*o+j;
         Cmid[ind] = C1[ind] + 0.5*(C2[ind] - C1[ind]);
      }
   }

   delete[] C2;

   double* Crad = new double[m*o];
   for(int i=0 ; i<m ; i++) {
      for(int j=0 ; j<o ; j++) {
         ind = i*o+j;
         Crad[ind] = Cmid[ind] - C1[ind];
      }
   }

   delete[] C1;

   double* Babs = new double[n*o];
   //Compute abs of B
   for(int i=Lb(B,1) ; i<=Ub(B,1) ; i++) {
      for(int j=Lb(B,2) ; j<=Ub(B,2) ; j++) {
         ind = (i-Lb(B,1))*o+(j-Lb(B,2));
         Babs[ind] = fabs(DB[ind]);
      }
   }

   beta = 1.0;

   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, o,
               n, alpha, Arad, n, Babs, o, beta, Crad, o);    

   delete[] Arad;
   delete[] Babs;

   Resize(C, lb1_C, lb1_C+m-1, lb2_C, lb2_C+o-1);
   //C = imatrix(lb1_C, lb1_C+m-1, lb2_C, lb2_C+o-1);

   for(int i=0 ; i<m ; i++) {
      for(int j=0 ; j<o ; j++) {
         ind = i*o+j;
         UncheckedSetSup(C[i+Lb(C,1)][j+Lb(C,2)], Cmid[ind] + Crad[ind]);
      }
   }

   setround(-1);

   for(int i=0 ; i<m ; i++) {
      for(int j=0 ; j<o ; j++) {
         ind = i*o+j;
         UncheckedSetInf(C[i+Lb(C,1)][j+Lb(C,2)], Cmid[ind] - Crad[ind]);
      }
   }

   setround(rnd);

   delete[] Crad;
   delete[] Cmid;
}
#endif
#endif

#ifdef _CXSC_BLAS_CMATRIX
#ifndef _CXSC_BLAS_CMATRIX_INC
#define _CXSC_BLAS_CMATRIX_INC
inline void blasmatmul(const cmatrix &A, const cmatrix &B, cmatrix &C) {
   
   double* DA = (double*) A.to_blas_array();
   double* DB = (double*) B.to_blas_array();
   double* DC = (double*) C.to_blas_array();       

   const int m = Ub(A,1) - Lb(A,1) + 1;
   const int n = Ub(A,2) - Lb(A,2) + 1;
   const int o = Ub(C,2) - Lb(C,2) + 1;
   const complex alpha(1.0,0.0);
   const complex beta(0.0,0.0);
   
   cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, o,
               n, &alpha, DA, n, DB, o, &beta, DC, o);    
}


inline void blasmatmul(const cmatrix &A, const rmatrix &B, cmatrix &C) {
   cmatrix tmp(B);

   double* DA = (double*) A.to_blas_array();
   double* DB = (double*) tmp.to_blas_array();
   double* DC = (double*) C.to_blas_array();       

   const int m = Ub(A,1) - Lb(A,1) + 1;
   const int n = Ub(A,2) - Lb(A,2) + 1;
   const int o = Ub(C,2) - Lb(C,2) + 1;
   const complex alpha(1.0,0.0);
   const complex beta(0.0,0.0);

   cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, o,
               n, &alpha, DA, n, DB, o, &beta, DC, o);    
}

inline void blasmatmul(const rmatrix &A, const cmatrix &B, cmatrix &C) {
   cmatrix tmp(A);

   double* DA = (double*) tmp.to_blas_array();
   double* DB = (double*) B.to_blas_array();
   double* DC = (double*) C.to_blas_array();       

   const int m = Ub(A,1) - Lb(A,1) + 1;
   const int n = Ub(A,2) - Lb(A,2) + 1;
   const int o = Ub(C,2) - Lb(C,2) + 1;
   const complex alpha(1.0,0.0);
   const complex beta(0.0,0.0);

   cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, o,
               n, &alpha, DA, n, DB, o, &beta, DC, o);    
}
#endif
#endif

#ifdef _CXSC_BLAS_CIMATRIX
#ifndef _CXSC_BLAS_CIMATRIX_INC
#define _CXSC_BLAS_CIMATRIX_INC
inline void blasmatmul(const cmatrix &A, const cmatrix &B, cimatrix &C) {
   const int m = Ub(A,1) - Lb(A,1) + 1;
   const int o = Ub(C,2) - Lb(C,2) + 1;
   imatrix tmp1(m,o), tmp2(m,o);
   
   blasmatmul(Re(A),Re(B),tmp1);
   blasmatmul(Im(A),Im(B),tmp2);
   SetRe(C,tmp1-tmp2);
   blasmatmul(Re(A),Im(B),tmp1);
   blasmatmul(Im(A),Re(B),tmp2);
   SetIm(C,tmp1+tmp2);
}


inline void blasmatmul(const cmatrix &A, const imatrix &B, cimatrix &C) {
   const int m = Ub(A,1) - Lb(A,1) + 1;
   const int o = Ub(C,2) - Lb(C,2) + 1;

   imatrix T(m,o);

   blasmatmul(Re(A),B,T);

   SetRe(C, T);

   blasmatmul(Im(A),B,T);

   SetIm(C, T);
}

inline void blasmatmul(const imatrix &A, const cmatrix &B, cimatrix &C) {
   const int m = Ub(A,1) - Lb(A,1) + 1;
   const int o = Ub(C,2) - Lb(C,2) + 1;

   imatrix T(m,o);

   blasmatmul(A,Re(B),T);

   SetRe(C, T);

   blasmatmul(A,Im(B),T);

   SetIm(C, T);
}

inline void blasmatmul(const rmatrix &A, const cimatrix &B, cimatrix &C) {
   const int m = Ub(A,1) - Lb(A,1) + 1;
   const int o = Ub(C,2) - Lb(C,2) + 1;

   imatrix T(m,o);

   blasmatmul(A,Re(B),T);

   SetRe(C, T);

   blasmatmul(A,Im(B),T);

   SetIm(C, T);
}

inline void blasmatmul(const cimatrix &A, const rmatrix &B, cimatrix &C) {
   const int m = Ub(A,1) - Lb(A,1) + 1;
   const int o = Ub(C,2) - Lb(C,2) + 1;

   imatrix T(m,o);

   blasmatmul(Re(A),B,T);

   SetRe(C, T);

   blasmatmul(Im(A),B,T);

   SetIm(C, T);
}

inline void blasmatmul(const imatrix &A, const cimatrix &B, cimatrix &C) {
   const int m = Ub(A,1) - Lb(A,1) + 1;
   const int o = Ub(C,2) - Lb(C,2) + 1;

   imatrix T(m,o);

   blasmatmul(A,Re(B),T);

   SetRe(C, T);

   blasmatmul(A,Im(B),T);

   SetIm(C, T);
}

inline void blasmatmul(const cimatrix &A, const imatrix &B, cimatrix &C) {
   const int m = Ub(A,1) - Lb(A,1) + 1;
   const int o = Ub(C,2) - Lb(C,2) + 1;

   imatrix T(m,o);

   blasmatmul(Re(A),B,T);

   SetRe(C, T);

   blasmatmul(Im(A),B,T);

   SetIm(C, T);
}


inline void blasmatmul(const cimatrix &A, const cmatrix &B, cimatrix &C) {

   const int m = Ub(A,1) - Lb(A,1) + 1;
   const int o = Ub(C,2) - Lb(C,2) + 1;
   int rnd = getround();

   imatrix T(m,o);

   blasmatmul(Re(A),Re(B),T);

   SetRe(C, T);

   blasmatmul(-Im(A),Im(B),T);

   setround(-1);

   for(int i=Lb(C,1) ; i<=Ub(C,1) ; i++) {
      for(int j=Lb(C,2) ; j<=Ub(C,2) ; j++) {
         InfRe(C[i][j]) += Inf(T[i][j]);
      }
   }

   setround(1);

   for(int i=Lb(C,1) ; i<=Ub(C,1) ; i++) {
      for(int j=Lb(C,2) ; j<=Ub(C,2) ; j++) {
         SupRe(C[i][j]) += Sup(T[i][j]);
      }
   }

   setround(0);

   blasmatmul(Re(A),Im(B),T);

   SetIm(C, T);

   blasmatmul(Im(A),Re(B),T);

   setround(-1);

   for(int i=Lb(C,1) ; i<=Ub(C,1) ; i++) {
      for(int j=Lb(C,2) ; j<=Ub(C,2) ; j++) {
         InfIm(C[i][j]) += Inf(T[i][j]);
      }
   }

   setround(1);

   for(int i=Lb(C,1) ; i<=Ub(C,1) ; i++) {
      for(int j=Lb(C,2) ; j<=Ub(C,2) ; j++) {
         SupIm(C[i][j]) += Sup(T[i][j]);
      }
   }

   setround(rnd);
}

inline void blasmatmul(const cmatrix &A, const cimatrix &B, cimatrix &C) {
   int rnd = getround();

   imatrix T(Lb(C,1),Ub(C,1),Lb(C,2),Ub(C,2));

   blasmatmul(Re(A),Re(B),T);

   SetRe(C, T);

   blasmatmul(-Im(A),Im(B),T);

   setround(-1);

   for(int i=Lb(C,1) ; i<=Ub(C,1) ; i++) {
      for(int j=Lb(C,2) ; j<=Ub(C,2) ; j++) {
         InfRe(C[i][j]) += Inf(T[i][j]);
      }
   }

   setround(1);

   for(int i=Lb(C,1) ; i<=Ub(C,1) ; i++) {
      for(int j=Lb(C,2) ; j<=Ub(C,2) ; j++) {
         SupRe(C[i][j]) += Sup(T[i][j]);
      }
   }

   setround(0);

   blasmatmul(Re(A),Im(B),T);

   SetIm(C, T);

   blasmatmul(Im(A),Re(B),T);

   setround(-1);

   for(int i=Lb(C,1) ; i<=Ub(C,1) ; i++) {
      for(int j=Lb(C,2) ; j<=Ub(C,2) ; j++) {
         InfIm(C[i][j]) += Inf(T[i][j]);
      }
   }

   setround(1);

   for(int i=Lb(C,1) ; i<=Ub(C,1) ; i++) {
      for(int j=Lb(C,2) ; j<=Ub(C,2) ; j++) {
         SupIm(C[i][j]) += Sup(T[i][j]);
      }
   }

   setround(rnd);
}


inline void blasmatmul(const cimatrix &A, const cimatrix &B, cimatrix &C) {
   const int m = Ub(A,1) - Lb(A,1) + 1;
   const int o = Ub(C,2) - Lb(C,2) + 1;
   int rnd = getround();

   imatrix T(m,o);

   blasmatmul(Re(A),Re(B),T);

   SetRe(C, T);

   blasmatmul(-Im(A),Im(B),T);

   setround(-1);

   for(int i=Lb(C,1) ; i<=Ub(C,1) ; i++) {
      for(int j=Lb(C,2) ; j<=Ub(C,2) ; j++) {
         InfRe(C[i][j]) += Inf(T[i][j]);
      }
   }

   setround(1);

   for(int i=Lb(C,1) ; i<=Ub(C,1) ; i++) {
      for(int j=Lb(C,2) ; j<=Ub(C,2) ; j++) {
         SupRe(C[i][j]) += Sup(T[i][j]);
      }
   }

   setround(0);

   blasmatmul(Re(A),Im(B),T);

   SetIm(C, T);

   blasmatmul(Im(A),Re(B),T);

   setround(-1);

   for(int i=Lb(C,1) ; i<=Ub(C,1) ; i++) {
      for(int j=Lb(C,2) ; j<=Ub(C,2) ; j++) {
         InfIm(C[i][j]) += Inf(T[i][j]);
      }
   }

   setround(1);

   for(int i=Lb(C,1) ; i<=Ub(C,1) ; i++) {
      for(int j=Lb(C,2) ; j<=Ub(C,2) ; j++) {
         SupIm(C[i][j]) += Sup(T[i][j]);
      }
   }

   setround(rnd);
}

inline void blasmatmul(const rmatrix &A, const cmatrix &B, cimatrix &C) {
  int m = Ub(A,1)-Lb(A,1)+1;
  int n = Ub(A,2)-Lb(A,2)+1;
  imatrix re(m,n),im(m,n);

  blasmatmul(A,Re(B),re);
  blasmatmul(A,Im(B),im);

  SetRe(C,re);
  SetIm(C,im);
}

inline void blasmatmul(const cmatrix &A, const rmatrix &B, cimatrix &C) {
  int m = Ub(A,1)-Lb(A,1)+1;
  int n = Ub(A,2)-Lb(A,2)+1;
  imatrix re(m,n),im(m,n);

  blasmatmul(Re(A),B,re);
  blasmatmul(Im(A),B,im);

  SetRe(C,re);
  SetIm(C,im);
}

//=========================================================================

/*
//Computes B=A+B
inline void blasadd(const rvector& A, rvector& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();

  cblas_daxpy(VecLen(A), 1.0, DA, 1, DB, 1 );
}

//Computes B=A+B
inline void blasadd(const rvector& A, ivector& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
  int rnd = getround();
 
  setround(-1);
  cblas_daxpy(VecLen(A), 1.0, DA, 1, DB, 2 );
  setround(1);
  cblas_daxpy(VecLen(A), 1.0, DA, 1, DB+1, 2 );
  setround(rnd);
}

//Computes B=A+B
inline void blasadd(const rvector& A, cvector& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
 
  cblas_daxpy(VecLen(A), 1.0, DA, 1, DB, 2 );
}

//Computes B=A+B
inline void blasadd(const rvector& A, civector& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
  int rnd = getround();
 
  setround(-1);
  cblas_daxpy(VecLen(A), 1.0, DA, 1, DB, 4 );
  setround(1);
  cblas_daxpy(VecLen(A), 1.0, DA, 1, DB+1, 4 );
  setround(rnd);
}


//Computes B=A+B
inline void blasadd(const cvector& A, cvector& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();


  cblas_daxpy(VecLen(A), 1.0, DA, 2, DB, 2 );
  cblas_daxpy(VecLen(A), 1.0, DA+1, 2, DB+1, 2 );
}

//Computes B=A+B
inline void blasadd(const cvector& A, civector& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
  int rnd = getround();

  setround(-1);
  cblas_daxpy(VecLen(A), 1.0, DA, 2, DB, 4 );
  cblas_daxpy(VecLen(A), 1.0, DA, 2, DB+2, 4 );
  setround(1);
  cblas_daxpy(VecLen(A), 1.0, DA+1, 2, DB+1, 4 );
  cblas_daxpy(VecLen(A), 1.0, DA+1, 2, DB+3, 4 );
  setround(rnd);
}

//Computes B=A+B
inline void blasadd(const ivector& A, ivector& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
  int rnd = getround();

  setround(-1);

  cblas_daxpy(VecLen(A), 1.0, DA, 2, DB, 2 );

  setround(1);

  cblas_daxpy(VecLen(A), 1.0, DA+1, 2, DB+1, 2 );

  setround(rnd);
}

//Computes B=A+B
inline void blasadd(const ivector& A, civector& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
  int rnd = getround();

  setround(-1);
  cblas_daxpy(VecLen(A), 1.0, DA, 2, DB, 4 );
  cblas_daxpy(VecLen(A), 1.0, DA, 2, DB+2, 4 );

  setround(1);
  cblas_daxpy(VecLen(A), 1.0, DA+1, 2, DB+1, 4 );
  cblas_daxpy(VecLen(A), 1.0, DA+1, 2, DB+3, 4 );

  setround(rnd);
}

//Computes B=A+B
inline void blasadd(const civector& A, civector& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
  int rnd = getround();

  setround(-1);

  cblas_daxpy(VecLen(A), 1.0, DA, 4, DB, 4 );
  cblas_daxpy(VecLen(A), 1.0, DA+2, 4, DB+2, 4 );

  setround(1);

  cblas_daxpy(VecLen(A), 1.0, DA+1, 4, DB+1, 4 );
  cblas_daxpy(VecLen(A), 1.0, DA+3, 4, DB+3, 4 );

  setround(rnd);
}

//Computes C=A+B
inline void blasadd(const cvector& A, const ivector& B, civector& C) {
  C = A;
  double* DB = B.to_blas_array();
  double* DC = C.to_blas_array();
  int rnd = getround();

  setround(-1);

  cblas_daxpy(VecLen(A), 1.0, DB, 4, DC, 4 );

  setround(1);

  cblas_daxpy(VecLen(A), 1.0, DB+1, 4, DC+1, 4 );

  setround(rnd);
}

//Computes B=-A+B
inline void blassub(const rvector& A, rvector& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();

  cblas_daxpy(VecLen(A), -1.0, DA, 1, DB, 1 );
}

//Computes B=-A+B
inline void blassub(const rvector& A, ivector& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
  int rnd = getround();
 
  setround(-1);
  cblas_daxpy(VecLen(A), -1.0, DA, 1, DB, 2 );
  setround(1);
  cblas_daxpy(VecLen(A), -1.0, DA, 1, DB+1, 2 );
  setround(rnd);
}

//Computes B=-A+B
inline void blassub(const rvector& A, cvector& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
 
  cblas_daxpy(VecLen(A), -1.0, DA, 1, DB, 2 );
}

//Computes B=-A+B
inline void blassub(const rvector& A, civector& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
  int rnd = getround();
 
  setround(-1);
  cblas_daxpy(VecLen(A), -1.0, DA, 1, DB, 4 );
  setround(1);
  cblas_daxpy(VecLen(A), -1.0, DA, 1, DB+1, 4 );
  setround(rnd);
}


//Computes B=-A+B
inline void blassub(const cvector& A, cvector& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();


  cblas_daxpy(VecLen(A), -1.0, DA, 2, DB, 2 );
  cblas_daxpy(VecLen(A), -1.0, DA+1, 2, DB+1, 2 );
}

//Computes B=-A+B
inline void blassub(const cvector& A, civector& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
  int rnd = getround();

  setround(-1);
  cblas_daxpy(VecLen(A), -1.0, DA, 2, DB, 4 );
  cblas_daxpy(VecLen(A), -1.0, DA, 2, DB+2, 4 );
  setround(1);
  cblas_daxpy(VecLen(A), -1.0, DA+1, 2, DB+1, 4 );
  cblas_daxpy(VecLen(A), -1.0, DA+1, 2, DB+3, 4 );
  setround(rnd);
}

//Computes B=-A+B
inline void blassub(const ivector& A, ivector& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
  int rnd = getround();

  setround(-1);

  cblas_daxpy(VecLen(A), -1.0, DA+1, 2, DB, 2 );

  setround(1);

  cblas_daxpy(VecLen(A), -1.0, DA, 2, DB+1, 2 );

  setround(rnd);
}

//Computes B=-A+B
inline void blassub(const ivector& A, civector& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
  int rnd = getround();

  setround(-1);
  cblas_daxpy(VecLen(A), -1.0, DA+1, 2, DB, 4 );
  cblas_daxpy(VecLen(A), -1.0, DA+1, 2, DB+2, 4 );

  setround(1);
  cblas_daxpy(VecLen(A), -1.0, DA, 2, DB+1, 4 );
  cblas_daxpy(VecLen(A), -1.0, DA, 2, DB+3, 4 );

  setround(rnd);
}

//Computes B=-A+B
inline void blassub(const civector& A, civector& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
  int rnd = getround();

  setround(-1);

  cblas_daxpy(VecLen(A), -1.0, DA+1, 4, DB, 4 );
  cblas_daxpy(VecLen(A), -1.0, DA+3, 4, DB+2, 4 );

  setround(1);

  cblas_daxpy(VecLen(A), -1.0, DA, 4, DB+1, 4 );
  cblas_daxpy(VecLen(A), -1.0, DA+2, 4, DB+3, 4 );

  setround(rnd);
}

//==============================================================

//Computes B=A+B
inline void blasadd(const rmatrix& A, rmatrix& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();

  cblas_daxpy(RowLen(A)*ColLen(A), 1.0, DA, 1, DB, 1 );
}

//Computes B=A+B
inline void blasadd(const rmatrix& A, imatrix& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
  int rnd = getround();
 
  setround(-1);
  cblas_daxpy(RowLen(A)*ColLen(A), 1.0, DA, 1, DB, 2 );
  setround(1);
  cblas_daxpy(RowLen(A)*ColLen(A), 1.0, DA, 1, DB+1, 2 );
  setround(rnd);
}

//Computes B=A+B
inline void blasadd(const rmatrix& A, cmatrix& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
 
  cblas_daxpy(RowLen(A)*ColLen(A), 1.0, DA, 1, DB, 2 );
}

//Computes B=A+B
inline void blasadd(const rmatrix& A, cimatrix& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
  int rnd = getround();
 
  setround(-1);
  cblas_daxpy(RowLen(A)*ColLen(A), 1.0, DA, 1, DB, 4 );
  setround(1);
  cblas_daxpy(RowLen(A)*ColLen(A), 1.0, DA, 1, DB+1, 4 );
  setround(rnd);
}


//Computes B=A+B
inline void blasadd(const cmatrix& A, cmatrix& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();


  cblas_daxpy(RowLen(A)*ColLen(A), 1.0, DA, 2, DB, 2 );
  cblas_daxpy(RowLen(A)*ColLen(A), 1.0, DA+1, 2, DB+1, 2 );
}

//Computes B=A+B
inline void blasadd(const cmatrix& A, cimatrix& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
  int rnd = getround();

  setround(-1);
  cblas_daxpy(RowLen(A)*ColLen(A), 1.0, DA, 2, DB, 4 );
  cblas_daxpy(RowLen(A)*ColLen(A), 1.0, DA, 2, DB+2, 4 );
  setround(1);
  cblas_daxpy(RowLen(A)*ColLen(A), 1.0, DA+1, 2, DB+1, 4 );
  cblas_daxpy(RowLen(A)*ColLen(A), 1.0, DA+1, 2, DB+3, 4 );
  setround(rnd);
}

//Computes B=A+B
inline void blasadd(const imatrix& A, imatrix& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
  int rnd = getround();

  setround(-1);

  cblas_daxpy(RowLen(A)*ColLen(A), 1.0, DA, 2, DB, 2 );

  setround(1);

  cblas_daxpy(RowLen(A)*ColLen(A), 1.0, DA+1, 2, DB+1, 2 );

  setround(rnd);
}

//Computes B=A+B
inline void blasadd(const imatrix& A, cimatrix& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
  int rnd = getround();

  setround(-1);
  cblas_daxpy(RowLen(A)*ColLen(A), 1.0, DA, 2, DB, 4 );
  cblas_daxpy(RowLen(A)*ColLen(A), 1.0, DA, 2, DB+2, 4 );

  setround(1);
  cblas_daxpy(RowLen(A)*ColLen(A), 1.0, DA+1, 2, DB+1, 4 );
  cblas_daxpy(RowLen(A)*ColLen(A), 1.0, DA+1, 2, DB+3, 4 );

  setround(rnd);
}

//Computes B=A+B
inline void blasadd(const cimatrix& A, cimatrix& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
  int rnd = getround();

  setround(-1);

  cblas_daxpy(RowLen(A)*ColLen(A), 1.0, DA, 4, DB, 4 );
  cblas_daxpy(RowLen(A)*ColLen(A), 1.0, DA+2, 4, DB+2, 4 );

  setround(1);

  cblas_daxpy(RowLen(A)*ColLen(A), 1.0, DA+1, 4, DB+1, 4 );
  cblas_daxpy(RowLen(A)*ColLen(A), 1.0, DA+3, 4, DB+3, 4 );

  setround(rnd);
}

//Computes C=A+B
inline void blasadd(const cmatrix& A, const imatrix& B, cimatrix& C) {
  C = A;
  double* DB = B.to_blas_array();
  double* DC = C.to_blas_array();
  int rnd = getround();

  setround(-1);

  cblas_daxpy(RowLen(A)*ColLen(A), 1.0, DB, 4, DC, 4 );

  setround(1);

  cblas_daxpy(RowLen(A)*ColLen(A), 1.0, DB+1, 4, DC+1, 4 );

  setround(rnd);
}

//Computes B=-A+B
inline void blassub(const rmatrix& A, rmatrix& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();

  cblas_daxpy(RowLen(A)*ColLen(A), -1.0, DA, 1, DB, 1 );
}

//Computes B=-A+B
inline void blassub(const rmatrix& A, imatrix& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
  int rnd = getround();
 
  setround(-1);
  cblas_daxpy(RowLen(A)*ColLen(A), -1.0, DA, 1, DB, 2 );
  setround(1);
  cblas_daxpy(RowLen(A)*ColLen(A), -1.0, DA, 1, DB+1, 2 );
  setround(rnd);
}

//Computes B=-A+B
inline void blassub(const rmatrix& A, cmatrix& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
 
  cblas_daxpy(RowLen(A)*ColLen(A), -1.0, DA, 1, DB, 2 );
}

//Computes B=-A+B
inline void blassub(const rmatrix& A, cimatrix& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
  int rnd = getround();
 
  setround(-1);
  cblas_daxpy(RowLen(A)*ColLen(A), -1.0, DA, 1, DB, 4 );
  setround(1);
  cblas_daxpy(RowLen(A)*ColLen(A), -1.0, DA, 1, DB+1, 4 );
  setround(rnd);
}


//Computes B=-A+B
inline void blassub(const cmatrix& A, cmatrix& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();


  cblas_daxpy(RowLen(A)*ColLen(A), -1.0, DA, 2, DB, 2 );
  cblas_daxpy(RowLen(A)*ColLen(A), -1.0, DA+1, 2, DB+1, 2 );
}

//Computes B=-A+B
inline void blassub(const cmatrix& A, cimatrix& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
  int rnd = getround();

  setround(-1);
  cblas_daxpy(RowLen(A)*ColLen(A), -1.0, DA, 2, DB, 4 );
  cblas_daxpy(RowLen(A)*ColLen(A), -1.0, DA, 2, DB+2, 4 );
  setround(1);
  cblas_daxpy(RowLen(A)*ColLen(A), -1.0, DA+1, 2, DB+1, 4 );
  cblas_daxpy(RowLen(A)*ColLen(A), -1.0, DA+1, 2, DB+3, 4 );
  setround(rnd);
}

//Computes B=-A+B
inline void blassub(const imatrix& A, imatrix& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
  int rnd = getround();

  setround(-1);

  cblas_daxpy(RowLen(A)*ColLen(A), -1.0, DA+1, 2, DB, 2 );

  setround(1);

  cblas_daxpy(RowLen(A)*ColLen(A), -1.0, DA, 2, DB+1, 2 );

  setround(rnd);
}

//Computes B=-A+B
inline void blassub(const imatrix& A, cimatrix& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
  int rnd = getround();

  setround(-1);
  cblas_daxpy(RowLen(A)*ColLen(A), -1.0, DA+1, 2, DB, 4 );
  cblas_daxpy(RowLen(A)*ColLen(A), -1.0, DA+1, 2, DB+2, 4 );

  setround(1);
  cblas_daxpy(RowLen(A)*ColLen(A), -1.0, DA, 2, DB+1, 4 );
  cblas_daxpy(RowLen(A)*ColLen(A), -1.0, DA, 2, DB+3, 4 );

  setround(rnd);
}

//Computes B=-A+B
inline void blassub(const cimatrix& A, cimatrix& B) {
  double* DA = A.to_blas_array();
  double* DB = B.to_blas_array();
  int rnd = getround();

  setround(-1);

  cblas_daxpy(RowLen(A)*ColLen(A), -1.0, DA+1, 4, DB, 4 );
  cblas_daxpy(RowLen(A)*ColLen(A), -1.0, DA+3, 4, DB+2, 4 );

  setround(1);

  cblas_daxpy(RowLen(A)*ColLen(A), -1.0, DA, 4, DB+1, 4 );
  cblas_daxpy(RowLen(A)*ColLen(A), -1.0, DA+2, 4, DB+3, 4 );

  setround(rnd);
} */

#endif
#endif

}
