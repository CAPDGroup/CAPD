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

/* CVS $Id: sparseidot.inl,v 1.12 2014/01/30 17:23:49 cxsc Exp $ */

#include "sparseidot.hpp"
#include "idotk.inl"

namespace cxsc {
        //Rounding mode must be set to upwards!
	inline void sort(const interval &x, const interval &y, real& x_inf, real& y_inf, real &x_sup, real &y_sup) {
		
		if(Inf(x) >= 0  &&  Sup(x) >= 0) {
	
			if(Inf(y) >= 0  &&  Sup(y) >= 0) {
				x_inf = Inf(x);
				y_inf = Inf(y);
				x_sup = Sup(x);
				y_sup = Sup(y);
			} else if(Inf(y) < 0  &&  Sup(y) >= 0) {
				x_inf = Sup(x);
				y_inf = Inf(y);
				x_sup = Sup(x);
				y_sup = Sup(y);
			} else {
				x_inf = Sup(x);
				y_inf = Inf(y);
				x_sup = Inf(x);
				y_sup = Sup(y);
			}
	
		} else if(Inf(x) < 0  &&  Sup(x) >= 0) {
	
			if(Inf(y) >= 0  &&  Sup(y) >= 0) {
				x_inf = Inf(x);
				y_inf = Sup(y);
				x_sup = Sup(x);
				y_sup = Sup(y);
			} else if(Inf(y) < 0  &&  Sup(y) >= 0) {

#ifndef _CXSC_DOTK_ROUND_SOFT
				if((-Inf(x))*Sup(y) >= Sup(x)*(-Inf(y))) {
#else
				if( mulu(-Inf(x),Sup(y)) >= mulu(Sup(x),-Inf(y)) ) {
#endif
					x_inf = Inf(x);
					y_inf = Sup(y);
				} else {
					x_inf = Sup(x);
					y_inf = Inf(y);
				}

#ifndef _CXSC_DOTK_ROUND_SOFT
				if(Inf(x)*Inf(y) >= Sup(x)*Sup(y)) {
#else
				if( mulu(Inf(x),Inf(y)) >= mulu(Sup(x),Sup(y)) ) {
#endif
					x_sup = Inf(x);
					y_sup = Inf(y);
				} else {
					x_sup = Sup(x);
					y_sup = Sup(y);
				}
				
			} else {
				x_inf = Sup(x);
				y_inf = Inf(y);
				x_sup = Inf(x);
				y_sup = Inf(y);
			}
	
		} else {
	
			if(Inf(y) >= 0  &&  Sup(y) >= 0) {
				x_inf = Inf(x);
				y_inf = Sup(y);
				x_sup = Sup(x);
				y_sup = Inf(y);
			} else if(Inf(y) < 0  &&  Sup(y) >= 0) {
				x_inf = Inf(x);
				y_inf = Sup(y);
				x_sup = Inf(x);
				y_sup = Inf(y);
			} else {
				x_inf = Sup(x);
				y_inf = Sup(y);
				x_sup = Inf(x);
				y_sup = Inf(y);
			}
	
		}
	
		
	}

	//Sorts inf and sup of an interval vector and a real vector into two real vector for the computation of the infimum 
	//and two real vectors for the computation of the supremum. 
	inline void sort(const interval &x, const real &y, real& x_inf, real& y_inf, real &x_sup, real &y_sup) {
		
			if(Inf(x) >= 0  &&  Sup(x) >= 0) {
		
				if(y >= 0) {
					x_inf = Inf(x);
					y_inf = y;
					x_sup = Sup(x);
					y_sup = y;
				} else {
					x_inf = Sup(x);
					y_inf = y;
					x_sup = Inf(x);
					y_sup = y;
				}
		
			} else if(Inf(x) < 0  &&  Sup(x) >= 0) {
		
				if(y >= 0) {
					x_inf = Inf(x);
					y_inf = y;
					x_sup = Sup(x);
					y_sup = y;
				} else {
					x_inf = Sup(x);
					y_inf = y;
					x_sup = Inf(x);
					y_sup = y;
				}
		
			} else {
		
				if(y >= 0) {
					x_inf = Inf(x);
					y_inf = y;
					x_sup = Sup(x);
					y_sup = y;
				} else {
					x_inf = Sup(x);
					y_inf = y;
					x_sup = Inf(x);
					y_sup = y;
				}
		
			}
		
		
	}

	//Sorts inf and sup of an interval vector and a real vector into two real vector for the computation of the infimum 
	//and two real vectors for the computation of the supremum. Rounding mode must be set to upwards!
	inline void sort(const real &y, const interval &x, real& x_inf, real& y_inf, real &x_sup, real &y_sup) {
		
			if(Inf(x) >= 0  &&  Sup(x) >= 0) {
		
				if(y >= 0) {
					x_inf = Inf(x);
					y_inf = y;
					x_sup = Sup(x);
					y_sup = y;
				} else {
					x_inf = Sup(x);
					y_inf = y;
					x_sup = Inf(x);
					y_sup = y;
				}
		
			} else if(Inf(x) < 0  &&  Sup(x) >= 0) {
		
				if(y >= 0) {
					x_inf = Inf(x);
					y_inf = y;
					x_sup = Sup(x);
					y_sup = y;
				} else {
					x_inf = Sup(x);
					y_inf = y;
					x_sup = Inf(x);
					y_sup = y;
				}
		
			} else {
		
				if(y >= 0) {
					x_inf = Inf(x);
					y_inf = y;
					x_sup = Sup(x);
					y_sup = y;
				} else {
					x_inf = Sup(x);
					y_inf = y;
					x_sup = Inf(x);
					y_sup = y;
				}
		
			}
		
		
	}


    sparse_idot::sparse_idot(unsigned int p) : val(0.0),corr_inf(0.0), corr_sup(0.0), err_inf(0.0),err_sup(0.0),k(p),n(0) {
      if(k==0) dot = new idotprecision(0.0);
      else dot = NULL;
    }

    sparse_idot::sparse_idot(unsigned int p, int nnz) : val(0.0),corr_inf(0.0), corr_sup(0.0), err_inf(0.0),err_sup(0.0),k(p),n(0) {
      if(k==0) dot = new idotprecision(0.0);
      else dot = NULL;
      cm_x.reserve(nnz);
      cm_y.reserve(nnz);
      ca_x.reserve(nnz);
      ca_y.reserve(nnz);
    }

    sparse_idot::sparse_idot(const sparse_idot& s) : cm_x(s.cm_x), cm_y(s.cm_y), ca_x(s.ca_x), ca_y(s.ca_y), val(s.val), corr_inf(s.corr_inf), corr_sup(s.corr_sup), err_inf(s.err_inf), err_sup(s.err_sup), k(s.k), n(s.n) {
      if(this != &s) {
        if(s.dot == NULL) {
          dot = NULL;
        } else {
          dot = new idotprecision();
          *dot = *(s.dot);
        }
      }
    }

    sparse_idot::~sparse_idot() {
      if(dot != NULL) delete dot;
    }

    void sparse_idot::reset() {
      if(k==0) {
        *dot = 0.0;
      } else if(k==1) {
        val = 0.0;
        cm_x.clear();
        cm_y.clear();
        ca_x.clear();
        ca_y.clear();
      } else {
        cm_x.clear();
        cm_y.clear();
        ca_x.clear();
        ca_y.clear();
        val = 0.0;
        corr_inf = corr_sup = err_inf = err_sup = 0.0;
      }
      n = 0;
    }

    void sparse_idot::add_dot(const interval& x, const interval& y) {
      if(k==0) {

        accumulate(*dot, x, y);

      } else if(k==1) {

        real x_inf,y_inf,x_sup,y_sup;
#ifndef _CXSC_DOTK_ROUND_SOFT
	setround(1);
        sort(x,y,x_inf,y_inf,x_sup,y_sup);
        setround(0);
#else
        sort(x,y,x_inf,y_inf,x_sup,y_sup);
#endif
        cm_x.push_back(x_inf);
        ca_x.push_back(x_sup);
        cm_y.push_back(y_inf);
        ca_y.push_back(y_sup);

      } else if(k==2) {
     
        interval a;
        real b_inf, b_sup, c_inf, c_sup, t;
        TwoProduct(x,y,a,b_inf,b_sup);
        TwoSum(val,a,val,c_inf,c_sup);
        t = b_inf+c_inf;
        err_inf += abs(t);
        corr_inf += t;
        t = b_sup+c_sup;
        err_sup += abs(t);
        corr_sup += t;

      } else if(k>=3) {

        interval a;
        real b_inf,b_sup;
        TwoProduct(x,y,a,b_inf,b_sup);
        cm_x.push_back(b_inf);
        cm_y.push_back(b_sup);
        TwoSum(val,a,val,b_inf,b_sup);
        ca_x.push_back(b_inf);
        ca_y.push_back(b_sup);

      }

      n++;
    }

    void sparse_idot::add_dot(const interval& x, const real& y) {
      if(k==0) {

        accumulate(*dot, x, y);

      } else if(k==1) {

        real x_inf,y_inf,x_sup,y_sup;
        sort(x,y,x_inf,y_inf,x_sup,y_sup);
        cm_x.push_back(x_inf);
        ca_x.push_back(x_sup);
        cm_y.push_back(y_inf);
        ca_y.push_back(y_sup);

      } else if(k==2) {
     
        interval a;
        real b_inf, b_sup, c_inf, c_sup, t;
        TwoProduct(x,y,a,b_inf,b_sup);
        TwoSum(val,a,val,c_inf,c_sup);
        t = b_inf+c_inf;
        err_inf += abs(t);
        corr_inf += t;
        t = b_sup+c_sup;
        err_sup += abs(t);
        corr_sup += t;

      } else if(k>=3) {

        interval a;
        real b_inf,b_sup;
        TwoProduct(x,y,a,b_inf,b_sup);
        cm_x.push_back(b_inf);
        cm_y.push_back(b_sup);
        TwoSum(val,a,val,b_inf,b_sup);
        ca_x.push_back(b_inf);
        ca_y.push_back(b_sup);

      }

      n++;
    }

    void sparse_idot::add_dot(const real& x, const interval& y) {
      if(k==0) {

        accumulate(*dot, x, y);

      } else if(k==1) {

        real x_inf,y_inf,x_sup,y_sup;
        sort(x,y,x_inf,y_inf,x_sup,y_sup);
        cm_x.push_back(x_inf);
        ca_x.push_back(x_sup);
        cm_y.push_back(y_inf);
        ca_y.push_back(y_sup);

      } else if(k==2) {
     
        interval a;
        real b_inf, b_sup, c_inf, c_sup, t;
        TwoProduct(x,y,a,b_inf,b_sup);
        TwoSum(val,a,val,c_inf,c_sup);
        t = b_inf+c_inf;
        err_inf += abs(t);
        corr_inf += t;
        t = b_sup+c_sup;
        err_sup += abs(t);
        corr_sup += t;

      } else if(k>=3) {

        interval a;
        real b_inf,b_sup;
        TwoProduct(x,y,a,b_inf,b_sup);
        cm_x.push_back(b_inf);
        cm_y.push_back(b_sup);
        TwoSum(val,a,val,b_inf,b_sup);
        ca_x.push_back(b_inf);
        ca_y.push_back(b_sup);

      }

      n++;
    }

    void sparse_idot::add_dot_err(const interval& x, const interval& y) {
      add_dot(x,y);
    }

    void sparse_idot::add_dot_err(const interval& x, const real& y) {
      add_dot(x,y);
    }

    void sparse_idot::add_dot_err(const real& x, const interval& y) {
      add_dot(x,y);
    }

    interval sparse_idot::result() {
   
      if(n==0) return interval(0.0,0.0);

      if(k==0) {

        return rnd(*dot);

      } else if(k==1) {

       
#ifndef _CXSC_DOTK_ROUND_SOFT
        setround(-1);
        for(unsigned int i=0 ; i<cm_x.size() ; i++) {
          Inf(val) += cm_x[i]*cm_y[i];
        }

        setround(1);
        for(unsigned int i=0 ; i<ca_x.size() ; i++) {
          Sup(val) += ca_x[i]*ca_y[i];
        }

	setround(0);
#else
        for(unsigned int i=0 ; i<cm_x.size() ; i++) {
          Inf(val) = addd(Inf(val), muld(cm_x[i],cm_y[i]));
        }

        for(unsigned int i=0 ; i<ca_x.size() ; i++) {
          Sup(val) = addu(Sup(val), mulu(ca_x[i],ca_y[i]));
        }
#endif

        return val;

      } else if(k==2) {

	real alpha, delta, error_inf, error_sup;

        Inf(val) += corr_inf;
        Sup(val) += corr_sup;

	//Error infimum	
	delta = (n*Epsilon) / (1.0-2*n*Epsilon);
	alpha = (Epsilon*abs(Inf(val))) + (delta*err_inf+3*MinReal/Epsilon);
	error_inf = alpha / (1.0 - 2*Epsilon);
	
	//Error supremum
	alpha = (Epsilon*abs(Sup(val))) + (delta*err_sup+3*MinReal/Epsilon);
	error_sup = alpha / (1.0 - 2*Epsilon);

        return val + interval(-error_inf,error_sup);

      } else if(k>=3) {

	for(int j=1 ; j<k-1 ; j++) {
		for(int i=1 ; i<n ; i++) 
			TwoSum(cm_x[i],cm_x[i-1],cm_x[i],cm_x[i-1]);    
              	TwoSum(ca_x[0],cm_x[n-1],ca_x[0],cm_x[n-1]);    
		for(int i=1 ; i<n ; i++) 
			TwoSum(ca_x[i],ca_x[i-1],ca_x[i],ca_x[i-1]);    
		TwoSum(Inf(val),ca_x[n-1],Inf(val),ca_x[n-1]);
	}
		
	corr_inf = std::accumulate(cm_x.begin(),cm_x.end(),corr_inf);
	corr_inf = std::accumulate(ca_x.begin(),ca_x.end(),corr_inf);

        real tmperr = 0.0;

        for(unsigned int j=0 ; j<cm_x.size() ; j++) {
          tmperr += abs(cm_x[j]);
        }

        for(unsigned int j=0 ; j<ca_x.size() ; j++) {
          tmperr += abs(ca_x[j]);
        }
		
	real alpha, delta, error;
	
	delta = (n*Epsilon) / (1.0-2*n*Epsilon);
	alpha = (Epsilon*abs(Inf(val)+corr_inf)) + (delta*tmperr+3*MinReal/Epsilon);
	error = alpha / (1.0 - 2*Epsilon);

        err_inf = error;

	Inf(val) = subd(addd(Inf(val),corr_inf), err_inf);


	for(int j=1 ; j<k-1 ; j++) {
		for(int i=1 ; i<n ; i++) 
			TwoSum(cm_y[i],cm_y[i-1],cm_y[i],cm_y[i-1]);    
              	TwoSum(ca_y[0],cm_y[n-1],ca_y[0],cm_y[n-1]);    
		for(int i=1 ; i<n ; i++) 
			TwoSum(ca_y[i],ca_y[i-1],ca_y[i],ca_y[i-1]);    
		TwoSum(Sup(val),ca_y[n-1],Sup(val),ca_y[n-1]);
	}
		
	corr_sup = std::accumulate(cm_y.begin(),cm_y.end(),corr_sup);
	corr_sup = std::accumulate(ca_y.begin(),ca_y.end(),corr_sup);

        tmperr = 0.0;

        for(unsigned int j=0 ; j<cm_y.size() ; j++) {
          tmperr += abs(cm_y[j]);
        }

        for(unsigned int j=0 ; j<ca_y.size() ; j++) {
          tmperr += abs(ca_y[j]);
        }
		
	delta = (n*Epsilon) / (1.0-2*n*Epsilon);
	alpha = (Epsilon*abs(Sup(val)+corr_sup)) + (delta*tmperr+3*MinReal/Epsilon);
	error = alpha / (1.0 - 2*Epsilon);

        err_sup = error;

        Sup(val) = addu(Sup(val), addu(corr_sup,err_sup));		

        return val;

      }
      return val;
    }

    void sparse_idot::result(idotprecision& res_dot) {
   
      if(n==0) return;

      if(k==0) {

        res_dot += *dot;

      } else if(k==1) {
       
#ifndef _CXSC_DOTK_ROUND_SOFT
        setround(-1);
        for(unsigned int i=0 ; i<cm_x.size() ; i++) {
          Inf(val) += cm_x[i]*cm_y[i];
        }

        setround(1);
        for(unsigned int i=0 ; i<ca_x.size() ; i++) {
          Sup(val) += ca_x[i]*ca_y[i];
        }

	setround(0);
#else
        for(unsigned int i=0 ; i<cm_x.size() ; i++) {
          Inf(val) = addd(Inf(val), muld(cm_x[i],cm_y[i]));
        }

        for(unsigned int i=0 ; i<ca_x.size() ; i++) {
          Sup(val) = addu(Sup(val), mulu(ca_x[i],ca_y[i]));
        }
#endif
        res_dot += val;

      } else if(k==2) {

	real alpha, delta, error_inf(0.0), error_sup(0.0);
	delta = (n*Epsilon) / (1.0-2*n*Epsilon);

	//Error infimum	
	alpha = (Epsilon*abs(Inf(val))) + (delta*err_inf+3*MinReal/Epsilon);
	error_inf = alpha / (1.0 - 2*Epsilon);
	
	//Error supremum
	alpha = (Epsilon*abs(Sup(val))) + (delta*err_sup+3*MinReal/Epsilon);
	error_sup = alpha / (1.0 - 2*Epsilon);

        res_dot += val + interval(-error_inf,error_sup);
        Inf(res_dot) += corr_inf;
        Sup(res_dot) += corr_sup;

      } else if(k>=3) {

	for(int j=1 ; j<k-1 ; j++) {
		for(int i=1 ; i<n ; i++) 
			TwoSum(cm_x[i],cm_x[i-1],cm_x[i],cm_x[i-1]);    
              	TwoSum(ca_x[0],cm_x[n-1],ca_x[0],cm_x[n-1]);    
		for(int i=1 ; i<n ; i++) 
			TwoSum(ca_x[i],ca_x[i-1],ca_x[i],ca_x[i-1]);    
		TwoSum(Inf(val),ca_x[n-1],Inf(val),ca_x[n-1]);
	}
		
	corr_inf = std::accumulate(cm_x.begin(),cm_x.end(),corr_inf);
	corr_inf = std::accumulate(ca_x.begin(),ca_x.end(),corr_inf);

        real tmperr = 0.0;

        for(unsigned int j=0 ; j<cm_x.size() ; j++) {
          tmperr += abs(cm_x[j]);
        }

        for(unsigned int j=0 ; j<ca_x.size() ; j++) {
          tmperr += abs(ca_x[j]);
        }
		
	real alpha, delta, error;
	
	delta = (n*Epsilon) / (1.0-2*n*Epsilon);
	alpha = (Epsilon*abs(Inf(val)+corr_inf)) + (delta*tmperr+3*MinReal/Epsilon);
	error = alpha / (1.0 - 2*Epsilon);

        err_inf = error;

	Inf(res_dot) += corr_inf;
	Inf(res_dot) -= err_inf;

	for(int j=1 ; j<k-1 ; j++) {
		for(int i=1 ; i<n ; i++) 
			TwoSum(cm_y[i],cm_y[i-1],cm_y[i],cm_y[i-1]);    
              	TwoSum(ca_y[0],cm_y[n-1],ca_y[0],cm_y[n-1]);    
		for(int i=1 ; i<n ; i++) 
			TwoSum(ca_y[i],ca_y[i-1],ca_y[i],ca_y[i-1]);    
		TwoSum(Sup(val),ca_y[n-1],Sup(val),ca_y[n-1]);
	}


	corr_sup = std::accumulate(cm_y.begin(),cm_y.end(),corr_sup);
	corr_sup = std::accumulate(ca_y.begin(),ca_y.end(),corr_sup);

        tmperr = 0.0;

        for(unsigned int j=0 ; j<cm_y.size() ; j++) {
          tmperr += abs(cm_y[j]);
        }

        for(unsigned int j=0 ; j<ca_y.size() ; j++) {
          tmperr += abs(ca_y[j]);
        }
		
	delta = (n*Epsilon) / (1.0-2*n*Epsilon);
	alpha = (Epsilon*abs(Sup(val)+corr_sup)) + (delta*tmperr+3*MinReal/Epsilon);
	error = alpha / (1.0 - 2*Epsilon);

        err_sup = error;

        Sup(res_dot) += corr_sup;
        Sup(res_dot) += err_sup;

      }
    }


}
