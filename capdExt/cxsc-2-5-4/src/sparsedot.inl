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

/* CVS $Id: sparsedot.inl,v 1.13 2014/01/30 17:23:49 cxsc Exp $ */

#include <real.hpp>
#include <dot.hpp>
#include <vector>
#include <numeric>
#include <sparsedot.hpp>
#include <dotk.inl>

namespace cxsc {

    sparse_dot::sparse_dot(unsigned int p) : val(0.0),corr(0.0),err(0.0),n(0),k(p) {
      if(k==0) dot = new dotprecision(0.0);
      else dot = NULL;
    }

    sparse_dot::sparse_dot(const sparse_dot& s) : cm(s.cm), ca(s.ca), val(s.val), corr(s.corr), err(s.err), n(s.n), k(s.k) {
      if(this != &s) {
        if(s.dot == NULL) {
          dot = NULL;
        } else {
          dot = new dotprecision();
          *dot = *(s.dot);
        }
      }
    }

     sparse_dot::~sparse_dot() {
      if(dot != NULL) delete dot;
    }

     void sparse_dot::reset() {
      if(k==0) {
        *dot = 0.0;
      } else if(k==1) {
        val = err = 0.0;
      } else {
        cm.clear();
        ca.clear();
        val = corr = err = 0.0;
      }
      n = 0;
    }

     void sparse_dot::add_dot(const real& x, const real& y) {
      if(k==0) {
        accumulate(*dot, x, y);
      } else if(k==1) {
        val += x*y;
      } else if(k==2) {
        real a,b,c;
        TwoProduct(x,y,a,b);
        TwoSum(val,a,val,c);
        corr += (b+c);
      } else if(k>=3) {
        real a,b;
        TwoProduct(x,y,a,b);
        cm.push_back(b);
        TwoSum(val,a,val,b);
        ca.push_back(b);
      }
    }

     void sparse_dot::add_dot_err(const real& x, const real& y) {
      if(k==0) {
        accumulate(*dot, x, y);
      } else if(k==1) {
        ca.push_back(x);
        cm.push_back(y);
      } else if(k==2) {
        real a,b,c;
        TwoProduct(x,y,a,b);
        TwoSum(val,a,val,c);
        c += b;
        err += abs(c);
        corr += c;
        n++;
      } else if(k>=3) {
        real a,b;
        TwoProduct(x,y,a,b);
        cm.push_back(b);
        TwoSum(val,a,val,b);
        ca.push_back(b);
      }
    }

     real sparse_dot::result() {
      if(k==0) {
        return rnd(*dot);
      } else if(k==1) {
        return val;
      } else if(k==2) {
        return val+corr;
      } else if(k>=3) {
        n = cm.size();
        if(n == 0) return val;

	for(int j=1 ; j<k-1 ; j++) {
		for(int i=1 ; i<n ; i++) 
			TwoSum(cm[i],cm[i-1],cm[i],cm[i-1]);    
              	TwoSum(ca[0],cm[n-1],ca[0],cm[n-1]);    
		for(int i=1 ; i<n ; i++) 
			TwoSum(ca[i],ca[i-1],ca[i],ca[i-1]);    
		TwoSum(val,ca[n-1],val,ca[n-1]);
	}
		
	corr = std::accumulate(cm.begin(),cm.end(),corr);
	corr = std::accumulate(ca.begin(),ca.end(),corr);
		
	val += corr;

        return val;
      }
      return val;
    }

     void sparse_dot::result(dotprecision& res_dot) {
      if(k==0) {
        res_dot += *dot;
      } else if(k==1) {
        real resd = 0.0, resu = 0.0;
        
#ifndef _CXSC_DOTK_ROUND_SOFT
        setround(-1);
        for(unsigned int i=0 ; i<ca.size() ; i++)
          resd += ca[i]*cm[i];

        setround(1);
        for(unsigned int i=0 ; i<ca.size() ; i++)
          resu += ca[i]*cm[i];

        setround(0);
        val = resd+(resu-resd)*0.5;

        setround(1);
        res_dot.set_err(res_dot.get_err() + (resu-val));

        setround(0);
#else
        for(unsigned int i=0 ; i<ca.size() ; i++)
          resd = addd(resd, muld(ca[i],cm[i]));

        for(unsigned int i=0 ; i<ca.size() ; i++)
          resu = addu(resu, mulu(ca[i],cm[i]));

        val = resd+(resu-resd)*0.5;

        res_dot.set_err(addu(res_dot.get_err(), subu(resu,val)));
#endif

        res_dot += val;

      } else if(k==2) {

	real alpha, delta, error;
	delta = (n*Epsilon) / (1.0-2*n*Epsilon);
	alpha = (Epsilon*abs(val)) + (delta*err+3*MinReal/Epsilon);
	error = alpha / (1.0 - 2*Epsilon);
	res_dot.set_err(addu(res_dot.get_err(), error));      

        res_dot += val;
        res_dot += corr;

      } else if(k>=3) {

        n = cm.size();
        if(n==0) return;
 
	for(int j=1 ; j<k-1 ; j++) {
		for(int i=1 ; i<n ; i++) 
			TwoSum(cm[i],cm[i-1],cm[i],cm[i-1]);    
              	TwoSum(ca[0],cm[n-1],ca[0],cm[n-1]);    
		for(int i=1 ; i<n ; i++) 
			TwoSum(ca[i],ca[i-1],ca[i],ca[i-1]);    
		TwoSum(val,ca[n-1],val,ca[n-1]);
	}
		
	corr = std::accumulate(cm.begin(),cm.end(),corr);
	corr = std::accumulate(ca.begin(),ca.end(),corr);

        res_dot += val;
        res_dot += corr;

        real tmperr = 0.0;

        for(unsigned int j=0 ; j<cm.size() ; j++) {
          tmperr += abs(cm[j]);
        }

        for(unsigned int j=0 ; j<ca.size() ; j++) {
          tmperr += abs(ca[j]);
        }
		
	real alpha, delta, error;
	
	delta = (n*Epsilon) / (1.0-2*n*Epsilon);
	alpha = (Epsilon*abs(val+corr)) + (delta*tmperr+3*MinReal/Epsilon);
	error = alpha / (1.0 - 2*Epsilon);

        res_dot.set_err(addu(res_dot.get_err(),error));

      }
    }



}


 
