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

/* CVS $Id: sparsecdot.inl,v 1.12 2014/01/30 17:23:49 cxsc Exp $ */

#include <sparsecdot.hpp>
#include <cdotk.inl>

namespace cxsc {

    sparse_cdot::sparse_cdot(unsigned int p) : val(0.0),corr(0.0),err(0.0),n(0),k(p) {
      if(k==0) dot = new cdotprecision(0.0);
      else dot = NULL;
    }

    sparse_cdot::sparse_cdot(const sparse_cdot& s) : cm(s.cm), ca(s.ca), val(s.val), corr(s.corr), err(s.err), n(s.n), k(s.k) {
      if(this != &s) {
        if(s.dot == NULL) {
          dot = NULL;
        } else {
          dot = new cdotprecision();
          *dot = *(s.dot);
        }
      }
    }

    sparse_cdot::~sparse_cdot() {
      if(dot != NULL) delete dot;
    }

    void sparse_cdot::reset() {
      if(k==0) {
        *dot = 0.0;
      } else if(k==1) {
        val = 0.0;
        err = 0.0;   
      } else {
        cm.clear();
        ca.clear();
        val = 0.0;
        corr = 0.0;
        err = 0.0;   
      }
      n = 0;
    }

    void sparse_cdot::add_dot(const complex& x, const complex& y) {
      if(k==0) {
        accumulate(*dot, x, y);
      } else if(k==1) {
        val += complex(Re(x)*Re(y)-Im(x)*Im(y), Re(x)*Im(y)+Im(x)*Re(y));
      } else if(k==2) {
        real a,b,c;
        TwoProduct(Re(x),Re(y),a,b);
        TwoSum(Re(val),a,Re(val),c);
        Re(corr) += (b+c);
        TwoProduct(-Im(x),Im(y),a,b);
        TwoSum(Re(val),a,Re(val),c);
        Re(corr) += (b+c);
        TwoProduct(Re(x),Im(y),a,b);
        TwoSum(Im(val),a,Im(val),c);
        Im(corr) += (b+c);
        TwoProduct(Im(x),Re(y),a,b);
        TwoSum(Im(val),a,Im(val),c);
        Im(corr) += (b+c);
      } else if(k>=3) {
        real a;
        complex b,c;
        TwoProduct(Re(x),Re(y),a,Re(b));
        TwoSum(Re(val),a,Re(val),Re(c));
        TwoProduct(Re(x),Im(y),a,Im(b));
        TwoSum(Im(val),a,Im(val),Im(c));
        cm.push_back(b);
        ca.push_back(c);

        TwoProduct(-Im(x),Im(y),a,Re(b));
        TwoSum(Re(val),a,Re(val),Re(c));
        TwoProduct(Im(x),Re(y),a,Im(b));
        TwoSum(Im(val),a,Im(val),Im(c));
        cm.push_back(b);
        ca.push_back(c);
      }
    }

    void sparse_cdot::add_dot(const complex& x, const real& y) {
      if(k==0) {
        accumulate(*dot, x, y);
      } else if(k==1) {
        val += complex(Re(x)*y, Im(x)*y);
      } else if(k==2) {
        real a,b,c;
        TwoProduct(Re(x),y,a,b);
        TwoSum(Re(val),a,Re(val),c);
        Re(corr) += (b+c);
        TwoProduct(Im(x),y,a,b);
        TwoSum(Im(val),a,Im(val),c);
        Im(corr) += (b+c);
      } else if(k>=3) {
        real a;
        complex b,c;
        TwoProduct(Re(x),y,a,Re(b));
        TwoSum(Re(val),a,Re(val),Re(c));
        TwoProduct(Im(x),y,a,Im(b));
        TwoSum(Im(val),a,Im(val),Im(c));
        cm.push_back(b);
        ca.push_back(c);
      }
    }

    void sparse_cdot::add_dot(const real& y, const complex& x) {
      if(k==0) {
        accumulate(*dot, x, y);
      } else if(k==1) {
        val += complex(Re(x)*y, Im(x)*y);
      } else if(k==2) {
        real a,b,c;
        TwoProduct(Re(x),y,a,b);
        TwoSum(Re(val),a,Re(val),c);
        Re(corr) += (b+c);
        TwoProduct(Im(x),y,a,b);
        TwoSum(Im(val),a,Im(val),c);
        Im(corr) += (b+c);
      } else if(k>=3) {
        real a;
        complex b,c;
        TwoProduct(Re(x),y,a,Re(b));
        TwoSum(Re(val),a,Re(val),Re(c));
        TwoProduct(Im(x),y,a,Im(b));
        TwoSum(Im(val),a,Im(val),Im(c));
        cm.push_back(b);
        ca.push_back(c);
      }
    }

    void sparse_cdot::add_dot_err(const complex& x, const complex& y) {
      if(k==0) {
        accumulate(*dot, x, y);
      } else if(k==1) {
        cm.push_back(x);
        ca.push_back(y);
      } else if(k==2) {
        real a,b,c;
        TwoProduct(Re(x),Re(y),a,b);
        TwoSum(Re(val),a,Re(val),c);
        c += b;
        Re(corr) += c;
        Re(err) += abs(c);
        TwoProduct(-Im(x),Im(y),a,b);
        TwoSum(Re(val),a,Re(val),c);
        c += b;
        Re(corr) += c;
        Re(err) += abs(c);

        TwoProduct(Re(x),Im(y),a,b);
        TwoSum(Im(val),a,Im(val),c);
        c += b;    
        Im(corr) += c;
        Im(err) += abs(c);
        TwoProduct(Im(x),Re(y),a,b);
        TwoSum(Im(val),a,Im(val),c);
        c += b;
        Im(corr) += c;
        Im(err) += abs(c);
        n++;
      } else if(k>=3) {
        real a;
        complex b,c;
        TwoProduct(Re(x),Re(y),a,Re(b));
        TwoSum(Re(val),a,Re(val),Re(c));
        TwoProduct(Re(x),Im(y),a,Im(b));
        TwoSum(Im(val),a,Im(val),Im(c));
        cm.push_back(b);
        ca.push_back(c);

        TwoProduct(-Im(x),Im(y),a,Re(b));
        TwoSum(Re(val),a,Re(val),Re(c));
        TwoProduct(Im(x),Re(y),a,Im(b));
        TwoSum(Im(val),a,Im(val),Im(c));
        cm.push_back(b);
        ca.push_back(c);
      }
    }

    void sparse_cdot::add_dot_err(const complex& x, const real& y) {
      if(k==0) {
        accumulate(*dot, x, y);
      } else if(k==1) {
        cm.push_back(complex(y));
        ca.push_back(x);
      } else if(k==2) {
        real a,b,c;
        TwoProduct(Re(x),y,a,b);
        TwoSum(Re(val),a,Re(val),c);
        c+=b;
        Re(err) += abs(c);
        Re(corr) += c;
        TwoProduct(Im(x),y,a,b);
        TwoSum(Im(val),a,Im(val),c);
        c += b;
        Im(corr) += c;
        Im(err) += abs(b);
        n++;
      } else if(k>=3) {
        real a;
        complex b,c;
        TwoProduct(Re(x),y,a,Re(b));
        TwoSum(Re(val),a,Re(val),Re(c));
        TwoProduct(Im(x),y,a,Im(b));
        TwoSum(Im(val),a,Im(val),Im(c));
        cm.push_back(b);
        ca.push_back(c);
      }
    }

    void sparse_cdot::add_dot_err(const real& y, const complex& x) {
      if(k==0) {
        accumulate(*dot, x, y);
      } else if(k==1) {
        cm.push_back(complex(y));
        ca.push_back(x);
      } else if(k==2) {
        real a,b,c;
        TwoProduct(Re(x),y,a,b);
        TwoSum(Re(val),a,Re(val),c);
        c+=b;
        Re(err) += abs(c);
        Re(corr) += c;
        TwoProduct(Im(x),y,a,b);
        TwoSum(Im(val),a,Im(val),c);
        c += b;
        Im(corr) += c;
        Im(err) += abs(b);
        n++;
      } else if(k>=3) {
        real a;
        complex b,c;
        TwoProduct(Re(x),y,a,Re(b));
        TwoSum(Re(val),a,Re(val),Re(c));
        TwoProduct(Im(x),y,a,Im(b));
        TwoSum(Im(val),a,Im(val),Im(c));
        cm.push_back(b);
        ca.push_back(c);
      }
    }

    complex sparse_cdot::result() {
      if(k==0) {
        return rnd(*dot);
      } else if(k==1) {
        return val;
      } else if(k==2) {
        return val+corr;
      } else if(k>=3) {
        int n = cm.size();
        if(n==0) return val;

	for(int j=1 ; j<k-1 ; j++) {
		for(int i=1 ; i<n ; i++) {
			TwoSum(Re(cm[i]),Re(cm[i-1]),Re(cm[i]),Re(cm[i-1]));    
			TwoSum(Im(cm[i]),Im(cm[i-1]),Im(cm[i]),Im(cm[i-1]));    
                }
              	TwoSum(Re(ca[0]),Re(cm[n-1]),Re(ca[0]),Re(cm[n-1]));    
              	TwoSum(Im(ca[0]),Im(cm[n-1]),Im(ca[0]),Im(cm[n-1]));    
		for(int i=1 ; i<n ; i++) {
			TwoSum(Re(ca[i]),Re(ca[i-1]),Re(ca[i]),Re(ca[i-1]));    
			TwoSum(Im(ca[i]),Im(ca[i-1]),Im(ca[i]),Im(ca[i-1]));    
                }
		TwoSum(Re(val),Re(ca[n-1]),Re(val),Re(ca[n-1]));
		TwoSum(Im(val),Im(ca[n-1]),Im(val),Im(ca[n-1]));
	}
		
	corr = std::accumulate(cm.begin(),cm.end(),corr);
	corr = std::accumulate(ca.begin(),ca.end(),corr);
		
	val += corr;

        return val;
      }
      return val;
    }

    void sparse_cdot::result(cdotprecision& res_dot) {
      if(k==0) {

        res_dot += *dot;

      } else if(k==1) {

        complex resd(0.0), resu(0.0);

#ifndef _CXSC_DOTK_ROUND_SOFT
        setround(-1);
        for(unsigned int i=0 ; i<cm.size() ; i++)
           resd += complex(Re(cm[i])*Re(ca[i])-Im(cm[i])*Im(ca[i]), Re(cm[i])*Im(ca[i])+Im(cm[i])*Re(ca[i]));

        setround(1);
        for(unsigned int i=0 ; i<cm.size() ; i++)
           resu += complex(Re(cm[i])*Re(ca[i])-Im(cm[i])*Im(ca[i]), Re(cm[i])*Im(ca[i])+Im(cm[i])*Re(ca[i]));

        val = resd+(resu-resd)*0.5;
        err = val - resd;

        Re(res_dot).set_err(Re(res_dot).get_err() + Re(err)); 
        Im(res_dot).set_err(Im(res_dot).get_err() + Im(err)); 

        setround(0);
#else
        for(unsigned int i=0 ; i<cm.size() ; i++)
           resd = complex(addd(Re(resd), subd(muld(Re(cm[i]),Re(ca[i])) , muld(Im(cm[i]),Im(ca[i])) ) ), 
                          addd(Im(resd), addd(muld(Re(cm[i]),Im(ca[i])) , muld(Im(cm[i]),Re(ca[i])) ) ));

        for(unsigned int i=0 ; i<cm.size() ; i++)
           resu = complex(addu(Re(resu), subu(mulu(Re(cm[i]),Re(ca[i])) , mulu(Im(cm[i]),Im(ca[i])) ) ), 
                          addu(Im(resu), addu(mulu(Re(cm[i]),Im(ca[i])) , mulu(Im(cm[i]),Re(ca[i])) ) ));

        val = resd+(resu-resd)*0.5;
        err = subu(val,resd);

        Re(res_dot).set_err(addu(Re(res_dot).get_err() , Re(err))); 
        Im(res_dot).set_err(addu(Im(res_dot).get_err() , Im(err))); 
#endif
        
        res_dot += val;

      } else if(k==2) {

	real alpha_re, alpha_im, delta, error_re, error_im;
	delta = (2*n*Epsilon) / (1.0-4*n*Epsilon);

	alpha_re = (Epsilon*abs(Re(val))) + (delta*Re(err)+3*MinReal/Epsilon);
	error_re = alpha_re / (1.0 - 2*Epsilon);      
        Re(res_dot).set_err(addu(Re(res_dot).get_err(), error_re)); 

	alpha_im = (Epsilon*abs(Im(val))) + (delta*Im(err)+3*MinReal/Epsilon);      
	error_im = alpha_im / (1.0 - 2*Epsilon);            
        Im(res_dot).set_err(addu(Im(res_dot).get_err(), error_im)); 

        res_dot += val;
        res_dot += corr;

      } else if(k>=3) {

        n = cm.size();
        if(n==0) return;

	for(int j=1 ; j<k-1 ; j++) {
		for(int i=1 ; i<n ; i++) {
			TwoSum(Re(cm[i]),Re(cm[i-1]),Re(cm[i]),Re(cm[i-1]));    
			TwoSum(Im(cm[i]),Im(cm[i-1]),Im(cm[i]),Im(cm[i-1]));    
                }
              	TwoSum(Re(ca[0]),Re(cm[n-1]),Re(ca[0]),Re(cm[n-1]));    
              	TwoSum(Im(ca[0]),Im(cm[n-1]),Im(ca[0]),Im(cm[n-1]));    
		for(int i=1 ; i<n ; i++) {
			TwoSum(Re(ca[i]),Re(ca[i-1]),Re(ca[i]),Re(ca[i-1]));    
			TwoSum(Im(ca[i]),Im(ca[i-1]),Im(ca[i]),Im(ca[i-1]));    
                }
		TwoSum(Re(val),Re(ca[n-1]),Re(val),Re(ca[n-1]));
		TwoSum(Im(val),Im(ca[n-1]),Im(val),Im(ca[n-1]));
	}
		
	corr = std::accumulate(cm.begin(),cm.end(),corr);
	corr = std::accumulate(ca.begin(),ca.end(),corr);

        res_dot += val;
        res_dot += corr;
	
        complex tmperr(0.0);

        for(unsigned int j=0 ; j<cm.size() ; j++) {
          Re(tmperr) += abs(Re(cm[j]));
          Im(tmperr) += abs(Im(cm[j]));
        }

        for(unsigned int j=0 ; j<ca.size() ; j++) {
          Re(tmperr) += abs(Re(ca[j]));
          Im(tmperr) += abs(Im(ca[j]));
        }

        complex res = val + corr;
		
	real alpha, delta, error;
	
	delta = (2*n*Epsilon) / (1.0-4*n*Epsilon);

	alpha = (Epsilon*abs(Re(res))) + (delta*Re(tmperr)+3*MinReal/Epsilon);
	error = alpha / (1.0 - 2*Epsilon);
        Re(res_dot).set_err(addu(Re(res_dot).get_err(),error));

	alpha = (Epsilon*abs(Im(res))) + (delta*Im(tmperr)+3*MinReal/Epsilon);
	error = alpha / (1.0 - 2*Epsilon);
        Im(res_dot).set_err(addu(Im(res_dot).get_err(),error));
	
      }

    }



}
 
 
 
