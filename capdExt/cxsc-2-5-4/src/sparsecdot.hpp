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

/* CVS $Id: sparsecdot.hpp,v 1.8 2014/01/30 17:23:49 cxsc Exp $ */

#ifndef _CXSC_SPARSECDOT_HEADER
#define _CXSC_SPARSECDOT_HEADER

#include <complex.hpp>
#include <cdot.hpp>
#include <sparsedot.hpp>

namespace cxsc {

class sparse_cdot {
  private:
    cdotprecision* dot;
    std::vector<complex> cm;
    std::vector<complex> ca;
    complex val;
    complex corr; 
    complex err;
    int n;
    int k;

  public:
    
    sparse_cdot(unsigned int p);

    sparse_cdot(const sparse_cdot& s);

    ~sparse_cdot();

    void reset();

    void add_dot(const complex& x, const complex& y);

    void add_dot(const complex& x, const real& y);

    void add_dot(const real& y, const complex& x);

    void add_dot_err(const complex& x, const complex& y);

    void add_dot_err(const complex& x, const real& y);

    void add_dot_err(const real& y, const complex& x);

    complex result();

    void result(cdotprecision& res_dot);

};

}

#endif

 
 
