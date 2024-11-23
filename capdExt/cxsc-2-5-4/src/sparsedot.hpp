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

/* CVS $Id: sparsedot.hpp,v 1.10 2014/01/30 17:23:49 cxsc Exp $ */

#ifndef _CXSC_SPARSEDOT_HPP_INCLUDED
#define _CXSC_SPARSEDOT_HPP_INCLUDED

#include <real.hpp>
#include <dot.hpp>
#include <vector>
#include <numeric>

namespace cxsc {

class sparse_dot {
  private:
    dotprecision* dot;
    std::vector<real> cm;
    std::vector<real> ca;
    real val;
    real corr;
    real err;
    int n;
    int k;

  public:
    
    sparse_dot(unsigned int p);

    sparse_dot(const sparse_dot& s);

    ~sparse_dot();

    void reset();

    void add_dot(const real& x, const real& y);

    void add_dot_err(const real& x, const real& y);

    real result();

    void result(dotprecision& res_dot);

};

}

#endif

 
