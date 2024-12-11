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

/* CVS $Id: sparseidot.hpp,v 1.10 2014/01/30 17:23:49 cxsc Exp $ */

#ifndef _CXSC_SPARSEIDOT_HPP_INCLUDED
#define _CXSC_SPARSEIDOT_HPP_INCLUDED

#include <interval.hpp>
#include <idot.hpp>
#include <vector>
#include <numeric>
#include "sparsedot.hpp"

namespace cxsc {

class sparse_idot {
  private:
    idotprecision* dot;
    std::vector<real> cm_x;
    std::vector<real> cm_y;
    std::vector<real> ca_x;
    std::vector<real> ca_y;
    interval val;
    real corr_inf;
    real corr_sup;
    real err_inf;
    real err_sup;
    int k;
    int n;

  public:
    
    sparse_idot(unsigned int p);

    sparse_idot(unsigned int p, int nnz);

    sparse_idot(const sparse_idot& s);

    ~sparse_idot();

    void reset();

    void add_dot(const interval& x, const interval& y);

    void add_dot(const interval& x, const real& y);

    void add_dot(const real& x, const interval& y);

    void add_dot_err(const interval& x, const interval& y);

    void add_dot_err(const interval& x, const real& y);

    void add_dot_err(const real& x, const interval& y);

    interval result();

    void result(idotprecision& res_dot);

};

}

#endif

 
