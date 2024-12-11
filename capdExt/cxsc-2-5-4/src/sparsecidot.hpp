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

/* CVS $Id: sparsecidot.hpp,v 1.9 2014/01/30 17:23:49 cxsc Exp $ */

#ifndef _CXSC_SPARSECIDOT_HEADER
#define _CXSC_SPARSECIDOT_HEADER

#include <cinterval.hpp>
#include <vector>
#include <numeric>
#include "sparseidot.hpp"

namespace cxsc {

class sparse_cidot {
  private:
    sparse_idot re;
    sparse_idot im;

  public:
    
    sparse_cidot(unsigned int p) : re(p),im(p) {
    }

    sparse_cidot(unsigned int p, int nnz) : re(p,nnz), im(p,nnz) {
    }

    void reset() {
      re.reset();
      im.reset();
    }

    void add_dot(const cinterval& x, const cinterval& y) {
      re.add_dot(Re(x), Re(y));
      re.add_dot(-Im(x), Im(y));
      im.add_dot(Re(x), Im(y));
      im.add_dot(Im(x),Re(y));
    }

    void add_dot(const cinterval& x, const interval& y) {
      re.add_dot(Re(x), y);
      im.add_dot(Im(x), y);
    }

    void add_dot(const interval& x, const cinterval& y) {
      re.add_dot(x, Re(y));
      im.add_dot(x, Im(y));
    }

    void add_dot(const cinterval& x, const complex& y) {
      re.add_dot(Re(x), Re(y));
      re.add_dot(-Im(x), Im(y));
      im.add_dot(Re(x), Im(y));
      im.add_dot(Im(x),Re(y));
    }

    void add_dot(const complex& x, const cinterval& y) {
      re.add_dot(Re(x), Re(y));
      re.add_dot(-Im(x), Im(y));
      im.add_dot(Re(x), Im(y));
      im.add_dot(Im(x),Re(y));
    }

    void add_dot(const cinterval& x, const real& y) {
      re.add_dot(Re(x), y);
      im.add_dot(Im(x), y);
    }

    void add_dot(const real& x, const cinterval& y) {
      re.add_dot(x, Re(y));
      im.add_dot(x, Im(y));
    }

    void add_dot(const complex& x, const interval& y) {
      re.add_dot(Re(x), y);
      im.add_dot(Im(x), y);
    }

    void add_dot(const interval& x, const complex& y) {
      re.add_dot(x, Re(y));
      im.add_dot(x, Im(y));
    }

    void add_dot_err(const cinterval& x, const cinterval& y) {
      re.add_dot_err(Re(x), Re(y));
      re.add_dot_err(-Im(x), Im(y));
      im.add_dot_err(Re(x), Im(y));
      im.add_dot_err(Im(x),Re(y));
    }

    void add_dot_err(const cinterval& x, const interval& y) {
      re.add_dot_err(Re(x), y);
      im.add_dot_err(Im(x), y);
    }

    void add_dot_err(const interval& x, const cinterval& y) {
      re.add_dot_err(x, Re(y));
      im.add_dot_err(x, Im(y));
    }

    void add_dot_err(const cinterval& x, const complex& y) {
      re.add_dot_err(Re(x), Re(y));
      re.add_dot_err(-Im(x), Im(y));
      im.add_dot_err(Re(x), Im(y));
      im.add_dot_err(Im(x),Re(y));
    }

    void add_dot_err(const complex& x, const cinterval& y) {
      re.add_dot_err(Re(x), Re(y));
      re.add_dot_err(-Im(x), Im(y));
      im.add_dot_err(Re(x), Im(y));
      im.add_dot_err(Im(x),Re(y));
    }

    void add_dot_err(const cinterval& x, const real& y) {
      re.add_dot_err(Re(x), y);
      im.add_dot_err(Im(x), y);
    }

    void add_dot_err(const real& x, const cinterval& y) {
      re.add_dot_err(x, Re(y));
      im.add_dot_err(x, Im(y));
    }

    void add_dot_err(const complex& x, const interval& y) {
      re.add_dot_err(Re(x), y);
      im.add_dot_err(Im(x), y);
    }

    void add_dot_err(const interval& x, const complex& y) {
      re.add_dot_err(x, Re(y));
      im.add_dot_err(x, Im(y));
    }

    cinterval result() {
      return cinterval(re.result(), im.result());
    }

    void result(cidotprecision& dot) {
      idotprecision tmp1(0.0), tmp2(0.0);
      re.result(tmp1);
      im.result(tmp2);
      dot += cidotprecision(tmp1,tmp2);
    }

};

}

#endif

 
