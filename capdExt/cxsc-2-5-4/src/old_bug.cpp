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

/* CVS $Id: old_bug.cpp,v 1.21 2014/01/30 17:23:48 cxsc Exp $ */

#include <iostream.h>
#include <interval.hpp>
#include <l_interv.hpp>
#include <l_imath.hpp>
#include <l_rmath.hpp>
#include <rmath.hpp>

namespace cxsc {

int main(void)
{


   interval a(1,3);
   interval b(2,4);

   cout << "  interval: " << endl << a << " & " << b << " = " << (a&b) << endl; 

   // not possible: l_interval al(a),bl(b);
   l_interval al(Inf(a),Sup(a)),bl(Inf(b),Sup(b));

   cout << "l_interval: " << endl << al << " & " << bl << " = " << (al&bl) << endl;


   cout << endl;
   cout << "Setting Infimum of " << a << " to zero delivers:" << endl;
   Inf(a)=0;
   Inf(al)=0;
   cout << "  interval: " << a << endl;
   cout << "l_interval: " << al << endl;


   cout << endl;


   l_real q; // not possible: q(1);

   q=1; 

   cout << "Calculation of asinh(1) in l_real (or l_interval):" << endl;

   stagprec=1;
   cout << "Precision stagprec=1 [real]:" << endl;
   cout << asinh(q) << " " << asinh(_l_interval(q)) << endl;

   stagprec=2;
   cout << "Precision stagprec=2 [real]:" << endl;
   cout << asinh(q) << " " << asinh(_l_interval(q)) << endl;

   cout << "Result with rmath: asinh(1)=" << asinh(real(1)) << endl;


   return 0;
}

} // namespace cxsc
