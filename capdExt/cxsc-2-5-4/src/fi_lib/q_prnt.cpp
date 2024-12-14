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

/* CVS $Id: q_prnt.cpp,v 1.15 2014/01/30 17:23:55 cxsc Exp $ */

#ifndef Q_PRNT_CPP
#define Q_PRNT_CPP

#include <cmath>
#include "fi_lib.hpp" 

namespace fi_lib{

 using namespace std;
 using namespace cxsc;

 void printup(real x)
 {  
  cout << SetPrecision(23,15) << Scientific;
  if (x==(real)CUTINT(x))
    cout << x;
  else if(fabs(_double(x))>(1e17-1)*1e27) 
    cout << q_succ(q_succ(x));
  else 
    cout << q_succ(x);
 }

 void printdown(real x)
 {
  cout << SetPrecision(23,15) << Scientific;
  if (x==(real)CUTINT(x)) 
    cout << x;
  else if(fabs(_double(x))>(1e17-1)*1e27) 
    cout << q_pred(q_pred(x));
  else 
    cout << q_pred(x);
  
 }

 void printInterval(interval x)
 {
    cout << x; 
 }

} // Namespace

#endif





