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

/* CVS $Id: cdot.cpp,v 1.25 2014/01/30 17:23:43 cxsc Exp $ */

#include "cdot.hpp"
#include "ioflags.hpp"

namespace cxsc {

//cdotprecision cdotakku[MAXCDOTAKKU]; 

// ---- Ausgabefunkt. ---------------------------------------

std::ostream & operator << (std::ostream &s, const cdotprecision& a) throw()
{
   s << '('
     << a.re << ',' 
     << a.im        
     << ')';
   return s;
}
std::string & operator << (std::string &s, const cdotprecision& a) throw()
{
   s += '(';
   s << a.re;
   s += ',';
   s << a.im;
   s += ')';
   return s;
}

std::istream & operator >> (std::istream &s, cdotprecision &a) throw()
{
   char c;

   skipeolnflag = inpdotflag = true;
   c = skipwhitespacessinglechar (s, '(');
   if (inpdotflag)
      s.putback(c);

   s >> a.re;

   skipeolnflag = inpdotflag = true;
   c = skipwhitespacessinglechar (s, ',');
   if (inpdotflag) 
      s.putback(c);

   s >> a.im;

   if (!waseolnflag)
   {
      skipeolnflag = false, inpdotflag = true;
      c = skipwhitespaces (s);
      if (inpdotflag && c != ')')
         s.putback(c);
   }

   return s;
}
        


std::string & operator >> (std::string &s, cdotprecision &a) throw()
{
   s = skipwhitespacessinglechar (s, '(');
   s >> a.re;
   s = skipwhitespacessinglechar (s, ',');
   s >> a.im;
   s = skipwhitespaces (s);

   if (s[0] == ')')
      s.erase(0,1);

   return s;
}

void operator >>(const std::string &s,cdotprecision &a) throw()
{
   std::string r(s);
   r>>a;
}

void operator >>(const char *s,cdotprecision &a) throw()
{
   std::string r(s);
   r>>a;
}                          

void rnd(const cdotprecision &a,complex &b,rndtype r) throw()
{
   Re(b)=rnd(a.re,r);
   Im(b)=rnd(a.im,r);
}

void rnd(const cdotprecision &a,complex &b,complex &c) throw()
{
   rnd(a,b,RND_DOWN);
   rnd(a,c,RND_UP);
}

void rnd (const cdotprecision& d, cinterval& x) throw() 
{
    complex a,b;
    rnd(d,a,b);
    x = cinterval(a,b);
}

complex rnd(const cdotprecision &a,rndtype r) throw()
{
   complex b;
   rnd(a,b,r);
   return b;
}

void accumulate(cdotprecision & a, const complex & b, const complex & c) throw()
{
   c_padd(a.re.ptr(),a.im.ptr(), *(a_cmpx*)&b,*(a_cmpx*)&c);
}

} // namespace cxsc

