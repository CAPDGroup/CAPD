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

/* CVS $Id: idot.cpp,v 1.23 2014/01/30 17:23:45 cxsc Exp $ */

#include "idot.hpp"
#include "ioflags.hpp"

namespace cxsc {

//idotprecision idotakku[MAXIDOTAKKU]; 
// ---- Ausgabefunkt. ---------------------------------------

std::ostream & operator << (std::ostream &s, const idotprecision& a) throw()
{
   s << '['          << SaveOpt << RndDown  
     << a.inf << ',' << RndUp 
     << a.sup        << RestoreOpt 
     << ']';
   return s;
}
std::string & operator << (std::string &s, const idotprecision& a) throw()
{
   s += '[';
   s << SaveOpt << RndDown  
     << a.inf;
   s += ',';
   s << RndUp 
     << a.sup   << RestoreOpt;
   s += ']';
   return s;
}

std::istream & operator >> (std::istream &s, idotprecision &a) throw()
{
   char c;

   skipeolnflag = inpdotflag = true;
   c = skipwhitespacessinglechar (s, '[');
   if (inpdotflag)
      s.putback(c);

   s >> SaveOpt >> RndDown >> a.inf;

   skipeolnflag = inpdotflag = true;
   c = skipwhitespacessinglechar (s, ',');
   if (inpdotflag) s.putback(c);

   s >> RndUp >> a.sup >> RestoreOpt;

   if (!waseolnflag)
   {
      skipeolnflag = false, inpdotflag = true;
      c = skipwhitespaces (s);
      if (inpdotflag && c != ']')
         s.putback(c);
   }

   /*if (a.inf > a.sup) {
      errmon (ERR_INTERVAL(EMPTY));
   } */
   return s;
}
        


std::string & operator >> (std::string &s, idotprecision &a) throw()
{
   s = skipwhitespacessinglechar (s, '[');
   s >> SaveOpt >> RndDown >> a.inf;
   s = skipwhitespacessinglechar (s, ',');
   s >> RndUp >> a.sup >> RestoreOpt;
   s = skipwhitespaces (s);

   if (s[0] == ']')
      s.erase(0,1);

    /*if (a.inf > a.sup) {
      errmon (ERR_INTERVAL(EMPTY));
    } */
   return s;
}

void operator >>(const std::string &s,idotprecision &a) throw()
{
   std::string r(s);
   r>>a;
}

void operator >>(const char *s,idotprecision &a) throw()
{
   std::string r(s);
   r>>a;
}                          

void rnd(const idotprecision &a,interval &b) throw()
{
   Inf(b)=rnd(a.inf,RND_DOWN);
   Sup(b)=rnd(a.sup,RND_UP);
}

interval rnd(const idotprecision &a) throw()
{
   interval b;
   rnd(a,b);
   return b;
}

void accumulate(idotprecision & a, const interval & b, const interval & c) throw()
{
   i_padd(a.inf.ptr(),a.sup.ptr(), *(a_intv*)&b,*(a_intv*)&c);
}

} // namespace cxsc

