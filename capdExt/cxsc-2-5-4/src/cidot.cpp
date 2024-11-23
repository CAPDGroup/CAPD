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

/* CVS $Id: cidot.cpp,v 1.23 2014/01/30 17:23:43 cxsc Exp $ */

#include "cidot.hpp"
#include "ioflags.hpp"

namespace cxsc {

//cidotprecision cidotakku[MAXCIDOTAKKU];

// ---- Ausgabefunkt. ---------------------------------------

std::ostream & operator << (std::ostream &s, const cidotprecision& a) throw()
{
   s << '(' << SaveOpt         
     << '[' << RndDown << a.reinf << ',' << RndUp << a.resup << ']' << ','  
     << '[' << RndDown << a.iminf << ',' << RndUp << a.imsup << ']' 
     << ')' << RestoreOpt;
     
   return s;
}
std::string & operator << (std::string &s, const cidotprecision &a) throw()
{
   s+="([";
   s << SaveOpt << RndDown << a.reinf;
   s+=',';
   s << RndUp << a.resup;
   s+="],[";
   s << RndDown << a.iminf;
   s+=',';
   s << RndUp << a.imsup << RestoreOpt;
   s+="])";
   return s;
}

std::istream & operator >> (std::istream &s, cidotprecision &a) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   char c;

   skipwhitespacesandputback (s, '(');
   skipwhitespacesandputback (s, '[');
   s >> SaveOpt >> RndDown >> a.reinf;
   skipwhitespacesandputback (s, ',');
   s >> RndUp >>a.resup;
   skipwhitespacesandputback (s, ']');
   skipwhitespacesandputback (s, ',');
   skipwhitespacesandputback (s, '[');
   s >> RndDown >> a.iminf;
   skipwhitespacesandputback (s, ',');
   s >> RndUp >> a.imsup >> RestoreOpt;

   if (!waseolnflag) 
   {
      skipeolnflag = false, inpdotflag = true;
      c = skipwhitespaces (s);
      if (inpdotflag && c != ']') 
         s.putback(c);
   }
   if (!waseolnflag) 
   {
      skipeolnflag = false, inpdotflag = true;
      c = skipwhitespaces (s);
      if (inpdotflag && c != ')') 
         s.putback(c);
   }

   if (a.reinf > a.resup || a.iminf > a.imsup) 
      cxscthrow(ERROR_CIDOTPRECISION_EMPTY_INTERVAL("std::istream & operator >> (std::istream &s, cidotprecision &a)"));
      
   return s;
}

std::string & operator >> (std::string &s, cidotprecision &a) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   s = skipwhitespacessinglechar (s, '(');
   s = skipwhitespacessinglechar (s, '[');
   s = s >> SaveOpt >> RndDown >> a.reinf;
   s = skipwhitespacessinglechar (s, ',');
   s = s >> RndUp >> a.resup;
   s = skipwhitespacessinglechar (s, ']');
   s = skipwhitespacessinglechar (s, ',');
   s = skipwhitespacessinglechar (s, '[');
   s = s >> RndDown >> a.iminf;
   s = skipwhitespacessinglechar (s, ',');
   s = s >> RndUp >> a.iminf >> RestoreOpt;
   s = skipwhitespaces (s);
   if (s[0] == ']') 
      s.erase(0,1);
   s = skipwhitespaces (s);
   if (s[0] == ')') 
      s.erase(0,1);

   if (a.reinf > a.resup || a.iminf > a.imsup) 
      cxscthrow(ERROR_CIDOTPRECISION_EMPTY_INTERVAL("std::string & operator >> (std::string &s, cidotprecision &a)"));

   return s;
}

void operator >>(const std::string &s,cidotprecision &a) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   std::string r(s);
   r>>a;
}
void operator >>(const char *s,cidotprecision &a) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   std::string r(s);
   r>>a;
}


void accumulate(cidotprecision & a,const cinterval & b,const cinterval & c) throw()
{
   z_padd(
      a.reinf.ptr(),a.iminf.ptr(),
      a.resup.ptr(),a.imsup.ptr(),
      *(a_cinv*)&b,*(a_cinv*)&c);      
}

} // namespace cxsc

