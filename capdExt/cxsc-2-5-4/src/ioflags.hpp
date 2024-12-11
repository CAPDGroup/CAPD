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

/* CVS $Id: ioflags.hpp,v 1.27 2014/01/30 17:23:45 cxsc Exp $ */

#ifndef _CXSC_IOFLAGS_HPP_INCLUDED
#define _CXSC_IOFLAGS_HPP_INCLUDED

#include <iostream>

namespace cxsc {

#define CXSC_UNDEFINED -1

enum precisionflag { WIDTH, DIGITS };

class IOFlags
{
   private:
      unsigned long long int flagmask;
      int    width;
      int    digits;
   public:
      static const unsigned long long int dummy;
      static const unsigned long long int allmask;
      static const unsigned long long int typmask, hex, dec;
      static const unsigned long long int rndmask, rndnext, rndup, rnddown, rndnone;
      static const unsigned long long int realformat, dotformat;
      static const unsigned long long int saveopt, restoreopt;
      static const unsigned long long int rightjust, leftjust;
      static const unsigned long long int blank, noblank, forcesign;
      static const unsigned long long int fixed, scientific, variable, varfixwidth;
      static const unsigned long long int inoutmask, fullinout, sparseinout, matrixmarketinout;
  
      IOFlags();
      IOFlags(const unsigned long long int flag);
      IOFlags(const IOFlags&);
      
      void   setf(const unsigned long long int flag);
      void   setf(const unsigned long long int flag,const unsigned long long int mask);
      unsigned long long int getf(const unsigned long long int flag) const;
      int   isset(const unsigned long long int flag) const;

      void  setprecision (const int value = CXSC_UNDEFINED, const precisionflag sel = WIDTH);
      int   getprecision(const precisionflag sel=WIDTH) const;
      // setprecision
};

extern IOFlags ioflags;

enum rndtype { RND_DOWN = -1, RND_NEXT, RND_UP };

extern const IOFlags Hex, Dec, RndNext, RndUp, RndDown;
extern const IOFlags SaveOpt, RestoreOpt;
extern const IOFlags RealFormat, DotFormat;
extern const IOFlags RightJust, LeftJust;
extern const IOFlags Blank, NoBlank, ForceSign;
extern const IOFlags Fixed, Scientific, Variable, VarFixWidth;     
extern const IOFlags FullInOut, SparseInOut, MatrixMarketInOut;     

extern int realwidth, realdigits;
extern int dotwidth, dotdigits;
extern bool ScanDotRndFlag;

void SetOption (const IOFlags&);

IOFlags SetPrecision (
  const int width = CXSC_UNDEFINED, const int digits = CXSC_UNDEFINED);
int GetPrecision (const precisionflag sel = WIDTH);

IOFlags SetDotPrecision (
  const int width = CXSC_UNDEFINED, const int digits = CXSC_UNDEFINED);
int GetDotPrecision (const precisionflag sel = WIDTH);


  std::ostream& operator<< (std::ostream&, const IOFlags&);
  std::istream& operator>> (std::istream&, const IOFlags&);
  std::string&   operator<< (std::string&, const IOFlags&);
  std::string&   operator>> (std::string&, const IOFlags&);
       

// einige Hilfsroutinen fuer das einlesen
extern bool inpdotflag, skipeolnflag, waseolnflag;

   //------------------------------------------------------------------

    char skipwhitespaces (std::istream& s);
    char skipwhitespacessinglechar (std::istream& s, char);
    char skipleadingchars (std::istream& s, char beg, char end);
    void skipwhitespacesandputback (std::istream& s, char c);

    //------------------------------------------------------------------

    std::string skipwhitespaces (std::string);
    std::string skipwhitespacessinglechar (std::string, char);
    std::string skipleadingchars (std::string, char beg, char end);
    
    char * cskipwhitespaces (char *);
    char * cskipwhitespacessinglechar (char *, char);
    char * cskipleadingchars (char *, char beg, char end);
    
    //------------------------------------------------------------------

} // namespace cxsc 

#endif // _CXSC_IOFLAGS_HPP_INCLUDED
