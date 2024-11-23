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

/* CVS $Id: realio.cpp,v 1.31 2014/01/30 17:23:48 cxsc Exp $ */

#include <cstring>
#include "real.hpp"
#include "ioflags.hpp"
#include "RtsFunc.h"
#include "dot.hpp"

namespace cxsc {

int d_init_dm (void);
void d_outp(char *buffer, Dotprecision c,
              int FormatFlag, int FracDigits, int rnd,
              int *length);

#if _WIN32
extern __declspec(thread) char *dm;
#elif __APPLE__ && !CXSC_FORCE_TLS
extern char *dm; 
#else
extern __thread char *dm;
#endif


char* GetHexDigits (char* s, a_btyp& x, int count);

int r_outpx (char *buffer, a_real s,a_intg FormatFlag, a_intg FracDigits, a_intg rnd, a_intg *length);

#define WASGOOD          0 
#define ISINFINITYPLUS  -1
#define ISINFINITYMINUS -2
#define ISQUIETNAN      -3
#define ISSIGNALINGNAN -4

std::string realToHex(const real& a) 
{
   a_btyp* b = (a_btyp*) &a;
   int i;
   char cs[256];

   for (i=0; i < realwidth-19; i++) cs[i] = ' ';
   cs[i] = 0;
   sprintf (&cs[strlen(cs)], "%s", ((b[HIGHREAL] & 0x80000000L) ? "-" : "+"));
   sprintf (&cs[strlen(cs)],"%c",   '1');
   sprintf (&cs[strlen(cs)],"%05lX", (b[HIGHREAL] & 0xFFFFFL));
   sprintf (&cs[strlen(cs)],"%08lX",(long unsigned int)(b[LOWREAL]));
   sprintf (&cs[strlen(cs)],"e%03X",(unsigned int)((b[HIGHREAL] >> 20) & 0x7FF));
   return string(cs);
}

std::string & operator<<(std::string &s, const real& a) throw()
{
#if _WIN32
   static __declspec(thread) char cs[256];
#elif __APPLE__ && !CXSC_FORCE_TLS
    static char cs[256];
#else
   static __thread char cs[256];
#endif
   
   if (ioflags.isset(IOFlags::hex))
   {
      a_btyp* b = (a_btyp*) &a;
      int i;
      char cs[256];

      for (i=0; i < realwidth-19; i++) cs[i] = ' ';
      cs[i] = 0;
      sprintf (&cs[strlen(cs)], "%s", ((b[HIGHREAL] & 0x80000000L) ? "-" : "+"));
      sprintf (&cs[strlen(cs)],"%c",   '1');
      sprintf (&cs[strlen(cs)],"%05lX", (b[HIGHREAL] & 0xFFFFFL));
      sprintf (&cs[strlen(cs)],"%08lX",(long unsigned int)(b[LOWREAL]));
      sprintf (&cs[strlen(cs)],"e%03X",(unsigned int)((b[HIGHREAL] >> 20) & 0x7FF));
      s+=cs;
   } else if (ioflags.isset(IOFlags::rndnone))
   {
      if (IsSignalingNaN(a)) s+="<SignallingNaN>";
      else if (IsQuietNaN(a)) s+="<QuietNaN>";
      else if (IsInfinity(a)) s+="<Infinity>";
      else 
      {
         if (realdigits && realwidth)
            sprintf (cs, "%*.*g", realwidth, realdigits, a.w); // no need for "lg"
         else if (realwidth)
            sprintf (cs, "%*g", realwidth, a.w);
         else
            sprintf (cs, "%g", a.w);
         s+=cs;
      }
   } else
   {
      rndtype rnd;
      a_intg length, formatflag, addblanks;
      a_intg digits = realdigits;
      char *str;

      if (d_init_dm () == -1) 
      { 
         // throw
         // errmon (ERR_ALL(NOMOREMEMORY));
         // errmon (ERR_ALL(NOCONTINUEPOSSIBLE));
      }

      if (ioflags.isset(IOFlags::rndup)) 
         rnd = RND_UP;
      else if (ioflags.isset(IOFlags::rnddown)) 
         rnd = RND_DOWN;
      else 
         rnd = RND_NEXT;

      if (ioflags.isset(IOFlags::variable)) 
         formatflag = realwidth;
      else if (ioflags.isset(IOFlags::varfixwidth)) 
         formatflag = realwidth, digits = -digits;
      else 
         formatflag = (ioflags.isset(IOFlags::fixed)) ? 0 : -1;

      switch (r_outpx (dm, a.w, formatflag, digits, rnd, &length)) 
      { 
         case WASGOOD:
            dm[length] = 0;
            str = dm;
            if (*str == '+') 
            {
               if (ioflags.isset(IOFlags::blank))
                  *str = ' ';
               else if (ioflags.isset(IOFlags::noblank)) 
                  str++;
            }
            break;
         case ISINFINITYPLUS:  
            str = (char*)"<+Infinity>";     
            break;
         case ISINFINITYMINUS: 
            str = (char*)"<-Infinity>";     
            break;
         case ISQUIETNAN:      
            str = (char*)"<QuietNaN>";      
            break;
         case ISSIGNALINGNAN:  
            str = (char*)"<SignalingNaN>";  
            break;
         default:              
            str = (char*)"<ERROR>";
            break;
      }
      length = strlen(str);
      addblanks = (length < realwidth) ? realwidth - length : 0;

      if (ioflags.isset(IOFlags::rightjust)) 
         for (;addblanks; addblanks--) 
            s+= ' ';

      s+=str;
      
      for (;addblanks; addblanks--) 
         s+= ' ';
   }
   return s;
}

std::ostream & operator <<(std::ostream &o,const real &a) throw()
{
   std::string s="";
   s << a;
   o << s;
   return o;
}

std::string & operator>> (std::string & str, real& a) throw()
{
   char *s=new char[str.size()+1];
   char *orgs=s;
   strcpy(s,str.c_str());  

   if (ioflags.isset(IOFlags::hex))
   {
      a_btyp* b = (a_btyp*) &a;
      a_btyp x;

      b[0] = b[1] = 0;
      s = cskipwhitespaces (s);
      if (*s == '-') 
      {
         b[HIGHREAL] |= 0x80000000L; 
         s++; 
      } else if (*s == '+') 
         s++;

    // ----------------------------------------------------
    // es wird ohne Pruefung folgendes Format vorausgestzt :
    //    1xxxxxxxxxxxxxeXXX
    // wobei x = HexDigits der Mantisse
    //       X = HexDigits des Exponents
                      // saemmtliche HexDigits muessen gross geschrieben sein !

      if (*s) 
         s++;                                      // skip '1'
      s = GetHexDigits (s, x, 5);  
      b[HIGHREAL] |= x;    // get mantissa
      s = GetHexDigits (s, x, 8);  
      b[LOWREAL]   = x;
      if (*s) 
         s++;                                      // skip 'e'     
      s = GetHexDigits (s, x, 3);
      b[HIGHREAL] |= x << 20;                 // get exponent
      if (*s) 
         s++;                            // skip at least one more char
   } else
   {
      rndtype rndfl;

      if (ioflags.isset(IOFlags::rndup)) 
         rndfl = RND_UP;
      else if (ioflags.isset(IOFlags::rnddown)) 
         rndfl = RND_DOWN;
      else 
         rndfl = RND_NEXT;

      str=s;
      dotprecision dot;
      str >> dot;
      strcpy(s,str.c_str()); // Ooooooooohhhh... :((
      a = rnd (dot, rndfl);
   }
   str=s;
   delete [] orgs;
   return str;
}   
void operator >>(const char *a,real &b) throw()
{
   std::string c(a);
   c>>b;
}
void operator >>(const string &a,real &b) throw()
{
   std::string c(a);
   c>>b;
}
std::istream& operator>> (std::istream& s, real& a) throw()
{
   if (ioflags.isset(IOFlags::hex))
   {
      char inp[20];
      int i;
      char c;

      // !! There are no checks about the right input-format, it is assumed
      //    the fixed format : S1xxxxxxxxxxxxxeXXX
      //    whereby S   means the sign
      //            '1' is the leading implicit binary one of the mantissa
      //            x   are the mantissa hexdigits
      //            'e' signals the beginning of the exponent
      //            X   are the exponent hexdigits

      // skip white spaces and
      // get the sign, the 14 Mantissadigits and the exponent
      // (total 19 chars)
                                               
      c = skipwhitespaces (s);
      for (i=0; i < 19; i++) 
      {
         inp[i] = c;
         if (s.good()) s.get(c); else c = '\0';
      }
      inp[i] = '\0';

      inp >> a;
   } else
   {
      rndtype rndfl;

      if (ioflags.isset(IOFlags::rndup)) 
         rndfl = RND_UP;
      else if (ioflags.isset(IOFlags::rnddown)) 
         rndfl = RND_DOWN;
      else 
         rndfl = RND_NEXT;

      dotprecision dot;
      s >> dot;
      a = rnd (dot, rndfl);
   }
   return s;
}                                                                                                                    

//----------------------------------------------------------------------------
// GetHexDigits
//
//  Interpretiert die naechsten "count" Zeichen aus dem String "s"
//  als HexDigits, der binaere Wert wird dann in "x" zurueckgegeben.
//  Als Returnwert wird ein Zeiger auf die Stringposition nach den
//  HexDigits zurueckgegeben.

char* GetHexDigits (char* s, a_btyp& x, int count)
{
  int i, c;

  for (x=0,i=0; i < count && *s; s++,i++) {
    if ((c = *s) >= 'A') c -= 'A' - 10; else c -= '0';
    if (c < 0 || c > 0xF) c = 0;
    x = (x << 4) | c;
  }
  return s;
}        

} // namespace cxsc

/****************************************************************/
/*                                                              */
/*      Filename        : r_outpx.c                             */
/*                                                              */
/*      Entries         : void r_outp                           */
/*                         (buffer,s,FormatFlag,                */
/*                          FracDigits,rnd,length)              */
/*                        char *buffer;                         */
/*                        a_real s;                             */
/*                        a_intg rnd,FormatFlag,FracDigits;     */
/*                        a_intg *length;                       */
/*                                                              */
/*      Arguments       : buffer - output buffer holding string */
/*                        s - IEEE value                        */
/*                        FormatFlag - format selection         */
/*                              -1 = scientific format          */
/*                               0 = fixed format               */
/*                             > 0 = variable format            */
/*                                   value is the total field   */
/*                                   width                      */
/*                        FracDigits - number of fraction digits*/
/*                        rnd - rounding mode                   */
/*                              -1 = round downwards            */
/*                               0 = round to nearest           */
/*                               1 = round upwards              */
/*                        length - size of buffer string        */
/*                                                              */
/*      Description     : Convert an IEEE double format number  */
/*                        to a character string.                */
/*                                                              */
/****************************************************************/

/*#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
#endif*/

namespace cxsc {

#define WASGOOD          0 
#define ISINFINITYPLUS  -1
#define ISINFINITYMINUS -2
#define ISQUIETNAN      -3
#define ISSIGNALINGNAN -4

#define MANT_INFINITY(a)        ((a)[0]==HIDDEN_BIT && (a)[1]==ZERO)
#define SIGNALING(a)            ((a) & SIGNAL_BIT)
#define SIGNAL_BIT              ((a_btyp)0x00080000L)  

int r_outpx(char *buffer, a_real s, a_intg FormatFlag,
            a_intg FracDigits, a_intg rnd, a_intg *length)
{
   a_intg ActWidth,DecPlaces,expo,IntDigits,MinNumChars;
   a_intg dexpo,digits,bdp,k,l,addpoint;
   a_intg HoldWidth = (FracDigits < 0);
   a_bool vz,zero;
   a_btyp mant[BSIZE];

   if (HoldWidth) 
      FracDigits = -FracDigits;
   *length = 0;

   zero = b_deko(s,&expo,mant,&vz);

/*                                                                   */
/* infinity or NaN                                                   */
/*                                                                   */
   if (expo>EXPO_MAX) 
   {
      if (MANT_INFINITY(mant)) 
      {                         /* Infinity   */
         k = (vz) ? ISINFINITYMINUS : ISINFINITYPLUS;
      } else 
      {                                             /* NaN        */
         k = (SIGNALING(mant[0])) ? ISSIGNALINGNAN : ISQUIETNAN;
      }
      return k;
   }

   for (k=D_U_RATIO;k<BSIZE;k++) 
      mant[k] = ZERO;

   if (vz) 
      rnd = -rnd;

   if (FormatFlag > 0 && FormatFlag-FracDigits <= 2) 
   {
      FormatFlag = FracDigits+3;
   }
   if (!zero && FormatFlag > 0) 
   {
      dexpo = (a_intg) ((expo*30103L)/100000L);
      if (dexpo < -((FracDigits+1)/2) ||
         (dexpo > 0 && dexpo >= FormatFlag-FracDigits-2)) 
      {
         FormatFlag = -1;
         if (HoldWidth) 
         {
            FracDigits = (FracDigits < 5) ? 0 : FracDigits-5;
         }
      }
   }

   if (expo>800) 
      bdp = (BUFFERSIZE-8)-(FormatFlag == -1 ? 0 : FracDigits);
   else if (expo<-800) 
      bdp = 8;
   else 
      bdp = B_D_P;

   addpoint = (FracDigits == 0) ? 0 : 1;
/*                                                                   */
/* floating-point representation                                     */
/*                                                                   */
   if (FormatFlag == -1) 
   {
      DecPlaces = FracDigits;
      ActWidth =  DecPlaces+ExpDigits+4+addpoint;
/*                                                                   */
/* number is zero                                                    */
/*                                                                   */
      if (zero) 
      {
         *buffer++ = (vz) ? '-' : '+';
         *buffer++ = '0';
         if (FracDigits) 
         {
            *buffer++ = XSC_DECIMAL_POINT;
            for (k=0;k<DecPlaces;k++) *buffer++ = '0';
         }
         *buffer++ = EXPONENT_E;
         *buffer++ = PLUS_SIGN;
         for (k=0;k<ExpDigits;k++) 
            *buffer++ = '0';

         *length = ActWidth;

         return WASGOOD;
      }
/*                                                                   */
/* determine output string                                           */
/*                                                                   */
      digits = DecPlaces+(1+2);
      dexpo = -1;
      b_out(mant,expo,digits,buffer,&bdp,&dexpo);

      if (dexpo>0 && dexpo>DecPlaces+2) 
         digits = dexpo+1;

      b_rnd(rnd,buffer,digits,DecPlaces+1,&bdp,&dexpo);

      if (FracDigits) 
      {
         buffer[bdp-dexpo-1] = buffer[bdp-dexpo];
         buffer[bdp-dexpo] = XSC_DECIMAL_POINT;
      }
      buffer[bdp-dexpo+DecPlaces+1] = EXPONENT_E;
      buffer[bdp-dexpo+DecPlaces+2] = (dexpo<0) ? MINUS_SIGN : PLUS_SIGN;
      expo = (dexpo<0) ? -dexpo : dexpo;
      for (k=ExpDigits;k>0;k--) 
      {
         buffer[bdp-dexpo+DecPlaces+2+k] = expo%10+'0';
         expo /= 10;
      }
      *length = 2+addpoint+DecPlaces+2+ExpDigits;

      l = bdp - 1 - addpoint - dexpo;
      buffer[l] = (vz) ? '-' : '+';
      for (k=0;k<*length;k++,l++) 
         buffer[k] = buffer[l];
   } else
/*                                                                   */
/* fixed-point representation                                        */
/*                                                                   */
   {
/*                                                                   */
/* number is zero                                                    */
/*                                                                   */
      if (zero) 
      {
         *length = FracDigits+2+addpoint;
         *buffer++ = (vz) ? '-' : '+';
         *buffer++ = '0';
         if (FracDigits) 
         {
            *buffer++ = XSC_DECIMAL_POINT;
            for (k=0;k<FracDigits;k++) 
               *buffer++ = '0';
         }
         return WASGOOD;
      }

      /* estimate number of decimal digits */
      if (expo>=0) 
      {
         IntDigits = ((expo+1)*61)/200+1;
      } else 
      {
         IntDigits = 0;
         dexpo = 0;
      }

      /* fill fractional part with zeros                   */
      for (k=0;k<=FracDigits+2;k++) 
         buffer[bdp+k] = '0';

      digits = IntDigits+FracDigits+2;
      b_out(mant,expo,digits,buffer,&bdp,&dexpo);

      /* correct setting of IntDigits */
      if (expo>=0) 
      {
         IntDigits = dexpo+1;
      } else 
      {
         IntDigits = 1;
         digits++;
      }

      b_rnd(rnd,buffer,digits,IntDigits+FracDigits,&bdp,&dexpo);

      /* correct setting of IntDigits after rounding */
      IntDigits = dexpo+1;

      /* value is zero after rounding;                     */
      if (vz) 
      {
         for (k=(bdp+1)-IntDigits;k<(bdp+1)+FracDigits;k++) 
         {
            if (buffer[k]!='0') 
               break;
         }
         if (k==(bdp+1)+FracDigits) 
            vz = FALSE;
      }

      MinNumChars = IntDigits+FracDigits+1+addpoint;

      if (FracDigits) 
      {
         /* generate decimal point                            */
         for (k=bdp-IntDigits;k<bdp;k++) 
            buffer[k] = buffer[k+1];
         buffer[bdp] = XSC_DECIMAL_POINT;
      }

      *length = MinNumChars;
/*                                                                   */
/* sign of number                                                    */
/*                                                                   */
      l = bdp - 1 - addpoint - dexpo;
      buffer[l] = (vz) ? '-' : '+';
      for (k=0;k<*length;k++,l++) 
         buffer[k] = buffer[l];
   }

   return WASGOOD;
}

} // namespace cxsc

