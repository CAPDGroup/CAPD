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

/* CVS $Id: dotio.cpp,v 1.29 2014/01/30 17:23:45 cxsc Exp $ */

#include <iostream>
#include <string>
#include <cstring>

#include "dot.hpp"
#include "ioflags.hpp"
#include "RtsFunc.h"

#include "dot_defs.hpp"

namespace cxsc {

// #include "dot_defs.hpp"


int d_init_dm (void);
void d_outp(char *buffer, Dotprecision c,
                  int FormatFlag, int FracDigits, int rnd,
                  int *length);

#if _WIN32
extern  __declspec(thread) char *dm;
#elif __APPLE__ && !CXSC_FORCE_TLS
extern char *dm;
#else
extern __thread char *dm;
#endif


std::string & operator <<(std::string & s,const dotprecision &a) throw()
{
   if(ioflags.isset(IOFlags::realformat))
   {
      //char *sh = new char[1024];
      string sh;
      real rl,ru;
      rnd (a, rl, ru);    // Bem.: In rnd wird (*this) ggfl. entfernt

      sh="dot(";                   sh << SaveOpt << RndDown;
      /*sh << rl;*/ sh+=", ";      sh << RndUp;
      /*sh<< ru;*/  sh+=")";       sh << RestoreOpt;  
      s+=sh;    
      // delete [] sh;
   } else
   {
      rndtype rnd;
      int formatflag, addblanks, digits=dotdigits;
      int length;
      char *str;
      if (d_init_dm () == -1) 
      {
         // THROW!
         return s;
         //errmon (ERR_ALL(NOMOREMEMORY));
         //errmon (ERR_ALL(NOCONTINUEPOSSIBLE));
      }  
       
      if (ioflags.isset(IOFlags::rndup)) rnd = RND_UP;
      else if (ioflags.isset(IOFlags::rnddown)) rnd = RND_DOWN;
      else rnd = RND_NEXT;

      if (ioflags.isset(IOFlags::variable))
         formatflag = dotwidth;
      else if (ioflags.isset(IOFlags::varfixwidth))
         formatflag = dotwidth, digits = -digits;
      else
         formatflag = (ioflags.isset(IOFlags::fixed)) ? 0 : -1;
    
      d_outp (str = dm, a.akku, formatflag, digits, rnd, &length);
      if (*str == '+') 
      {
         if (ioflags.isset(IOFlags::blank))         *str = ' ';
         else if (ioflags.isset(IOFlags::noblank))  str++,length--;
      }
      addblanks = (length < dotwidth) ? dotwidth - length : 0;
      if (ioflags.isset(IOFlags::rightjust)) 
      {
         for (;addblanks; addblanks--) s+= ' ';
      }
      s+=str;
      if (!ioflags.isset(IOFlags::rightjust))
         for (;addblanks; addblanks--) s+= ' ';
      
   }
   return s;
}

std::ostream & operator <<(std::ostream & s,const dotprecision &a) throw()
{
   string str="";
   str << a;
   s << str;
   return s;
}

std::string & operator >>(std::string & s,dotprecision &a) throw()
{
   rndtype rnd;
   a_intg rndfl;

   if (ioflags.isset(IOFlags::rndup)) 
      rnd = RND_UP;
   else if (ioflags.isset(IOFlags::rnddown)) 
      rnd = RND_DOWN;
   else rnd = RND_NEXT;

   if (d_init_dm () == -1) 
   {
      // throw!
      //errmon (ERR_ALL(NOMOREMEMORY));
      //errmon (ERR_ALL(NOCONTINUEPOSSIBLE));
   }

   a= 0.0;   // AW! wg. Initialisierungsprobleme
   
   strcpy(dm,s.c_str());

   s = cxsc::d_scanp (a.akku,dm, rnd, &rndfl);
   if (rndfl) 
      ScanDotRndFlag = true;
  
   return s;
}
void operator >>(const std::string &s,dotprecision &a) throw()
{
   string s2(s);
   s2 >> a;
}
void operator >>(const char *s,dotprecision &a) throw()
{
   string s2(s);
   s2 >> a;
}

std::istream & operator >>(std::istream & s,dotprecision &a) throw()
{
   char c;
   string d="";

   skipeolnflag = inpdotflag = true;

   /* - skip white spaces ----------------------------------- */
   c = skipwhitespaces (s);

   /* - get sign and skip following white spaces ------------ */
   if (c == '+' || c == '-') 
   {
      d+=c;
      c = skipwhitespaces (s);
   }
   /* - skip leading zeros -------------------------------- */
   if (c == '0') 
      c = skipleadingchars (s, '0', '0');

   /* - get digits of integer part ----------------------- */
   do {
      if (c >= '0' && c <= '9') 
         d+= c;
      else 
         break;

      if (s.good()) 
         s.get(c);
      else 
         inpdotflag = false, 
         c = '\0';

   } while (s.good());

   /* - get point --------------------------------------- */
   if (c == '.') 
   {
      d+= '.';
      if (s.good()) 
         s.get(c);
      else 
         inpdotflag = false, 
         c = '\0';
   }
   /* - get digits of fractional part ------------------- */
   do {
      if (c >= '0' && c <= '9') 
         d+= c;
      else 
         break;

      if (s.good()) 
         s.get(c);
      else 
         inpdotflag = false, 
         c = '\0';

   } while (s.good());

   /* - get Exponent ------------------------------------ */
   if (c == 'E' || c == 'e')
   {
      d+= c;
      if (s.good()) 
         s.get(c);
      else inpdotflag = false, 
           c = '\0';

      /* - get sign of Exponent -------------------------- */
      if (c == '+' || c == '-') 
      {
         d+= c;
         if (s.good()) 
            s.get(c);
         else 
            inpdotflag = false, 
            c = '\0';
      }                                     
 
      /* - get Exponent digits --------------------------- */
      do {
         if (c >= '0' && c <= '9') 
            d+= c;
         else 
            break;

         if (s.good()) 
            s.get(c);
         else 
            inpdotflag = false, 
            c = '\0';

      } while (s.good());
   }

   /* Fehler auf zu langen Inputstring pruefen ---- mr???? */
   // Braucht bei Strings nicht mehr beruecksichtigt werden

   // --------------------------------------------------------------
                    // mindestens 1 Trennzeichen wird gefordert und uebergangen
   waseolnflag = (c == '\n');
        
   // --------------------------------------------------------------
   // erzeugten String scannen

   d>>a;

   return s;
}                                     
} // namespace cxsc

/****************************************************************/
/*                                                              */
/*      Filename        : d_outp.c                              */
/*                                                              */
/*      Entries         : void d_out                            */
/*                        (buffer,c,FormatFlag,FracDigits,      */
/*                        rnd,length)                           */
/*                        char *buffer;                         */
/*                        dotprecision c;                       */
/*                        a_intg FormatFlag,FracDigits,rnd;     */
/*                        a_intg *length                        */
/*                                                              */
/*      Arguments       :                                       */
/*                        buffer - output string                */
/*                        c - accu for output                   */
/*                        FormatFlag - format selection         */
/*                              -1 = scientific format          */
/*                               0 = fixed format               */
/*                             > 0 = variable format            */
/*                                   value is the total field   */
/*                                   width                      */
/*                        FracDigits - number of fraction digits*/
/*                        rnd - rounding monde (-1,0,1)         */
/*                        length - size of buffer string        */
/*                                                              */
/*                                                              */
/*      Description     : Decimal representation determined from*/
/*                        dotprecision akku                     */
/*                                                              */
/*      External        :                                       */
/*                        d_out  - conversion of akku           */
/*                                                              */
/*      Globals         : dm - I/O-buffer                       */
/*                                                              */
/*      Author          : M.Rauch                               */
/*      Date            : 1990-09-30                            */
/*                                                              */
/****************************************************************/



#include "dot_defs.hpp"
#include <stdlib.h>

namespace cxsc {

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* Stringbreich zur Aufnahme eines Akku                             */

/* char dm[A_DIGITS]; */
#if _WIN32
__declspec(thread) char *dm = NULL;
__declspec(thread) char *dmhlp = NULL;
#elif __APPLE__ && !CXSC_FORCE_TLS
char *dm = NULL;
char *dmhlp = NULL;
#else
__thread char *dm = NULL;
__thread char *dmhlp = NULL;
#endif

/* ---------------------------------------------------------------- */

int d_init_dm(void) 
{
 
   if (dm) 
      return 1;

   dmhlp = (char*) malloc (A_DIGITS);
   dm = (char*) malloc (A_DIGITS);

   return (dm && dmhlp) ? 0 : -1;
}

/* ---------------------------------------------------------------- */

void d_outp(char *buffer, Dotprecision c,
                  int FormatFlag, int FracDigits, int rnd,
                  int *length)

{
   a_intg dexpo,bdp,len;
   a_intg expo,i,digits,IntDigits,DecPlaces,vz;
   a_intg HoldWidth = (FracDigits < 0);
   char *s,*p;

#if TEST
   printf("d_outp\n");
#endif

   if (HoldWidth) FracDigits = -FracDigits;

/*                                                      */
/* Kehre n�tigenfalls Rundungsrichtung um               */
/*                                                      */

   if ((vz = (a_intg)c[A_SIGN]) != 0) {
     rnd = -rnd;
   }

   bdp = 2+A_I_DIGITS;
   len = bdp+1;
   if (c[A_END] > A_D_P) len += B_LENGTH * ((a_intg)c[A_END] - A_D_P);

   d_out (&dexpo, dm, &bdp, &len, c);
   dm[len] = '0';

/*                                                                      */
/* floating-point reprensentation                                       */
/*                                                                      */

   /* erzwinge n�tigenfalls Gleitkommadarstellung */

   if (FormatFlag > 0 && FormatFlag-FracDigits <= 2) {
     FormatFlag = FracDigits+3;
   }
   if (FormatFlag > 0) {
      if (dexpo < -((FracDigits+1)/2) ||
          (dexpo > 0 && dexpo >= FormatFlag-FracDigits-2)) {
        FormatFlag = -1;
        if (HoldWidth) {
          FracDigits = (FracDigits < 6) ? 0 : FracDigits-6;
        }
      }
   }

   if (FormatFlag == -1)
   {
     DecPlaces = FracDigits;
     digits = (len - bdp - 1) + (dexpo + 1);
     if (digits < DecPlaces) DecPlaces = digits-1;

     b_rnd (rnd, dm, digits+1, DecPlaces+1, &bdp, &dexpo);

     p = buffer; 
     *p++ = (vz) ? '-' : '+';           /* Vorzeichen */
     s = &dm[bdp-dexpo];
     *p++ = *s++;                               /* Digit vor Dezimalpunkt */
     if (FracDigits) {
       *p++ = '.';                              /* Dezimalpunkt */
       for (i=0;i<DecPlaces;i++) *p++ = *s++;   /* Nachkommadigits */
       for (;i<FracDigits;i++) *p++ = '0';      /* Nachkommadigits mit '0' */
     }
     *p++ = 'E';                                /* Exponentzeichen 'E' */
     *p++ = (dexpo<0) ? '-' : '+';              /* Vorzeichen Exponent */
     expo = (dexpo<0) ? -dexpo : dexpo;
     for (i=A_E_DIGITS-1; i >= 0; i--) {        /* Exponent */
       p[i] = expo%10 + '0';
       expo /= 10;
     }
     *length = 3+FracDigits+2+A_E_DIGITS;
   }

/*                                                                      */
/* fixed-point reprensentation                                          */
/*                                                                      */
 
   else
   {
     DecPlaces = (FracDigits < (len-bdp-1)) ? FracDigits : (len-bdp-1);
     if (dexpo >= 0) {
       IntDigits = dexpo+1;
     }
     else {
       IntDigits = 1;
       for (i=0; i < -dexpo; i++) dm[bdp+i] = '0';
       dexpo = 0;
     }
     digits = (len - bdp - 1) + (dexpo + 1);

     b_rnd (rnd, dm, digits+1, IntDigits+DecPlaces, &bdp, &dexpo);

     p = buffer; 
     *p++ = (vz) ? '-' : '+';                   /* Vorzeichen */
     s = &dm[bdp-IntDigits+1];
     for(i=0;i<IntDigits;i++) *p++ = *s++;      /* Digits vor Dezimalpunkt */
     if (FracDigits) {
       *p++ = '.';                              /* Dezimalpunkt */
       for (i=0;i<DecPlaces;i++) *p++ = *s++;   /* Nachkommadigits */
       for (;i<FracDigits;i++) *p++ = '0';      /* Nachkommadigits mit 0 */
     }
     *length = IntDigits+FracDigits+2;
   }
   if (FracDigits == 0)  (*length)--;
   buffer[*length] = '\0';
}

/****************************************************************/
/*                                                              */
/*      Filename        : d_out.c                               */
/*                                                              */
/*      Entries         : void d_out                            */
/*                        (dexpo,buffer,bdp,len,c)   4          */
/*                        a_intg *dexpo,*bdp,*len;              */
/*                        char *buffer;                         */
/*                        dotprecision c;                       */
/*                                                              */
/*      Arguments       :                                       */
/*                        dexpo - decimal exponent of first     */
/*                                digit                         */
/*                        buffer - output string                */
/*                        bdp - position of decimal point       */
/*                        len - input: usable size of buffer    */
/*                             output: totaly produced digits   */
/*                        c - accu for output                   */
/*                                                              */
/*                                                              */
/*      Description     : Decimal representation determined from*/
/*                        dotprecision akku                     */
/*                                                              */
/*      External        :                                       */
/*                        b_outf - conversion of fractional part*/
/*                        b_outi - conversion of integer part   */
/*                                                              */
/*      Globals         : b_cm__ - I/O-buffer                   */
/*                                                              */
/*      Author          : M.Rauch                               */
/*      Date            : 1990-09-30                            */
/****************************************************************/

}
// #include "o_defs.h"
#if WINDOWS_X86_32
#include "dot.hpp"
Dotprecision b_cm__;
#else
extern Dotprecision b_cm__;
#endif

namespace cxsc {
void d_out(a_intg *dexpo, char *buffer, a_intg *bdp, a_intg *len,
                 Dotprecision c)

{
   a_intg i,digits,cont;

#if TEST
   printf("d_out\n");
#endif

/* copy akku c to temporary akku b_cm__ */

   b_cm__[A_BEGIN] = c[A_BEGIN];
   b_cm__[A_END] = c[A_END];
   for (i=(a_intg)c[A_BEGIN]; i <= (a_intg)c[A_END]; i++) 
     b_cm__[i] = c[i];

/* test if akku is zero */

   if (b_cm__[A_BEGIN]==ZERO || b_cm__[A_END]==ZERO || b_cm__[A_BEGIN]>b_cm__[A_END]) {
     buffer[*bdp] = '0';
     for (i=*bdp+1; i < *len; i++) buffer[i] = '0';
     *dexpo = 0;
     return;
   }

   /* clear accu contents between number and decimal point      */
   for (i=(a_intg)b_cm__[A_END]+1; i <= A_D_P; i++) 
       b_cm__[i] = ZERO;
   for (i=A_D_P+1; i < (a_intg)b_cm__[A_BEGIN]; i++)
       b_cm__[i] = ZERO;

/*                                                                   */
/* conversion of integer part                                        */
/*                                                                   */

   *dexpo = -1;
   if (b_cm__[A_BEGIN] <= A_D_P) {
     digits = *len;
     b_outi(&digits,buffer,bdp,dexpo,b_cm__);
   }

/*                                                                   */
/* conversion of fraction part                                       */
/*                                                                   */

   digits = (*len < *bdp+1) ? 0 : *len - *bdp - 1;
   if (digits>0) {
     cont = 0;
     b_outf(&digits,buffer,bdp,&cont,b_cm__);
     if (*dexpo < 0) {
       for (cont=*bdp+1; cont < *len-1; cont++) {
         if (buffer[cont] != '0') break;
         (*dexpo)--;
       }
     }
   }

   return;
}

/****************************************************************/
/*                                                              */
/*      Filename        : d_scanp.c                             */
/*                                                              */
/*      Entries         : char* d_scanp (c,inpbuf,rnd,rndfl)    */
/*                        dotprecision c;                       */
/*                        char *inpbuf;                         */
/*                        a_intg rnd,*rndfl;                    */
/*                                                              */
/*      Arguments       :                                       */
/*                        c       - accu holding IEEE value     */
/*                        inpbuf  - buffer with input digits    */
/*                        rnd - direction of rounding           */
/*                        rndfl - flag 'rounding occurred'      */
/*                                                              */
/*      Description     : Convert a character string            */
/*                        to the corresponding IEEE value       */
/*                                                              */
/*      Author          : M.Rauch, University of Karlsruhe      */
/*      Date            : 1990-10-20                            */
/****************************************************************/


char* d_scanp(Dotprecision c, char *inpbuf, a_intg rnd, a_intg *rndfl)
{
   char *s;
   a_intg sign, dexpo, bdp, len;
   a_btyp *p,*pe;
/*                                                              */
/*      convert input string into a convertible form            */
/*                                                              */

   s = d_scan (inpbuf, &sign, &dexpo, dm, &bdp, &len);

/*                                                              */
/*      change rnd-direction if input string is negativ         */
/*                                                              */

   c[A_SIGN] = sign;
   if (sign && rnd) 
      rnd = -rnd;

/*                                                              */
/*      convert integer and fractional parts                    */
/*                                                              */

   d_scani (c, dm, &dexpo, &bdp, &len);
   *rndfl = d_scanf (c, dm, &dexpo, &bdp, &len, rnd);

/*                                                              */
/*      adjust akku and test if akku is ZERO                    */
/*                                                              */

   for (p=&c[(a_intg)c[A_BEGIN]],pe=&c[(a_intg)c[A_END]]; p <= pe; p++)
   {
      if (*p != ZERO) 
         break;
      c[A_BEGIN]++;
   }
   for (; pe >= p; pe--)
   {
      if (*pe != ZERO) 
         break;
      c[A_END]--;
   }
   if (p > pe) 
   {
      c[A_BEGIN] = c[A_END] = 0;
   }

   c[A_STATUS] |= A_PZERO+A_MZERO;
   return s;
}

/****************************************************************/
/*                                                              */
/*      Filename        : d_scani.c                             */
/*                                                              */
/*      Entries         : d_scani (c, buffer, dexpo, bdp, len)  */
/*                        Dotprecision c;                       */
/*                        char *buffer;                         */
/*                        int *dexpo,*bdp,*len;                 */
/*                                                              */
/*      Arguments       :                                       */
/*                        c       - accu holding IEEE value     */
/*                        buffer  - buffer with input digits    */
/*                        bdp - position of decimal point       */
/*                        dexpo   - exponent of first non-zero  */
/*                                  digit                       */
/*                        len - position behind last input digit*/
/*                                                              */
/*      Description     : Convert a character string    */
/*        to the integer part of IEEE value     */
/*                                                              */
/*      Author          : M.Rauch, University of Karlsruhe      */
/*      Date            : 1990-10-19                            */
/****************************************************************/


void d_scani(Dotprecision c, char *buffer, a_intg *dexpo,
                   a_intg *bdp, a_intg *len)

{
   a_intg i,j;
   a_btyp carry;
   a_btyp *s,*p,h,hh;
   char *q,*qe;

/*                                                                   */
/* initialize                                                        */
/*                                                                   */
   c[A_BEGIN] = c[A_END] = A_D_P,
   c[A_D_P] = ZERO;
   if (*dexpo < 0) 
      return;

   i = (*dexpo+1) % B2D_LOG10;
   if (i) 
   {
      q = buffer + *bdp - *dexpo - 1;
      for (; i < B2D_LOG10; i++,(*dexpo)++,q--) 
         *q = '0';
   }
   for (i=*len; i <= *bdp; i++) 
      buffer[i] = '0';

   q = buffer + *bdp - *dexpo;
   qe = buffer + *bdp;
   s = &c[(a_intg)c[A_BEGIN]];

/*                                                                   */
/* convert by repeated multiplication                                */
/*                                                                   */

   while (q < qe)
   {
      /* get decimal digits                                */
      for (j=0,i=B2D_LOG10;i>0;i--,q++) 
      {
         /* carry = 10*carry + *q - '0'; */
         j = j*10 + *q - '0';
      }
      carry = j;

      for (p=&c[A_D_P]; p >= s; p--)
      {
         hh = GETLOW(*p)*B2D_POWER+carry,
         h = GETHIGH(*p)*B2D_POWER+GETHIGH(hh),
         carry = GETHIGH(h),
         *p = MOVEHIGH(h) | GETLOW(hh);
      }

      /* adjust digits                               */
      if (carry) 
      {
         c[A_BEGIN]--, s--;
         *s = carry;
      }

   }
   return;
}

/****************************************************************/
/*                                                              */
/*      Filename        : d_scanf.c                             */
/*                                                              */
/*      Entries         : a_intg d_scanf (c,buffer,dexpo,bdp,len,rnd) */
/*                        Dotprecision c;                       */
/*                        char *buffer;                         */
/*                        a_intg *dexpo,*bdp,*len,rnd;          */
/*                                                              */
/*      Arguments       :                                       */
/*                        c       - accu holding IEEE value     */
/*                        buffer  - buffer with input digits    */
/*                        bdp - position of decimal point       */
/*                        dexpo   - exponent of first non-zero  */
/*                                  digit                       */
/*                        len - position behind last input digit*/
/*                        rnd - direction of rounding           */
/*                                                              */
/*      Description     : Convert a character string            */
/*            to the fractional part of IEEE value              */
/*                                                              */
/*      Author          : M.Rauch, University of Karlsruhe      */
/*      Date            : 1990-10-19                            */
/****************************************************************/


a_intg d_scanf(Dotprecision c, char *buffer, a_intg *dexpo,
                     a_intg *bdp, a_intg *len, a_intg rnd)

{
   a_intg i,j,rndflag = 0;
   a_btyp mod, carry;
   a_btyp *s,*p,*pe,h,hh;
   char *q,*qe;


/*                                                                   */
/* initialize                                                        */
/*                                                                   */

   if (*len+1 <= *bdp) 
      return rndflag;

   for (i=*bdp-*dexpo; *dexpo < 0; (*dexpo)++,i--) 
      buffer[i] = '0';

   i = (*len - *bdp - 1) % B2D_LOG10;
   if (i) 
   {
      q = buffer + *len;
      for (; i < B2D_LOG10; i++,(*len)++,q++) 
         *q = '0';
   }

   qe = buffer + *bdp + 1;
   q = buffer + *len;
   s = &c[(a_intg)c[A_END]];

/*                                                                   */
/* convert by repeated division                                      */
/*                                                                   */

   carry = 0;
   while (q > qe)
   {
      /* get decimal digits                                */
      for (q-=B2D_LOG10,j=i=0;i<B2D_LOG10;i++) 
      {
         /* mod = 10*mod + q[i] - '0'; */
         j = j*10 + q[i] - '0';
      }
      mod = j + carry;
      if (mod == B2D_POWER) 
      {
         carry = 1, mod = 0;
      } else 
         carry = 0;

      p=&c[A_D_P+1];

      do {
         for (;p<=s;p++)
         {
            h = GETHIGH(*p) | MOVEHIGH(mod),
            hh = GETLOW(*p) | MOVEHIGH(h%B2D_POWER),
            mod = hh%B2D_POWER,
            *p = MOVEHIGH(h/B2D_POWER) | (hh/B2D_POWER);
         }
         if (mod && c[A_END] < A_LENGTH-1)
         {
            c[A_END]++, s++;
            *s = ZERO;
         }
      } while (p <= s);

      if (mod) 
         rndflag = 1;

      if (rnd < 0) 
         mod = 0;
      else if (rnd == 0 && mod < B2D_POWER/2) 
         mod = 0;

      if (mod) 
      {
         /* Runden nach oben : c um 1 Inkermentieren */

         p  = &c[A_LENGTH-1],        /* == s */
         pe = &c[A_D_P+1];
         for(; p >= pe; p--) 
         {
            (*p)++;
            if (*p) 
               break;
         }
         if (p < pe) 
            carry = 1;

      }
   }

   if (carry)
   {
      p = &c[A_D_P];
      pe = &c[(a_intg)c[A_BEGIN]];
      for(; p >= pe; p--)
      {
         (*p)++;
         if (*p) 
            break;
      }
      if (p < pe) 
      {
         --c[A_BEGIN];
         *p = 1;
      }
   }

   return rndflag;
}

/****************************************************************/
/*                                                              */
/*      Filename        : d_scan.c                              */
/*                                                              */
/*      Entries         : char* d_scan.c                        */
/*                        (inpbuf, sign,dexpo,outbuf,bdp,len)   */
/*                        a_intg *sign,*dexpo,*bdp,*len;        */
/*                        char *inpbuf, *outbuf;                */
/*                                                              */
/*      Arguments       :                                       */
/*                        inpbuf - string to scan               */
/*      Returnvalues    :                                       */
/*                        sign   - sign flag                    */
/*                        dexpo  - decimal exponent of first    */
/*                                 digit                        */
/*                        outbuf - output string                */
/*                                 min.length: A_DIGITS         */
/*                        bdp - position of decimal point       */
/*                        len - totaly produced digits          */
/*                             output: totaly produced digits   */
/*                        c - accu for output                   */
/*                                                              */
/*                                                              */
/*      Description     : the input string is scanned and       */
/*                        transferred into an ordered form      */
/*                                                              */
/*      Author          : M.Rauch                               */
/*      Date            : 1990-10-14                            */
/****************************************************************/


char* d_scan (char *inpbuf, a_intg *sign, a_intg *dexpo,
                    char *outbuf, a_intg *bdp, a_intg *len)

{
   a_intg  start, point, mantend, expo, end;
   a_intg  digits,fl,i,j;
   char c;

   /* intialiaze parameters                                             */

   *bdp = 1 + A_I_DIGITS;
   *len = (*bdp)+1;
   *dexpo = 0;

   /*                                                                   */
   /* determine structure of input string                               */
   /*                                                                   */

   fl = 0;
   for (start=0; (c = inpbuf[start]) != 0; start++) 
   { /* skip white spaces */
      if (c > ' ') 
         break;
   }

   if (c == '-' || c == '+')                          /* determine sign */
      *sign = ((c == '-') ? 1 : 0), 
      start++;
   else 
      *sign = 0;

   for (; (c = inpbuf[start]) != 0; start++) 
   {        /* skip white spaces */
      if (c > ' ') 
         break;
   }
   for (; (c = inpbuf[start]) != 0; start++) 
   {        /* skip leading zero */
      if (c != '0') 
         break;
   }

   for (end=start; (c = inpbuf[end]) != 0; end++) 
   {   /* skip intdigits */
      if (c < '0' || c > '9') 
         break;
   }
   if (c == '.')                                  /* pos. of decimalpoint */
      point = end, 
      end++, 
      fl++; 
   else 
      point = -1;
   
   for (; (c = inpbuf[end]) != 0; end++) 
   {            /* skip fracdigits */
      if (c < '0' || c > '9') 
         break;
   }
   mantend = end;                                     /* end of mantissa */

   if (c == 'E' || c == 'e')                       /* determine exponent */
   {
      c = inpbuf[++end];
      if (c == '+' || c == '-')                        /* sign of exponent */
         expo = ((c == '-') ? 1 : 0), 
         end++;
      else 
         expo = 0;
      for (i=0; (c = inpbuf[end]) != 0; end++) 
      {       /* value of exponent */
         if (c < '0' || c > '9') 
            break;
         if (i >= (A_E_MAX/10)) 
         {
            break;     /* Error Exponent To big ---- mr???? */
         }
         i = 10*i + c - '0';
      }
      expo = expo ? -i : i;
   } else 
   expo = 0;

   if (c) 
      end++;                    /* skip at least one delimeter char */

   /*                                                                   */
   /* check if there are some digits in the mantissa                    */
   /*                                                                   */

   if (start+fl == mantend) 
   {
      outbuf[*bdp] = '0';
      return &inpbuf[end];
   }

   if (point != -1) 
   {
      digits = mantend - start - 1;
      *dexpo = (point - start - 1) + expo;
      if (start == point) digits++, 
      (*dexpo)++;
   } else 
   {
      digits = mantend - start;
      *dexpo = (mantend - start - 1) + expo;
   }

   if (*dexpo >= (A_I_DIGITS-10)) 
   {
      /* Error Input Value Too big ---- mr ???? */
   }
   if (*dexpo <= - (A_F_DIGITS-10)) 
   {
      /* Error Input Value Too small ---- mr ???? */
   }
    
   /*                                                                   */
   /* take known digits into output string                              */
   /*                                                                   */

   *len = *bdp - *dexpo + digits;
   if (*len >= A_DIGITS) 
   {
      /* Error Input String Too long ---- mr ???? */
   }
   for (j=*len-1,i=mantend-1; i >= start; i--) 
   {
      if ((c = inpbuf[i]) != '.') 
         outbuf[j--] = c;
   }
   if (start == point) 
      outbuf[j] = '0';

   /*                                                                   */
   /* be shure that the decimal point is in the output string           */
   /*                                                                   */

   if (*dexpo < 0) 
   {
      for (j=-*dexpo,i=*bdp-*dexpo-1; j > 0; j--,i--) 
         outbuf[i] = '0';
      *dexpo = 0;
   }

   if (*len <= *bdp) 
   {
      for (j=(*bdp)-(*len)+1,i=*len; j > 0; j--,i++) 
         outbuf[i] = '0';
      *len = (*bdp)+1;  
   }

   return &inpbuf[end];
}

} // namespace cxsc

