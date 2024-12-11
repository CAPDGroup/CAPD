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

/* CVS $Id: ioflags.cpp,v 1.26 2014/01/30 17:23:45 cxsc Exp $ */

// Modul ioflags.cpp
//
//      Hier wird die Formatierung der Ein-/Ausgabe ueber Streams
//      verwaltet, hierfuer gibt es folgende Moeglichkeiten :
//
//      Zahlenformat der Ausgabe:
//        Hex           - Alle Zahlen werden in Hexdarstellung ausgegeben;
//                        hierbei enstehen keine Rundungsfehler, da die
//                        interne Zahlendarstellung binaer ist
//        Dec           - Alle Zahlen werden in Dezimalschreibweise ausgegeben;
//                        da hierbei normalerweise gerundet werden muss kann
//                        eine spezielle Rundung gewaehlt werden
//
//      Rundungsoptiionen:
//        RndNone       - keine spezifische Rundung (Standardeinstellung),
//                        es wird die im C-Paket eingebaute Standardausgabe
//                        verwendet
//        RndUp         - es wird auf 1 Ulp genau nach unten gerundet
//        RndDown       - es wird auf 1 Ulp genau nach oben gerundet
//        RndNext       - es wird auf 1 Ulp genau nach oben/unten gerundet
//
//      Ausgabefeldaufbau:
//        SetPrecision  - Breite des Ausgabefeldes, Anzahl der Nachkomma-
//                        stellen
//        GetPrecision  - abfragen der Breite des Ausgabefeldes oder
//                        der Anzahl der Nachlommastellen
//
//      Feldaufbau:
//        RightJust
//            Die Ausgabe in dezimaler Darstellung erfolgt
//            im aktuellen Format rechtsb�ndig. Werden f�r die
//            Ausgabe der Zahl weniger Zeichen ben�tigt als der
//            aktuelle Wert von width angibt, so werden
//            Leerzeichen vor der Zahl eingef�gt. width bezeichnet dabei
//            f�r reelle Gleitpunktzahlen den zuletzt an SetPrecision(), f�r
//            Dotprecision-Akkumulatoren den zuletzt an SetDotPrecision()
//            �bergebenen Parameter. Dieser Wert ist voreingestellt.
//       
//        LeftJust
//            Die Ausgabe in dezimaler Darstellung erfolgt im aktuellen
//            Format linksb�ndig. Werden f�r die Ausgabe der Zahl weniger
//            Zeichen ben�tigt als der aktuelle Wert von width angibt, so
//            werden Leerzeichen nach der Zahl eingef�gt.
//       
//       
//        Blank
//            Die Ausgabe einer positiven Zahl in dezimaler Darstellung
//            beginnt mit einem f�hrenden Leerzeichen ' '.
//       
//        ForceSign
//            Die Ausgabe einer positiven Zahl in dezimaler Darstellung
//            beginnt mit einem f�hrenden Pluszeichen '+'.
//       
//        NoBlank
//            Die Ausgabe einer positiven Zahl in dezimaler Darstellung
//            beginnt ohne ein f�hrendes Leerzeichen. Dieser Wert ist
//            voreingestellt.
//       
//        Fixed
//            Die Ausgabe in dezimaler Darstellung erfolgt in
//            Fixpunktnotation:
//       
//             [-]ddd.dddddd
//       
//            in einem Feld der Mindestbreite width. Es werden digits
//            Nachkommastellen ausgegeben. Ist digits gleich 0, so wird kein
//            Dezimalpunkt ausgegeben. width und digits bezeichnen dabei f�r
//            reelle Gleitpunktzahlen den zuletzt an SetPrecision(), f�r
//            Dotprecision-Akkumulatoren den zuletzt an SetDotPrecision()
//            �bergebenen Parameter.
//       
//        Scientific
//            Die Ausgabe in dezimaler Darstellung erfolgt in
//            wissenschaftlicher Notation
//       
//             [-]d.ddddddE+-eee[e]
//       
//            in einem Feld der Mindestbreite width.
//            Es werden eine Vorkomma- und digits Nachkommastellen ausgegeben.
//            Ist digits gleich 0, so entf�llt der Dezimalpunkt. F�r reelle
//            Gleitpunktzahlen werden drei, f�r Dotprecision-Akkumulatoren
//            vier Exponentenstellen ausgegeben.
//       
//        Variable
//            Die Ausgabe in dezimaler Darstellung erfolgt in Abh�ngigkeit
//            des auszugebenden Wertes in wissenschaftlicher oder in 
//            Fixpunktnotation. 
//            Die wissenschaftliche Notation wird gew�hlt, falls der Exponent 
//            der Zahl kleiner als -(digts+1)/2 oder gr��er als die Anzahl 
//            der darstellbaren Vorkommastellen (width-digits-'.'-'+/-/ ')
//            ist. Andernfalls wird in Fixpunktnotation ausgegeben. 
//       
//        VarFixWidth
//            Wie Variable, mit dem Unterschied, da� bei Umschaltung
//            auf Scientific die Anzahl der digits um bis zu 5 bzw. 6 Stellen
//            verringert wird, damit der Exponent in die �bergebenene
//            Feldbreite pa�t
//       
//      Diverses:
//        SaveOpt       - sichern der Ein-/Ausgabeformatierung auf einem
//                        eigenen Stack
//        RestoreOpt    - regenerieren der zuletzt ueber SaveOpt gesicherten
//                        Ein-/Ausgabeformatierung von einem eigenen Stack
//----------------------------------------------------------------------------
//
//      Aufruf/Verwendung der angegebenen Optionen :
//
//      - alle Optionen (ausser GetPrecision) koennen direkt mittels den
//        Ein-/Ausgabeoperatoren >> und << auf istream/ostream's (und string's)
//        angewendet werden
//      - GetPrecision hat einen Parameter (WIDTH oder DIGIT) und liefert
//        als int den entsprechenden Wert zurueck
//      - SetPrecision hat zwei Parameter (width, digit), hat einer dieser
//        Parameter einen Wert != CXSC_UNDEFINED (= -1) so wird entsprechend
//        realwidth und realdigits gesetzt, dabei haben realwidth und digits
//        folgende Eigenschaften :
// ????   realwidth  = 0     Feldbreite wird dynamisch vergeben
// ????              > 0     Feldbreite betraegt mindestens realwidth Zeichen
// ????   realdigits = 0     Gleitkommaausgabe
// ????              > 0     fixe Anzahl Nachkommastellen
//----------------------------------------------------------------------------
#include <string>

#include "ioflags.hpp"

namespace cxsc {

#define IOFLAGS_TYPMASK          0x00000000FLL
#define IOFLAGS_HEX              0x000000001LL
#define IOFLAGS_DEC              0x000000002LL

#define IOFLAGS_RNDMASK          0x0000000F0LL
#define IOFLAGS_RNDNEXT          0x000000010LL
#define IOFLAGS_RNDUP            0x000000020LL
#define IOFLAGS_RNDDOWN          0x000000040LL
#define IOFLAGS_RNDNONE          0x000000080LL

#define IOFLAGS_FORMATMASK       0x00000FF00LL
#define IOFLAGS_REALFORMAT       0x000000100LL
#define IOFLAGS_DOTFORMAT        0x000000200LL

#define IOFLAGS_ADJUSTMASK       0x000030000LL
#define IOFLAGS_RIGHTJUST        0x000010000LL
#define IOFLAGS_LEFTJUST         0x000020000LL

#define IOFLAGS_BLANKMASK        0x000F00000LL
#define IOFLAGS_BLANK            0x000100000LL
#define IOFLAGS_FORCESIGN        0x000200000LL
#define IOFLAGS_NOBLANK          0x000400000LL

#define IOFLAGS_FORMMODEMASK     0x00F000000LL
#define IOFLAGS_VARIABLE         0x001000000LL
#define IOFLAGS_FIXED            0x002000000LL
#define IOFLAGS_SCIENTIFIC       0x004000000LL
#define IOFLAGS_VARFIXWIDTH      0x008000000LL

#define IOFLAGS_SAVEOPT          0x080000000LL
#define IOFLAGS_RESTOREOPT       0x040000000LL

//Sparse output flags
#define IOFLAGS_INOUTMASK        0xF00000000LL
#define IOFLAGS_FULL_INOUT       0x100000000LL
#define IOFLAGS_SPARSE_INOUT     0x200000000LL
#define IOFLAGS_MATRIX_MARKET    0x400000000LL



const unsigned long long int IOFlags::dummy   = 0x000000000LL;
const unsigned long long int IOFlags::allmask = 0xFFFFFFFFFLL;

const unsigned long long int IOFlags::typmask = IOFLAGS_TYPMASK,
               IOFlags::hex     = IOFLAGS_HEX,
               IOFlags::dec     = IOFLAGS_DEC;

const unsigned long long int IOFlags::rndmask = IOFLAGS_RNDMASK,
               IOFlags::rndnext = IOFLAGS_RNDNEXT,
               IOFlags::rndup   = IOFLAGS_RNDUP,
               IOFlags::rnddown = IOFLAGS_RNDDOWN,
               IOFlags::rndnone = IOFLAGS_RNDNONE;

const unsigned long long int IOFlags::realformat = IOFLAGS_REALFORMAT,
               IOFlags::dotformat  = IOFLAGS_DOTFORMAT;
               
const unsigned long long int IOFlags::rightjust  = IOFLAGS_RIGHTJUST,
               IOFlags::leftjust   = IOFLAGS_LEFTJUST;

const unsigned long long int IOFlags::blank      = IOFLAGS_BLANK,
               IOFlags::forcesign  = IOFLAGS_FORCESIGN,
               IOFlags::noblank    = IOFLAGS_NOBLANK;

const unsigned long long int IOFlags::variable    = IOFLAGS_VARIABLE, 
               IOFlags::varfixwidth = IOFLAGS_VARFIXWIDTH,
               IOFlags::fixed       = IOFLAGS_FIXED,
               IOFlags::scientific  = IOFLAGS_SCIENTIFIC;
               
const unsigned long long int IOFlags::saveopt    = IOFLAGS_SAVEOPT,
               IOFlags::restoreopt = IOFLAGS_RESTOREOPT;

const unsigned long long int IOFlags::inoutmask = IOFLAGS_INOUTMASK,
               IOFlags::fullinout    = IOFLAGS_FULL_INOUT,
               IOFlags::sparseinout = IOFLAGS_SPARSE_INOUT,
               IOFlags::matrixmarketinout = IOFLAGS_MATRIX_MARKET;               


const IOFlags Hex         ( IOFLAGS_HEX );
const IOFlags Dec         ( IOFLAGS_DEC );
const IOFlags RndNext     ( IOFLAGS_RNDNEXT );
const IOFlags RndUp       ( IOFLAGS_RNDUP );
const IOFlags RndDown     ( IOFLAGS_RNDDOWN );
const IOFlags RndNone     ( IOFLAGS_RNDNONE );
const IOFlags RealFormat  ( IOFLAGS_REALFORMAT );
const IOFlags DotFormat   ( IOFLAGS_DOTFORMAT );
const IOFlags RightJust   ( IOFLAGS_RIGHTJUST );
const IOFlags LeftJust    ( IOFLAGS_LEFTJUST );
const IOFlags Blank       ( IOFLAGS_BLANK );
const IOFlags ForceSign   ( IOFLAGS_FORCESIGN );
const IOFlags NoBlank     ( IOFLAGS_NOBLANK );
const IOFlags Fixed       ( IOFLAGS_FIXED );
const IOFlags Scientific  ( IOFLAGS_SCIENTIFIC );
const IOFlags Variable    ( IOFLAGS_VARIABLE );
const IOFlags VarFixWidth ( IOFLAGS_VARFIXWIDTH );
const IOFlags SaveOpt     ( IOFLAGS_SAVEOPT );
const IOFlags RestoreOpt  ( IOFLAGS_RESTOREOPT );
const IOFlags FullInOut   ( IOFLAGS_FULL_INOUT );
const IOFlags SparseInOut ( IOFLAGS_SPARSE_INOUT );
const IOFlags MatrixMarketInOut ( IOFLAGS_MATRIX_MARKET );


                                                             
//----------------------------------------------------------------------------

IOFlags ioflags (
  IOFLAGS_DEC | IOFLAGS_RNDNEXT | IOFLAGS_DOTFORMAT |
  IOFLAGS_RIGHTJUST | IOFLAGS_NOBLANK | IOFLAGS_VARIABLE | IOFLAGS_FULL_INOUT);
               
//----------------------------------------------------------------------------

int    realwidth   = 10;
int    realdigits  = 6;
int    dotwidth    = 15;
int    dotdigits   = 10;

//----------------------------------------------------------------------------
// ScanDotRndFlag: zeigt an, ob beim Einlesen eines dot's gerundet wurde.
//                 Kann vom Programmierer vor ein Inputoperation auf
//                 FALSE gesetzt werden, und nach dem Input ausgewertet
//                 werden.
bool    ScanDotRndFlag = false;
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------

//----------------------------------------------------------------------------

class OptionList {
  friend void _SaveOpt();
  friend void _RestoreOpt();

  IOFlags opt;
  int realwidth, realdigits;
  int dotwidth, dotdigits;
  void* next;
};
static void* LastOption = NULL;

void _SaveOpt ()
{
  OptionList *tmp = new OptionList;

  tmp->next = LastOption;
  tmp->opt  = ioflags;
  tmp->realwidth = realwidth;
  tmp->realdigits = realdigits;
  tmp->dotwidth = dotwidth;
  tmp->dotdigits = dotdigits;
  LastOption = (void *)tmp;
}

void _RestoreOpt ()
{
  OptionList *tmp = (OptionList *) LastOption;

  if (tmp == NULL) return;
  dotdigits = tmp->dotdigits;
  dotwidth = tmp->dotwidth;
  realdigits = tmp->realdigits;
  realwidth = tmp->realwidth;
  ioflags = tmp->opt;
  LastOption = (void *) (tmp->next);
  delete tmp;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void SetOption (const IOFlags &r)
{
  // IOFlags r = rin;
  int i;

  if (r.isset(IOFlags::saveopt))    _SaveOpt ();
  if (r.isset(IOFlags::restoreopt)) _RestoreOpt ();

  if (r.isset(IOFlags::realformat)) {
    i = r.getprecision(WIDTH);  if (i != CXSC_UNDEFINED) realwidth = i;
    i = r.getprecision(DIGITS); if (i != CXSC_UNDEFINED) realdigits = i;
  }
  if (r.isset(IOFlags::dotformat)) {
    i = r.getprecision(WIDTH);  if (i != CXSC_UNDEFINED) dotwidth = i;
    i = r.getprecision(DIGITS); if (i != CXSC_UNDEFINED) dotdigits = i;
  }

  if (r.isset(IOFlags::dec))
    ioflags.setf (IOFlags::dec, IOFlags::typmask);
  else if (r.isset(IOFlags::hex))
    ioflags.setf (IOFlags::hex, IOFlags::typmask);

  if (r.isset(IOFlags::rndnext))
    ioflags.setf (IOFlags::rndnext, IOFlags::rndmask);
  else if (r.isset(IOFlags::rndup))
    ioflags.setf (IOFlags::rndup, IOFlags::rndmask);
  else if (r.isset(IOFlags::rnddown))
    ioflags.setf (IOFlags::rnddown, IOFlags::rndmask);
  else if (r.isset(IOFlags::rndnone))
    ioflags.setf (IOFlags::rndnone, IOFlags::rndmask);

  if (r.isset(IOFlags::rightjust))
    ioflags.setf (IOFlags::rightjust, IOFLAGS_ADJUSTMASK);
  else if (r.isset(IOFlags::leftjust))
    ioflags.setf (IOFlags::leftjust, IOFLAGS_ADJUSTMASK);

  if (r.isset(IOFlags::blank))
    ioflags.setf (IOFlags::blank, IOFLAGS_BLANKMASK);
  else if (r.isset(IOFlags::forcesign))
    ioflags.setf (IOFlags::forcesign, IOFLAGS_BLANKMASK);
  else if (r.isset(IOFlags::noblank))
    ioflags.setf (IOFlags::noblank, IOFLAGS_BLANKMASK);
  
  if (r.isset(IOFlags::variable))
    ioflags.setf (IOFlags::variable, IOFLAGS_FORMMODEMASK);
  else if (r.isset(IOFlags::varfixwidth))
    ioflags.setf (IOFlags::varfixwidth, IOFLAGS_FORMMODEMASK);
  else if (r.isset(IOFlags::fixed))
    ioflags.setf (IOFlags::fixed, IOFLAGS_FORMMODEMASK);
  else if (r.isset(IOFlags::scientific))
    ioflags.setf (IOFlags::scientific, IOFLAGS_FORMMODEMASK);

  if (r.isset(IOFlags::fullinout))
    ioflags.setf (IOFlags::fullinout, IOFlags::inoutmask);
  else if (r.isset(IOFlags::sparseinout))
    ioflags.setf (IOFlags::sparseinout, IOFlags::inoutmask);
  else if (r.isset(IOFlags::matrixmarketinout))
    ioflags.setf (IOFlags::matrixmarketinout, IOFlags::inoutmask);

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
  
IOFlags::IOFlags() {
  flagmask = 0; width = digits = CXSC_UNDEFINED;
}
IOFlags::IOFlags(const unsigned long long int flag) {
  flagmask = flag; width = digits = CXSC_UNDEFINED;
}
IOFlags::IOFlags(const IOFlags &r) {
  flagmask = r.flagmask, width = r.width, digits = r.digits;
}
  
void IOFlags::setf(const unsigned long long int flag) {
  flagmask = flag;
}
void IOFlags::setf(const unsigned long long int flag, const unsigned long long int mask) {
  flagmask = (flagmask & ~mask) | flag;
}

unsigned long long int IOFlags::getf(const unsigned long long int mask) const {
  return flagmask & mask;
} 

int IOFlags::isset(const unsigned long long int  flag) const {
  return ((flagmask & flag) == flag);
}
  
void IOFlags::setprecision (const int value, const precisionflag sel)
{
  switch (sel) {
    case WIDTH : width  = value; break;
    case DIGITS: digits = value; break;
  }
}


int  IOFlags::getprecision (const precisionflag sel) const
{ 
  switch (sel) {
    case WIDTH : return width;
    case DIGITS: return digits;
  }
  return CXSC_UNDEFINED;
}  

//----------------------------------------------------------------------------

IOFlags SetPrecision (const int width, const int digits)
{
  IOFlags tmp = RealFormat;
  tmp.setprecision (width, WIDTH);
  tmp.setprecision (digits, DIGITS);
  return tmp;
}
void _SetPrecision (const int width, const int digits)
{
  SetOption (SetPrecision(width,digits));
}

int GetPrecision (const precisionflag sel)
{
  switch (sel) {
    case WIDTH : return realwidth;
    case DIGITS: return realdigits;
  }
  return CXSC_UNDEFINED;
}

//----------------------------------------------------------------------------

IOFlags SetDotPrecision (const int width, const int digits)
{
  IOFlags tmp = DotFormat;
  tmp.setprecision (width, WIDTH);
  tmp.setprecision (digits, DIGITS);
  return tmp;
}
void _SetDotPrecision (const int width, const int digits)
{
  SetOption (SetDotPrecision(width,digits));
}

int GetDotPrecision (const precisionflag sel)
{
  switch (sel) {
    case WIDTH : return dotwidth;
    case DIGITS: return dotdigits;
  }
  return CXSC_UNDEFINED;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

std::ostream& operator<< (std::ostream& s, const IOFlags &r)
{
  SetOption (r);
  return s;
}

std::istream& operator>> (std::istream& s, const IOFlags &r)
{
  SetOption (r);
  return s;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

std::string & operator<< (std::string &s, const IOFlags &r)
{
  SetOption (r);
  return s;
}

std::string & operator>> (std::string &s, const IOFlags &r)
{
  SetOption (r);
  return s;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

bool inpdotflag;
bool skipeolnflag;
bool waseolnflag;

//--------------------------------------------------------------------

char skipwhitespaces (std::istream& s)
{
  char c;
  while (s.good()) {
    s.get(c);
    if (!c || c > ' ' || (c == '\n' && !skipeolnflag)) return c;
  }                  
  inpdotflag = false;
  return '\0';
}

// besser:
std::istream & eatwhite(std::istream &is) // B.Stroustrup: C++ 21.3
{
   char c;
   while(is.get(c))
   {
      if(!isspace(c))
      {
         is.putback(c);
         break;
      }
   }
   return is;
}

//--------------------------------------------------------------------

char skipwhitespacessinglechar (std::istream& s, char c)
{
  char ch;
  ch = skipwhitespaces (s);
  if (ch == c) {
    if (s.good()) s.get(ch);
    else inpdotflag = false, ch = '\0';
  }
  if (ch && ch <= ' ') ch = skipwhitespaces (s);

  return ch;
}
    
//--------------------------------------------------------------------

char skipleadingchars (std::istream& s, char beg, char end)
{
  char c;
  while (s.good()) {
    s.get(c);
    if (!c || c < beg || c > end) return c;
  }
  inpdotflag = false;
  return '\0';
}           
//--------------------------------------------------------------------

std::string skipwhitespaces (std::string s)
{
  while (s.length() && s[0] <= ' ' && (s[0] != '\n' || skipeolnflag))
     s.erase(0,1);
  return s;
}
 //--------------------------------------------------------------------

std::string skipwhitespacessinglechar (std::string s, char c)
{
  s = skipwhitespaces (s);
  if (s[0] == c) s.erase(0,1);
  s = skipwhitespaces (s);
  return s;
}

//--------------------------------------------------------------------

std::string skipleadingchars (std::string s, char beg, char end)
{
  while (s.length() && s[0] >= beg && s[0] <= end) s.erase(0,1);
  return s;
}
char* cskipwhitespaces (char *s)
{
  while (*s && *s <= ' ' && (*s != '\n' || skipeolnflag)) s++;
  return s;
}

//--------------------------------------------------------------------

char* cskipwhitespacessinglechar (char *s, char c)
{
  s = cskipwhitespaces (s);
  if (*s == c) s++;
  s = cskipwhitespaces (s);
  return s;
}

//--------------------------------------------------------------------

char* cskipleadingchars (char *s, char beg, char end)
{
  while (*s && *s >= beg && *s <= end) s++;
  return s;
}              
//--------------------------------------------------------------------

//--------------------------------------------------------------------

void skipwhitespacesandputback (std::istream& s, char c) {
  skipeolnflag = inpdotflag = true;
  c = skipwhitespacessinglechar (s, '(');
  if (inpdotflag) s.putback(c);
}

//--------------------------------------------------------------------
                                   
} // namespace cxsc

