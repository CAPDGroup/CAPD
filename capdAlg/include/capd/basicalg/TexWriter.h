//////////////////////////////////////////////////////////////////////////////
///
///  @file TexWriter.h
///  
///  @author kapela  @date   Apr 17, 2011
//////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2011 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_BASICALG_TEXWRITER_H_
#define _CAPD_BASICALG_TEXWRITER_H_

#include "capd/intervals/Interval.h"
#include "capd/rounding/RoundingTraits.h"
#include "capd/vectalg/Vector.h"
#include <complex>
#include <deque>
#include <iostream>
#include <cassert>
#include "capd/intervals/lib.h"
#ifdef __USE_FILIB__
#include "capd/filib/Interval.h"
#endif

namespace capd {

using capd::rounding::setRounding;

template<typename T>
std::string printToString(const T & data, int precision, bool fixed = false){
  std::ostringstream buf;
  buf.precision(precision);
  if(fixed) 
    buf << std::fixed;
  buf << data;
  return buf.str();
}

//-------------------------------------------------------------------------------

class TexWriter { //: public std::ostream

protected:
  std::ostream & out;
public:
  enum FloatStyle { FloatSci, FloatFix, FloatFix2 , FloatTxt};
  enum ImStyle {ImExplicit,  // i symbol is given explicitely
                ImSymbolii   // i symbol is defined as command \ii
  }; 
  enum ExponentStyle {  
    CStyle,        // C convention for exponents e.g 1.0e-5
    TexStyle       // TeX style e.g. 1.0\cdot 10^{-5}      
  };
  enum VectorStyle { HorizontalVector, CAPDStyle, VerticalVector};
  enum EquationStyle {GlobalEquationStyle=-1, NoEquation, InlineEquation, DisplayEquation };
  enum PlusSymbolStyle { StandardPlus, ImplicitPlus, NoPlus };
  FloatStyle floatStyle;
  std::string baseSymbol;
  std::string iSymbol;
  std::string vectorStartSymbol;
  std::string vectorEndSymbol;
  std::string vectorSeparator;
  std::string equationBeginSymbol;
  std::string equationEndSymbol;
  std::string plusSymbol;


  TexWriter(std::ostream & output)
  : out(output),
    floatStyle(FloatFix) {
    setBaseStyle(CStyle);
    setISymbol(ImExplicit);
    setVectorStyle(HorizontalVector);
    setEquationSymbol(NoEquation);
    setPlusSymbol(StandardPlus);
  }
  // Returns number of digits before decimal point in given number n
  static inline int numberOfDigits(double n){  return abs(n)>1 ? floor( log10( abs( n ) ) ) + 1 : 1; }
  
  template <typename T>
  T putDigits(std::deque<int> & container, int & sign, T number, int n, int precision);
  void writeSciInterval(const std::deque<int> & dl, const std::deque<int> & dr,
                        int p, int l_sign, int len );
  void writeFixInterval(const std::deque<int> & dl, const std::deque<int> & dr,
                        int p, int l_sign, int intDigits, int prec);
  void writeTxtInterval(const std::deque<int> & dl, const std::deque<int> & dr,
                        int p, int l_sign, int intDigits, int prec);
  template<typename Float>
  static void computeNumberOfDigitsForInterval(Float & l, Float & r, int prec, FloatStyle style, int & n, int & len, int & exponent);
 
  template<typename IntervalType>
  void writeInterval(const IntervalType& intv);
  template <typename VectorType>
  void writeVector(const VectorType & x);

  int precision(){
    return out.precision();
  }
  TexWriter &  precision(int newPrecision){
    out.precision(newPrecision);
    return *this;
  }
  TexWriter & setFloatStyle(FloatStyle style){
    floatStyle = style;
    return *this;
  }
  /// Sets style of a base.
  TexWriter &  setExponentStyle(ExponentStyle i){
    baseSymbol = (i==CStyle)? "\\mbox{e}" :
                 (i==TexStyle)? "\\cdot 10^" :
                   "";
    return *this;
  }
  /// Sets style of a base.
  /// i=0 C convention for exponents e.g 1.0e-5
  /// i=1 TeX e.g. 1.0\cdot 10^{-5}
  TexWriter & setBaseStyle(int i){
    return setExponentStyle(ExponentStyle(i));
  }
  
  /// Sets style of imaginary symbol i.
  TexWriter &  setISymbol(ImStyle i){
      iSymbol = (i!=ImExplicit)? "\\ii " : "\\mbox{i\\,}";
      return *this;
  }
  /// sets imaginary symbol i.
  /// @deprecated
  TexWriter &  setISymbol(int i){ return setISymbol(ImStyle(i)); }
  
  /// sets vector style
  TexWriter &  setVectorStyle(const std::string & start, const std::string & end, const std::string & separator){
    vectorStartSymbol = start;
    vectorEndSymbol = end;
    vectorSeparator = separator;
    return *this;
  }
  /**
   *  sets vector style.
   */
  TexWriter &  setVectorStyle(const VectorStyle & i){
    switch(i){
      case 0:
        setVectorStyle("(", ")", ",\\;");
        break;
      case 1:
        setVectorStyle("\\{", "\\}", ",\\;");
        break;
      case 2: // vertical vector
        setVectorStyle("\\left(\\begin{array}{c}\n", "\n\\end{array}\\right)", "\\\\\n");
        break;
    }
    return *this;
  }
  /**
   *  sets vector style.
   *  * i=0  sequence : (1, 3, 4)
   *  * i=1  CAPD vectors : {1, 3, 4}
   *  * i=2  vertical vector
   * @deprecated
   */
  TexWriter &  setVectorStyle(int i){
    return setVectorStyle(VectorStyle(i));
  }
  std::string getEquationBeginSymbol(EquationStyle i){
    return (i==NoEquation)? "":
           (i==InlineEquation)? "$":
           (i==DisplayEquation)? "\\[" : "";
  }
  std::string getEquationEndSymbol(EquationStyle i){
    return (i==NoEquation)? "":
           (i==InlineEquation)? "$":
           (i==DisplayEquation)? "\\]" : "";
  }
  TexWriter &  setEquationSymbol(EquationStyle i){
    equationBeginSymbol = getEquationBeginSymbol(i);
    equationEndSymbol = getEquationEndSymbol(i);
    return *this;
  }
  TexWriter &  setEquationSymbol(int i){
    return setEquationSymbol(EquationStyle(i));
  }
  TexWriter &  setPlusSymbol(int i){
    plusSymbol = (i==0) ? "+" :
                 (i==1) ? "\\;\\;\\;" :
                          "";
    return *this;
  }

  template <typename T>
  TexWriter & write(const T & x, EquationStyle equationStyle = GlobalEquationStyle);

  void writeDocumentHeader(std::string parameters = ""){
    out << "\\documentclass[a4paper,10pt]{article} \n"
        << "\\usepackage[utf8x]{inputenc} \n"
        << parameters
        << "\\begin{document} \n";
  }
  void writeDocumentFooter(){
    out <<  "\n\\end{document} \n";
  }

  template<typename T>
  friend TexWriter & operator<<(TexWriter & o, const T & x);
  friend TexWriter & operator<<(TexWriter & o, std::ostream& (*fn)(std::ostream&));

};

template<typename T>
TexWriter & operator<< (TexWriter & o, const T & x){
  o.out << x;
  return o;
}

template<typename T, typename R>
TexWriter & operator<< (TexWriter & o, const capd::intervals::Interval<T, R> & x){
  o.writeInterval(x);
  return o;
}

#ifdef __USE_FILIB__
template <typename T, filib::RoundingStrategy R, filib::IntervalMode M>
TexWriter & operator<< (TexWriter & o, const capd::filib::Interval<T, R, M> & x){
  o.writeInterval(x);
  return o;
}
#endif

TexWriter & operator<< (TexWriter & o, const capd::DInterval & x){
  o.writeInterval(x);
  return o;
}

template<typename T>
TexWriter & operator<< (TexWriter & o, const std::complex<T> & x){
  o << x.real();
  if (x.imag() < 0 )
    o<< " - " << o.iSymbol << -x.imag();
  else
    o <<" + " << o.iSymbol << x.imag();
  return o;
}

template<typename T, capd::vectalg::__size_type dim>
TexWriter & operator<< (
    TexWriter & o,
    const capd::vectalg::Vector<T, dim> & x){
  o.writeVector(x);
  return o;
}

inline TexWriter & operator<<(TexWriter & o, std::ostream& (*fn)(std::ostream&)) {
    fn(o.out);
    return o;
}



//inline std::ostream & operator << (std::ostream & o , const std::deque<int> & c){
//  for(unsigned int i=0; i < c.size(); i++)
//    o << c[i];
//  return o;
//}

///////////////////////////////////////////////////////////////////////////////
//  putDigits
///
/// It parses \b number and puts its \b nDigits digits in \b container.
///
/// @return part of number that remained after parsing
///////////////////////////////////////////////////////////////////////////////
template <typename T>
T TexWriter::putDigits(
    std::deque<int> & container,   ///< [out] container where digits will be stored
    int & sign,                    ///< [out] returns sign of number
    T number,                      ///< [in]  number to be parsed
    int nDecDigits,                ///< [in]  number of decimal digits (before decimal point) of a given number to be stored in container
    int nDigits                    ///< [in]  total number of digits to be stored in container
){

//  std::cout << "putDigits  number = " << number << "\n nDecDigits : " << nDecDigits << "\n nDigits : " << nDigits << "\n"; 
  if(nDecDigits<0 || nDigits<0)
    throw std::range_error("putDigits: number of digits is negative.");
  
  sign = (number==0)? 0 : (number < 0)? -1 : 1;

  if(sign == -1)
    number = -number;

  long long intPart = static_cast<long long>(toDouble(number));
//  std::cout <<"\n number : "<< number <<  "\n  intPart " << intPart;

  number -= (double)intPart;
//    std::cout <<"\n number : "<< number ;

  int i=0;
  for(; i<nDecDigits; ++i){
    container.push_front(intPart % 10);
    intPart /= 10;
//    std::cout <<"\n i : "<< i <<  "  intPart " << intPart << " digit " << container[0];
  }

  for(; i<nDigits; ++i){
    number *= 10;
    int digit = toInt(number);
    container.push_back(digit);
    number -= digit;
//    std::cout << "\n i : "<< i <<  "  number " << number << " digit " << container[i];
  }
  container.resize(nDigits);
  return (nDecDigits<=nDigits)? number : T(1.0);
}

/// Rounds up parsed number stored in d
void roundUp(
    std::deque<int> & d    ///< container with consecutive digits of a number
 ){
  int carry = 1;
  for(int i=d.size()-1; i>=0 && carry; --i){
    ++d[i];
    if(d[i]>=10){
      d[i]=0;
    } else {
      carry = 0;
    }
  }
}

void TexWriter::writeSciInterval(
    const std::deque<int> & dl,
    const std::deque<int> & dr,
    int p,
    int l_sign,
    int len
){

int z=0;
 if(l_sign < 0){   // r_sign is always positive
   out << "{}_{-";
 } else {
   while(z<len && dl[z]==dr[z]){
     out << dl[z];
     if(z==0)
       out << ".";
     ++z;
   }
   if(!z)
     out << "{}";
   out << "_{";
 }

 for(int s=z; s<len; ++s){
   out << dl[s];
   if(s==0)
     out << ".";
 }
 out << "}";

 out << "^{";
 if(!z && l_sign<0)
   out << plusSymbol;
   //out << "\\;\\;\\;";
 for(int s=z; s<len; ++s){
   out << dr[s];
   if(s==0)
     out << ".";
 }
 out << "}";

 if(p!=0){
     out << this->baseSymbol << "{" << p << "}";
 }
}

void TexWriter::writeFixInterval(
    const std::deque<int> & dl,
    const std::deque<int> & dr,
    int /*p*/,
    int l_sign,
    int intDigits,
    int len
){

  int z=0;
  bool l_started = false,  r_started=false;


  if(l_sign < 0){   // r_sign is always positive
    out << "{}_{";
    if(intDigits==0){
      out <<"-0.";
      l_started = true;
    }
  } else {
    if(intDigits==0){
      out <<"0.";
      l_started = true; r_started = true;
    }
    while(z<len && dl[z]==dr[z]){
      out << dl[z];
      l_started = true; r_started = true;
      ++z;
      if(z==intDigits)
        out << ".";
    }
    if(!z)
      out << "{}";
    out << "_{";
  }
  for(int s=z; s<len; ++s){
    if(l_started) {
           out << dl[s];
      }else{
        if( dl[s]!=0){
           l_started = true;
           if(l_sign < 0)
             out << "-";
           out << dl[s];
         } else {
           if(s!=intDigits-1)
             out << "\\;\\;";
         }
      }

    if(s==intDigits-1){
      out << ((l_started)? "." : (l_sign < 0)? "-0.":"0.");
      l_started = true;
    }
  }
 out << "}";

 out << "^{";
 if(l_sign<0){
   out << plusSymbol;
   if(intDigits==0){
     out <<"0.";
     r_started = true;
   }
 }
   //out << "\\;\\;\\;";

 for(int s=z; s<len; ++s){

   if(r_started) {
        out << dr[s];
   }else{
     if( dr[s]!=0){
        r_started = true;
        out << dr[s];
      } else {
        if(s!=intDigits-1)
          out << "\\;\\;";
      }
   }
      if(s==intDigits-1){
        out << ((r_started)? "." : "0.");
        r_started = true;
      }
 }
 out << "}";

// if(p!=0){
//     out << this->baseSymbol << "{" << p << "}";
// }
}

void TexWriter::writeTxtInterval(
    const std::deque<int> & dl,
    const std::deque<int> & dr,
    int /*p*/,
    int l_sign,
    int intDigits,
    int len
){

  int z=0;
  bool l_started = false,  r_started=false;


  if(l_sign < 0){   // r_sign is always positive
    out << "[";
    if(intDigits==0){
      out <<"-0.";
      l_started = true;
    }
  } else {
    if(intDigits==0){
      out <<"0.";
      l_started = true; r_started = true;
    }
    while(z<len && dl[z]==dr[z]){
      out << dl[z];
      l_started = true; r_started = true;
      ++z;
      if(z==intDigits)
        out << ".";
    }
    out << "[";
  }
  for(int s=z; s<len; ++s){
    if(l_started) {
           out << dl[s];
      }else{
        if( dl[s]!=0){
           l_started = true;
           if(l_sign < 0)
             out << "-";
           out << dl[s];
         }
      }

    if(s==intDigits-1){
      out << ((l_started)? "." : (l_sign < 0)? "-0.":"0.");
      l_started = true;
    }
  }

 out << ",";
 if(l_sign<0){
   out << plusSymbol;
   if(intDigits==0){
     out <<"0.";
     r_started = true;
   }
 }

 for(int s=z; s<len; ++s){

   if(r_started) {
        out << dr[s];
   }else{
     if( dr[s]!=0){
        r_started = true;
        out << dr[s];
      }
   }
      if(s==intDigits-1){
        out << ((r_started)? "." : "0.");
        r_started = true;
      }
 }
 out << "]";
}
/**
 * For given left and right end of interval and expected precision depending on style chosen
 * it computes needed total number of digits, number of digits before decimal point. 
 * Additionally for FloatSci style it returns common exponent and changes endpoints so 
 * that at least one of them is in the form 1.xxxxxx... (second can be 0.00xxxx)  
 */
template<typename Float>
void TexWriter::computeNumberOfDigitsForInterval(
    Float & leftEnd, 
    Float & rightEnd, 
    int precision, 
    FloatStyle style,
    int & numberOfDecDigits, 
    int & totalNumberOfDigits, 
    int & exponent
){
 
  int nl = numberOfDigits(toDouble(leftEnd));
  int nr = numberOfDigits(toDouble(rightEnd));
  numberOfDecDigits = std::max(nl, nr);
  
  switch(style){
    case FloatSci:
      totalNumberOfDigits = precision + 1;
      if(capd::abs(leftEnd)>=1 || capd::abs(rightEnd)>=1){
        exponent = numberOfDecDigits - 1;
        Float divisor = power(Float(10.), exponent);
        setRounding<Float>(capd::rounding::RoundDown);
        leftEnd/=divisor;
        setRounding<Float>(capd::rounding::RoundUp);
        rightEnd/=divisor;
        
      } else {
        // we determine negative exponent
        exponent = 0;
        while(capd::abs(leftEnd)<1 && capd::abs(rightEnd)<1) {
          --exponent;
          setRounding<Float>(capd::rounding::RoundDown);
          leftEnd*=10;
          setRounding<Float>(capd::rounding::RoundUp);
          rightEnd*=10;
          std::cout << "\n e " << exponent << " l " << leftEnd << " r " << rightEnd << "\n";
        }
      }
      numberOfDecDigits=1;
      break;
    case FloatFix:
    case FloatTxt:
      totalNumberOfDigits = numberOfDecDigits + precision;
      break;
    case FloatFix2:
      totalNumberOfDigits = std::max(numberOfDecDigits, precision);
      break;
  }
}
template<typename IntervalType>
void TexWriter::writeInterval(const IntervalType& intv){

  typedef typename IntervalType::BoundType Float;
  Float l = intv.leftBound();
  Float r = intv.rightBound();
  if(l==0 && r==0){
    out << printToString(l, out.precision());
    return;
  }
  bool isNegative = false;
  if(r<=0){
    Float t = r;  r = -l;   l = -t;
    isNegative = true;
  }
//  int nl = numberOfDigits(toInt(l));
//  int nr = numberOfDigits(toInt(r));
  
  int n=0;
  int len=0;
  int exponent=0;
  computeNumberOfDigitsForInterval(l, r, precision(), floatStyle, n, len, exponent);
  
  if(n > std::numeric_limits<long long>::digits10 ){
    out << intv;
    return;
  }


  std::deque<int> dl, dr;
  int l_sign, r_sign;

  setRounding<Float>(capd::rounding::RoundDown);
  Float l_rem = putDigits(dl, l_sign, l, n, len);

  setRounding<Float>(capd::rounding::RoundUp);
  Float r_rem = putDigits(dr, r_sign, r, n, len);

  assert(r_sign > 0);

  if(isNegative) 
    out << "-";
 
  if(l_sign < 0 && l_rem != 0)
    roundUp(dl);
  if(r_sign > 0 && r_rem != 0)
    roundUp(dr);
  switch(floatStyle){
    case FloatSci:
      writeSciInterval(dl, dr, exponent, l_sign, len );
      break;
    case FloatFix:
    case FloatFix2:
      writeFixInterval(dl, dr, exponent, l_sign, n, len );
      break;
    case FloatTxt:
      writeTxtInterval(dl, dr, exponent, l_sign, n, len );
  }
}

//template<typename IntervalType>
//void TexWriter::writeInterval(const IntervalType& intv){
//
//
//  typedef typename IntervalType::BoundType Float;
//  Float l = intv.leftBound();
//  Float r = intv.rightBound();
//  if(l==0 && r==0){
//    out << printToString(l, out.precision(), capd::rounding::RoundCut);
//    return;
//  }
//  if(r<=0){
//    Float t = r;  r = -l;   l = -t;
//    out << "-";
//  }
//  int nl = numberOfDigits(int(l));
//  int nr = numberOfDigits(int(r));
//
//  std::deque<int> dl, dr;
//
//  // we determine exponent
//  int p=0;
//  while((capd::abs(l)>10) || (capd::abs(r)>10)) {
//    ++p;
//    setRounding<Float>(capd::rounding::RoundDown);
//    l/=10;
//    setRounding<Float>(capd::rounding::RoundUp);
//    r/=10;
//  }
//
//  while(capd::abs(l)<1 && capd::abs(r)<1) {
//    --p;
//    setRounding<Float>(capd::rounding::RoundDown);
//    l*=10;
//    setRounding<Float>(capd::rounding::RoundUp);
//    r*=10;
//  }
//
//  int prec = out.precision();
//  std::string t1 = printToString(l, prec, capd::rounding::RoundDown),
//              t2 = printToString(r, prec, capd::rounding::RoundUp);
//  int len1 = t1.size(),
//      len2 = t2.size(),
//      len = std::min(len1,len2);
//  int z=0, s;
//  while(z<len && t1[z]==t2[z]){
//    out << t1[z];
//    ++z;
//  }
//  if(!z)
//    out << "{}";
//
//  out << "_{";
//  for(s=z; s<len1; ++s)
//    out << t1[s];
//  for(;s<=prec;++s)
//    out << "0";
//  out << "}";
//
//  out << "^{";
//  if(!z && l<0)
//    out << "\\;\\;\\;";
//  for(s=z;s<len2;++s)
//    out << t2[s];
//  for(;s<=prec;++s)
//    out << "0";
//  out << "}";
//
//  if(p!=0){
//      out << this->base << "{" << p << "}";
//  }
//}

template <typename VectorType>
void TexWriter::writeVector(const VectorType & v){

  out << vectorStartSymbol;
  if(v.dimension()>0){
    (*this) << v[0];
  }
  for(unsigned i=1; i<v.dimension(); i++) {
    out << vectorSeparator;
    (*this) << v[i];
  }
  out << vectorEndSymbol;
}

template <typename T>
TexWriter & TexWriter::write(const T & x, EquationStyle eqStyle ){
  out << ((eqStyle<0) ? equationBeginSymbol : getEquationBeginSymbol(eqStyle) );
  (*this)<< x;
  out << ((eqStyle<0) ? equationEndSymbol : getEquationEndSymbol(eqStyle) );
  return (*this);
}


} // end of namespace capd

#endif /* _CAPD_BASICALG_TEXWRITER_H_ */
