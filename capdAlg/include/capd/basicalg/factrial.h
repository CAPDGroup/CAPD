
/////////////////////////////////////////////////////////////////////////////
/// @file factrial.h
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.
#include <vector>
#ifndef _CAPD_CAPD_FACTRIAL_H_ 
#define _CAPD_CAPD_FACTRIAL_H_ 
#include <stdexcept>
namespace capd{

class Newton{
public:
  static Newton& getInstance(){
    return instance;
  }
  unsigned long long factorial(unsigned n);
//  unsigned long long newton(unsigned n, unsigned k);
  unsigned long long newton(unsigned n, unsigned k)
{
  unsigned first_undefined_index=index(first_unknown_newton_level,0);
	if(n>=first_unknown_newton_level){
      newton_storage.resize(index(n+1,0));
		if(first_undefined_index == 0){
			newton_storage[first_undefined_index++]=1;
			first_unknown_newton_level++;
		}
		for(unsigned m=first_unknown_newton_level;m<=n;m++){
			newton_storage[first_undefined_index]=newton_storage[first_undefined_index+m]=1;
			for(unsigned p=1;p<m;p++) newton_storage[first_undefined_index+p]=
					newton_storage[index(m-1,p-1)]+newton_storage[index(m-1,p)];
			first_undefined_index+=(m+1);
		}
		first_unknown_newton_level=n+1;
      }
	return newton_storage[index(n,k)];
}

private:
  Newton() : first_unknown_factorial(0), first_unknown_newton_level(0) {}
  std::vector<unsigned long long> factorial_storage;
  unsigned first_unknown_factorial;
  std::vector<unsigned long long> newton_storage;
  unsigned first_unknown_newton_level;

  inline unsigned index(unsigned n,unsigned k)
  {
    return n*(n+1)/2+k;
  }
  static Newton instance;
};

/// @addtogroup basicalg
/// @{
template <long N, long K>
struct Binomial
{
   static const long value = Binomial<N-1,K-1>::value + Binomial<N-1,K>::value;
};

template <long K>
struct Binomial<0,K>
{
   static const long value = 0;
};

template <long N>
struct Binomial<N,0>
{
   static const long value = 1;
};

template<long N>
struct Binomial<N,N>
{
   static const long value=1;
};

template <>
struct Binomial<0,0>
{
   static const long value = 1;
};

template<long N>
struct Factorial
{
   static const long value = N*Factorial<N-1>::value;
};

template<>
struct Factorial<1>
{
   static const long value = 1;
};

template<>
struct Factorial<0>
{
   static const long value = 1;
};

inline double realFactorial(unsigned n){
  static double factorials[] = {
    1.0,
    1.0,
    2.0,
    6.0,
    24.0,
    120.0,
    720.0,
    5040.0,
    40320.0,
    362880.0,
    3628800.0,
    39916800.0,
    479001600.0,
    6227020800.0,
    87178291200.0,
    1307674368000.0,
    20922789888000.0,
    355687428096000.0,
    6402373705728000.0,
    121645100408832000.0,
    2432902008176640000.0
  };
  if(n<=20) 
    return factorials[n];
  else
    throw std::runtime_error("Function realFactorial: n! exceeded double capacity.");
}
/// @}

} // namespace capd

///< compute and store n factorial
inline unsigned long long factorial(unsigned n){
  return capd::Newton::getInstance().factorial(n);
}


///< compute and store newton symbol (n \over k)
inline unsigned long long binomial(unsigned n, unsigned k){
  return capd::Newton::getInstance().newton(n,k);
}

#endif // _CAPD_CAPD_FACTRIAL_H_ 
