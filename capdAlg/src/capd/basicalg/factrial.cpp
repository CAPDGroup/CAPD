/// @addtogroup capd
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file factrial.cpp
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include <cstddef>
#include "capd/basicalg/factrial.h"

using namespace capd;

unsigned long long Newton::factorial(unsigned n)
{
	if(n<first_unknown_factorial){
		return factorial_storage[n];
	}else{
	  unsigned i;
    factorial_storage.resize(n+1);
		if(first_unknown_factorial == 0){
			factorial_storage[first_unknown_factorial++]=1;
		}
		unsigned long long result=factorial_storage[first_unknown_factorial-1];
		for(i=first_unknown_factorial;i<=n;i++) factorial_storage[i]=result*=i;
		first_unknown_factorial=n+1;
		return result;
	}
}
/*
unsigned long long Newton::newton(unsigned n, unsigned k)
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
*/
Newton Newton::instance;
/// @}
