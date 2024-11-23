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

/* CVS $Id: q_errm.cpp,v 1.15 2014/01/30 17:23:55 cxsc Exp $ */

#ifndef Q_ERRM_CPP
#define Q_ERRM_CPP

#include "fi_lib.hpp" 

namespace fi_lib {

 using namespace std;
 using cxsc::real;

/* ----------- error handling for argument == NaN ---------------*/
 real q_abortnan(int n, real *x, int fctn){
	cerr << endl << "*** Error in fi_lib (V1.3): Function: ";
	switch(fctn){
	    case  0: cerr << "q_sqrt" ; break; 
	    case  1: cerr << "q_sqr " ; break; 
	    case  2: cerr << "q_exp " ; break; 
	    case  3: cerr << "q_epm1" ; break; 
	    case  4: cerr << "q_exp2" ; break; 
    	    case  5: cerr << "q_ex10" ; break; 
	    case  6: cerr << "q_log " ; break; 
	    case  7: cerr << "q_lg1p" ; break; 
	    case  8: cerr << "q_log2" ; break; 
	    case  9: cerr << "q_lg10" ; break; 
	    case 10: cerr << "q_sin " ; break; 
	    case 11: cerr << "q_cos " ; break; 
	    case 12: cerr << "q_tan " ; break; 
	    case 13: cerr << "q_cot " ; break; 
	    case 14: cerr << "q_asin" ; break; 
	    case 15: cerr << "q_acos" ; break; 
	    case 16: cerr << "q_atan" ; break; 
	    case 17: cerr << "q_acot" ; break; 
	    case 18: cerr << "q_sinh" ; break; 
	    case 19: cerr << "q_cosh" ; break; 
	    case 20: cerr << "q_tanh" ; break; 
	    case 21: cerr << "q_coth" ; break; 
	    case 22: cerr << "q_asnh" ; break; 
	    case 23: cerr << "q_acnh" ; break; 
	    case 24: cerr << "q_atnh" ; break; 
	    case 25: cerr << "q_acth" ; break; 
	    case 27: cerr << "q_erf " ; break; 
	    case 28: cerr << "q_erfc" ; break; 
	}
	cerr << endl << "*** Error in fi_lib (V1.3): Argument == NaN ! ***" << endl;
	exit(n);
	return(*x); 
}

/* ------------ error handling for point arguments ---------------*/
 real q_abortr1(int n, real *x, int fctn){
	cerr << endl << "*** Error in fi_lib (V1.3): Function: ";
	switch(fctn){
	    case  0: cerr << "q_sqrt" ; break; 
	    case  1: cerr << "q_sqr " ; break; 
	    case  2: cerr << "q_exp " ; break; 
	    case  3: cerr << "q_epm1" ; break; 
	    case  4: cerr << "q_exp2" ; break; 
	    case  5: cerr << "q_ex10" ; break; 
	    case  6: cerr << "q_log " ; break; 
	    case  7: cerr << "q_lg1p" ; break; 
	    case  8: cerr << "q_log2" ; break; 
	    case  9: cerr << "q_lg10" ; break; 
	    case 10: cerr << "q_sin " ; break; 
	    case 11: cerr << "q_cos " ; break; 
	    case 12: cerr << "q_tan " ; break; 
	    case 13: cerr << "q_cot " ; break; 
	    case 14: cerr << "q_asin" ; break; 
	    case 15: cerr << "q_acos" ; break; 
	    case 16: cerr << "q_atan" ; break; 
	    case 17: cerr << "q_acot" ; break; 
	    case 18: cerr << "q_sinh" ; break; 
	    case 19: cerr << "q_cosh" ; break; 
	    case 20: cerr << "q_tanh" ; break; 
	    case 21: cerr << "q_coth" ; break; 
	    case 22: cerr << "q_asnh" ; break; 
	    case 23: cerr << "q_acnh" ; break; 
	    case 24: cerr << "q_atnh" ; break; 
	    case 25: cerr << "q_acth" ; break; 
	    case 26: cerr << "q_comp" ; break;
	    case 27: cerr << "q_erf " ; break; 
	    case 28: cerr << "q_erfc" ; break; 
	}
	if (n==INV_ARG) 
	    cerr << endl << "*** Error in fi_lib (V1.3): Invalid argument ! ***" << endl;
	else
	    cerr << endl << "*** Error in fi_lib (V1.3): Overflow (result) ! ***" << endl;
	//cerr << "*** Error in fi_lib (V1.3): Argument x = %24.15e \n",*x);
	cerr << "*** Error in fi_lib (V1.3): Argument x =  " << *x << endl;
	exit(n);
	return(*x);
}


/* ------------- error handling for interval arguments -------------*/
 interval q_abortr2(int n, real *x1, real *x2, int fctn){
	interval res;
	cerr << "*** Error in fi_lib (V1.3): Function: ";
	switch(fctn){
	    case  0: cerr << "j_sqrt" ; break; 
	    case  1: cerr << "j_sqr " ; break; 
	    case  2: cerr << "j_exp " ; break; 
	    case  3: cerr << "j_epm1" ; break; 
	    case  4: cerr << "j_exp2" ; break; 
	    case  5: cerr << "j_ex10" ; break; 
	    case  6: cerr << "j_log " ; break; 
	    case  7: cerr << "j_lg1p" ; break; 
	    case  8: cerr << "j_log2" ; break; 
	    case  9: cerr << "j_lg10" ; break; 
	    case 10: cerr << "j_sin " ; break; 
	    case 11: cerr << "j_cos " ; break; 
	    case 12: cerr << "j_tan " ; break; 
	    case 13: cerr << "j_cot " ; break; 
	    case 14: cerr << "j_asin" ; break; 
	    case 15: cerr << "j_acos" ; break; 
	    case 16: cerr << "j_atan" ; break; 
	    case 17: cerr << "j_acot" ; break; 
	    case 18: cerr << "j_sinh" ; break; 
	    case 19: cerr << "j_cosh" ; break; 
	    case 20: cerr << "j_tanh" ; break; 
	    case 21: cerr << "j_coth" ; break; 
	    case 22: cerr << "j_asnh" ; break; 
	    case 23: cerr << "j_acnh" ; break; 
	    case 24: cerr << "j_atnh" ; break; 
	    case 25: cerr << "j_acth" ; break; 
	    case 27: cerr << "q_erf " ; break; 
	    case 28: cerr << "q_erfc" ; break; 
	}
	if (n==INV_ARG) 
	    cerr << endl << "*** Error in fi_lib (V1.3): Invalid argument ! ***" << endl;
	else 
	    cerr << endl << "*** Error in fi_lib (V1.3): Overflow (result) ! ***" << endl;
	cerr << "*** Error in fi_lib (V1.3): Argument x.INF = %24.15e " << *x1 << endl;
	cerr << "*** Error in fi_lib (V1.3): Argument x.SUP = %24.15e " << *x2 << endl;
	exit(n);

	Inf(res)=*x1; Sup(res)=*x2;
	return(res);
    }

/* ------------ error handling for point arguments ---------------*/
 real q_abortdivd(int n, real *x){
	cerr << endl << "*** Error in fi_lib (V1.3): Function: div_id";
	cerr << endl << "*** Error in fi_lib (V1.3): Division by zero ! ***" << endl;
	cerr << "*** Error in fi_lib (V1.3): x = %24.15e \n" << *x << endl;;
	exit(n);
	return(*x);
 }

/* ------------- error handling for interval arguments -------------*/
 interval q_abortdivi(int n, real *x1, real *x2){
	interval res;
	cerr << endl << "*** Error in fi_lib (V1.3): Function: div_ii";
	cerr << endl << "*** Error in fi_lib (V1.3): Division by zero ! ***" << endl;
	cerr << "*** Error in fi_lib (V1.3): x.INF = %24.15e" << *x1 << endl;
	cerr << "*** Error in fi_lib (V1.3): x.SUP = %24.15e" << *x2 << endl;
	exit(n);
	Inf(res)=*x1; Sup(res)=*x2;
	return(res);
 } 

} // Namespace

#endif





