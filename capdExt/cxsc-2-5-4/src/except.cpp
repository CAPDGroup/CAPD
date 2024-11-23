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

/* CVS $Id: except.cpp,v 1.25 2014/01/30 17:23:45 cxsc Exp $ */

#include <string>
#include <iostream>
#include <xscclass.hpp>
#include "except.hpp"

/*
#if ROUND_C99_SAVE+ROUND_C99_QUICK
#include <fenv.h>
#endif

#if ROUND_C96_SAVE+ROUND_C96_QUICK
#include <fenv96.h>
#endif

#if SUN4_FORTE&&(ROUND_C96_SAVE+ROUND_C96_QUICK+ROUND_C99_SAVE+ROUND_C99_QUICK)
#if SUN_STUDIO_10
extern "C"{ void presub(){} }
#else
extern "C"{ void presub(int,fex_info_t*){} }
#endif
#endif
*/

namespace cxsc {

/*
//Definitions of member functions of class cxsc_status:

cxsc_status::cxsc_status()                     //constructor is private
{ 
//   std::cout << "  *** class C-XSC Status available *** " << std::endl;

#if ROUND_C96_SAVE+ROUND_C96_QUICK+ROUND_C99_SAVE+ROUND_C99_QUICK
#if SUN4_FORTE
  fex_set_handling(FEX_COMMON,FEX_NOHANDLER,presub);
#else
  feenableexcept(FE_DIVBYZERO|FE_OVERFLOW|FE_INVALID); //(***)
#endif
#endif
}

cxsc_status* cxsc_status::getObjectReference() //static member function
{
   if (cxsc_statusRef == 0)        //function is called for the first time?
   {
      cxsc_statusRef= new cxsc_status(); //create an object of class cxsc_status 
   }
   return cxsc_statusRef;          //return pointer to the object created
}


//definition of static attribut cxsc_statusRef declared in class cxsc_status:
cxsc_status* cxsc_status::cxsc_statusRef=0; 
*/


//----------------------  Error Modules --------------------------------------

//member functions of class ERROR_ALL:
int ERROR_ALL::errnum() const throw() { return 1000; }
string ERROR_ALL::errtext() const throw() { return fkt+": ERROR_ALL"; }
ERROR_ALL::ERROR_ALL() throw() : fkt("<unknown function>") { }
ERROR_ALL::ERROR_ALL(const string &f) throw():fkt(f) { }
ERROR_ALL::~ERROR_ALL() throw() { }

//member functions of class ERROR_DOT 
int ERROR_DOT::errnum() const throw() { return 2000; }
string ERROR_DOT::errtext() const throw() { return fkt+": ERROR_DOT "; }
ERROR_DOT::ERROR_DOT() throw() { fkt="<unknown function>"; }
ERROR_DOT::ERROR_DOT(const string &f) throw() { fkt=f; }

//member functions of class ERROR_REAL 
int ERROR_REAL::errnum() const throw() { return 3000; }
string ERROR_REAL::errtext() const throw() { return fkt+": ERROR_REAL "; }
ERROR_REAL::ERROR_REAL() throw() { fkt="<unknown function>"; }
ERROR_REAL::ERROR_REAL(const string &f) throw() { fkt=f; }

//member functions of class ERROR_INTERVAL 
int ERROR_INTERVAL::errnum() const throw() { return 4000; }
string ERROR_INTERVAL::errtext() const throw() { return fkt+": ERROR_INTERVAL "; }
ERROR_INTERVAL::ERROR_INTERVAL() throw() { fkt="<unknown function>"; }
ERROR_INTERVAL::ERROR_INTERVAL(const string &f) throw() { fkt=f; }

//member functions of class ERROR_COMPLEX 
int ERROR_COMPLEX::errnum() const throw() { return 5000; }
string ERROR_COMPLEX::errtext() const throw() { return fkt+": ERROR_COMPLEX "; }
ERROR_COMPLEX::ERROR_COMPLEX() throw() { fkt="<unknown function>"; }
ERROR_COMPLEX::ERROR_COMPLEX(const string &f) throw() { fkt=f; }

//member functions of class ERROR_CINTERVAL 
int ERROR_CINTERVAL::errnum() const throw() { return 6000; }
string ERROR_CINTERVAL::errtext() const throw() { return fkt+": ERROR_CINTERVAL "; }
ERROR_CINTERVAL::ERROR_CINTERVAL() throw() { fkt="<unknown function>"; }
ERROR_CINTERVAL::ERROR_CINTERVAL(const string &f) throw() { fkt=f; }

//member functions of class ERROR_VECTOR 
int ERROR_VECTOR::errnum() const throw() { return 51000; }
string ERROR_VECTOR::errtext() const throw() { return fkt+": ERROR_VECTOR "; }
ERROR_VECTOR::ERROR_VECTOR() throw() { fkt="<unknown function>"; }
ERROR_VECTOR::ERROR_VECTOR(const string &f) throw() { fkt=f; }

//member functions of class ERROR_RVECTOR 
int ERROR_RVECTOR::errnum() const throw() { return 7000; }
string ERROR_RVECTOR::errtext() const throw() { return fkt+": ERROR_RVECTOR "; }
ERROR_RVECTOR::ERROR_RVECTOR() throw() { fkt="<unknown function>"; }
ERROR_RVECTOR::ERROR_RVECTOR(const string &f) throw() { fkt=f; }

//member functions of class ERROR_IVECTOR 
int ERROR_IVECTOR::errnum() const throw() { return 8000; }
string ERROR_IVECTOR::errtext() const throw() { return fkt+": ERROR_IVECTOR "; }
ERROR_IVECTOR::ERROR_IVECTOR() throw() { fkt="<unknown function>"; }
ERROR_IVECTOR::ERROR_IVECTOR(const string &f) throw() { fkt=f; }

//member functions of class ERROR_CVECTOR 
int ERROR_CVECTOR::errnum() const throw() { return 9000; }
string ERROR_CVECTOR::errtext() const throw() { return fkt+": ERROR_CVECTOR "; }
ERROR_CVECTOR::ERROR_CVECTOR() throw() { fkt="<unknown function>"; }
ERROR_CVECTOR::ERROR_CVECTOR(const string &f) throw() { fkt=f; }

//member functions of class ERROR_CIVECTOR 
int ERROR_CIVECTOR::errnum() const throw() { return 10000; }
string ERROR_CIVECTOR::errtext() const throw() { return fkt+": ERROR_CIVECTOR "; }
ERROR_CIVECTOR::ERROR_CIVECTOR() throw() { fkt="<unknown function>"; }
ERROR_CIVECTOR::ERROR_CIVECTOR(const string &f) throw() { fkt=f; }

//member functions of class ERROR_MATRIX 
int ERROR_MATRIX::errnum() const throw() { return 52000; }
string ERROR_MATRIX::errtext() const throw() { return fkt+": ERROR_MATRIX "; }
ERROR_MATRIX::ERROR_MATRIX() throw() { fkt="<unknown function>"; }
ERROR_MATRIX::ERROR_MATRIX(const string &f) throw() { fkt=f; }

//member functions of class ERROR_RMATRIX 
int ERROR_RMATRIX::errnum() const throw() { return 11000; }
string ERROR_RMATRIX::errtext() const throw() { return fkt+": ERROR_RMATRIX "; }
ERROR_RMATRIX::ERROR_RMATRIX() throw() { fkt="<unknown function>"; }
ERROR_RMATRIX::ERROR_RMATRIX(const string &f) throw() { fkt=f; }

//member functions of class ERROR_IMATRIX 
int ERROR_IMATRIX::errnum() const throw() { return 12000; }
string ERROR_IMATRIX::errtext() const throw() { return fkt+": ERROR_IMATRIX "; }
ERROR_IMATRIX::ERROR_IMATRIX() throw() { fkt="<unknown function>"; }
ERROR_IMATRIX::ERROR_IMATRIX(const string &f) throw() { fkt=f; }

//member functions of class ERROR_CMATRIX 
int ERROR_CMATRIX::errnum() const throw() { return 13000; }
string ERROR_CMATRIX::errtext() const throw() { return fkt+": ERROR_CMATRIX "; }
ERROR_CMATRIX::ERROR_CMATRIX() throw() { fkt="<unknown function>"; }
ERROR_CMATRIX::ERROR_CMATRIX(const string &f) throw() { fkt=f; }

//member functions of class ERROR_CIMATRIX 
int ERROR_CIMATRIX::errnum() const throw() { return 14000; }
string ERROR_CIMATRIX::errtext() const throw() { return fkt+": ERROR_CIMATRIX "; }
ERROR_CIMATRIX::ERROR_CIMATRIX() throw() { fkt="<unknown function>"; }
ERROR_CIMATRIX::ERROR_CIMATRIX(const string &f) throw() { fkt=f; }

//member functions of class ERROR_LREAL
int ERROR_LREAL::errnum() const throw() { return 15000; }
string ERROR_LREAL::errtext() const throw() { return fkt+": ERROR_LREAL "; }
ERROR_LREAL::ERROR_LREAL() throw() { fkt="<unknown function>"; }
ERROR_LREAL::ERROR_LREAL(const string &f) throw() { fkt=f; }

//member functions of class ERROR_LINTERVAL 
int ERROR_LINTERVAL::errnum() const throw() { return 16000; }
string ERROR_LINTERVAL::errtext() const throw() { return fkt+": ERROR_LINTERVAL "; }
ERROR_LINTERVAL::ERROR_LINTERVAL() throw() { fkt="<unknown function>"; }
ERROR_LINTERVAL::ERROR_LINTERVAL(const string &f) throw() { fkt=f; }

//member functions of class ERROR_LRVECTOR 
int ERROR_LRVECTOR::errnum() const throw() { return 17000; }
string ERROR_LRVECTOR::errtext() const throw() { return fkt+": ERROR_LRVECTOR "; }
ERROR_LRVECTOR::ERROR_LRVECTOR() throw() { fkt="<unknown function>"; }
ERROR_LRVECTOR::ERROR_LRVECTOR(const string &f) throw() { fkt=f; }

//member functions of class ERROR_LIVECTOR 
int ERROR_LIVECTOR::errnum() const throw() { return 18000; }
string ERROR_LIVECTOR::errtext() const throw() { return fkt+": ERROR_LIVECTOR "; }
ERROR_LIVECTOR::ERROR_LIVECTOR() throw() { fkt="<unknown function>"; }
ERROR_LIVECTOR::ERROR_LIVECTOR(const string &f) throw() { fkt=f; }

//member functions of class ERROR_LRMATRIX 
int ERROR_LRMATRIX::errnum() const throw() { return 19000; }
string ERROR_LRMATRIX::errtext() const throw() { return fkt+": ERROR_LRMATRIX "; }
ERROR_LRMATRIX::ERROR_LRMATRIX() throw() { fkt="<unknown function>"; }
ERROR_LRMATRIX::ERROR_LRMATRIX(const string &f) throw() { fkt=f; }

//member functions of class ERROR_LIMATRIX
int ERROR_LIMATRIX::errnum() const throw() { return 20000; }
string ERROR_LIMATRIX::errtext() const throw() { return fkt+": ERROR_LIMATRIX "; }
ERROR_LIMATRIX::ERROR_LIMATRIX() throw() { fkt="<unknown function>"; }
ERROR_LIMATRIX::ERROR_LIMATRIX(const string &f) throw() { fkt=f; }

//member functions of class ERROR_INTVECTOR
int ERROR_INTVECTOR::errnum() const throw() { return 21000; }
string ERROR_INTVECTOR::errtext() const throw() { return fkt+": ERROR_INTVECTOR "; }
ERROR_INTVECTOR::ERROR_INTVECTOR() throw() { fkt="<unknown function>"; }
ERROR_INTVECTOR::ERROR_INTVECTOR(const string &f) throw() { fkt=f; }

//member functions of class ERROR_INTMATRIX 
int ERROR_INTMATRIX::errnum() const throw() { return 22000; }
string ERROR_INTMATRIX::errtext() const throw() { return fkt+": ERROR_INTMATRIX "; }
ERROR_INTMATRIX::ERROR_INTMATRIX() throw() { fkt="<unknown function>"; }
ERROR_INTMATRIX::ERROR_INTMATRIX(const string &f) throw() { fkt=f; }

//member functions of class ERROR_TAYLOR; Blomquist 22.12.2008;
int ERROR_TAYLOR::errnum() const throw() { return 23000; }
string ERROR_TAYLOR::errtext() const throw() { return fkt+": ERROR_TAYLOR "; }
ERROR_TAYLOR::ERROR_TAYLOR() throw() { fkt="<unknown function>"; }
ERROR_TAYLOR::ERROR_TAYLOR(const string &f) throw() { fkt=f; }

//-------------------------- Error Types --------------------------------

//member functions of class WRONG_ROUNDING
int WRONG_ROUNDING::errnum() const throw() { return 1; }
string WRONG_ROUNDING::errtext() const throw() { return fkt+": WRONG_ROUNDING"; }
WRONG_ROUNDING::WRONG_ROUNDING() throw() { fkt="<unknown function>"; }
WRONG_ROUNDING::WRONG_ROUNDING(const string &f) throw() { fkt=f; }

//member functions of class NO_MORE_MEMORY
int NO_MORE_MEMORY::errnum() const throw() { return 2; }
string NO_MORE_MEMORY::errtext() const throw() { return fkt+": NO_MORE_MEMORY"; }
NO_MORE_MEMORY::NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
NO_MORE_MEMORY::NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of class WRONG_DOT_TYPE
int WRONG_DOT_TYPE::errnum() const throw() { return 3; }
string WRONG_DOT_TYPE::errtext() const throw() { return fkt+": WRONG_DOT_TYPE"; }
WRONG_DOT_TYPE::WRONG_DOT_TYPE() throw() { fkt="<unknown function>"; }
WRONG_DOT_TYPE::WRONG_DOT_TYPE(const string &f) throw() { fkt=f; }

//member functions of class NOT_AVAILABLE
int NOT_AVAILABLE::errnum() const throw() { return 4; }
string NOT_AVAILABLE::errtext() const throw() { return fkt+": NOT_AVAILABLE"; }
NOT_AVAILABLE::NOT_AVAILABLE() throw() { fkt="<unknown function>"; }
NOT_AVAILABLE::NOT_AVAILABLE(const string &f) throw() { fkt=f; }

//member functions of class DIV_BY_ZERO
int DIV_BY_ZERO::errnum() const throw() { return 10; }
string DIV_BY_ZERO::errtext() const throw() { return fkt+": DIV_BY_ZERO"; }
DIV_BY_ZERO::DIV_BY_ZERO() throw() { fkt="<unknown function>"; }
DIV_BY_ZERO::DIV_BY_ZERO(const string &f) throw() { fkt=f; }

//member functions of class EMPTY_INTERVAL
int EMPTY_INTERVAL::errnum() const throw() { return 11; }
string EMPTY_INTERVAL::errtext() const throw() { return fkt+": EMPTY_INTERVAL"; }
EMPTY_INTERVAL::EMPTY_INTERVAL() throw() { fkt="<unknown function>"; }
EMPTY_INTERVAL::EMPTY_INTERVAL(const string &f) throw() { fkt=f; }

//member functions of class OVERFLOW_ERROR
int OVERFLOW_ERROR::errnum() const throw() { return 12; }
string OVERFLOW_ERROR::errtext() const throw() { return fkt+": OVERFLOW_ERROR"; }
OVERFLOW_ERROR::OVERFLOW_ERROR() throw() { fkt="<unknown function>"; }
OVERFLOW_ERROR::OVERFLOW_ERROR(const string &f) throw() { fkt=f; }

//member functions of class IN_EXACT_CH_OR_IS
int IN_EXACT_CH_OR_IS::errnum() const throw() { return 13; }
string IN_EXACT_CH_OR_IS::errtext() const throw() { return fkt+": IN_EXACT_CH_OR_IS"; }
IN_EXACT_CH_OR_IS::IN_EXACT_CH_OR_IS() throw() { fkt="<unknown function>"; }
IN_EXACT_CH_OR_IS::IN_EXACT_CH_OR_IS(const string &f) throw() { fkt=f; }

//member functions of class WRONG_BOUNDARIES
int WRONG_BOUNDARIES::errnum() const throw() { return 200; }
string WRONG_BOUNDARIES::errtext() const throw() { return fkt+": WRONG_BOUNDARIES"; }
WRONG_BOUNDARIES::WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
WRONG_BOUNDARIES::WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of class SUB_ARRAY_TOO_BIG
int SUB_ARRAY_TOO_BIG::errnum() const throw() { return 201; }
string SUB_ARRAY_TOO_BIG::errtext() const throw() { return fkt+": SUB_ARRAY_TOO_BIG"; }
SUB_ARRAY_TOO_BIG::SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
SUB_ARRAY_TOO_BIG::SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of class RES_OR_INP_OF_TEMP_OBJ
int RES_OR_INP_OF_TEMP_OBJ::errnum() const throw() { return 202; }
string RES_OR_INP_OF_TEMP_OBJ::errtext() const throw() { return fkt+": RES_OR_INP_OF_TEMP_OBJ"; }
RES_OR_INP_OF_TEMP_OBJ::RES_OR_INP_OF_TEMP_OBJ() throw() { fkt="<unknown function>"; }
RES_OR_INP_OF_TEMP_OBJ::RES_OR_INP_OF_TEMP_OBJ(const string &f) throw() { fkt=f; }

//member functions of class ELEMENT_NOT_IN_VEC
int ELEMENT_NOT_IN_VEC::errnum() const throw() { return 203; }
string ELEMENT_NOT_IN_VEC::errtext() const throw() { return fkt+": ELEMENT_NOT_IN_VEC"; }
ELEMENT_NOT_IN_VEC::ELEMENT_NOT_IN_VEC() throw() { fkt="<unknown function>"; }
ELEMENT_NOT_IN_VEC::ELEMENT_NOT_IN_VEC(const string &f) throw() { fkt=f; }

//member functions of class ROW_OR_COL_NOT_IN_MAT
int ROW_OR_COL_NOT_IN_MAT::errnum() const throw() { return 204; }
string ROW_OR_COL_NOT_IN_MAT::errtext() const throw() { return fkt+": ROW_OR_COL_NOT_IN_MAT"; }
ROW_OR_COL_NOT_IN_MAT::ROW_OR_COL_NOT_IN_MAT() throw() { fkt="<unknown function>"; }
ROW_OR_COL_NOT_IN_MAT::ROW_OR_COL_NOT_IN_MAT(const string &f) throw() { fkt=f; }

//member functions of class WRONG_ROW_OR_COL
int WRONG_ROW_OR_COL::errnum() const throw() { return 205; }
string WRONG_ROW_OR_COL::errtext() const throw() { return fkt+": WRONG_ROW_OR_COL "; }
WRONG_ROW_OR_COL::WRONG_ROW_OR_COL() throw() { fkt="<unknown function>"; }
WRONG_ROW_OR_COL::WRONG_ROW_OR_COL(const string &f) throw() { fkt=f; }

//member functions of class TYPE_CAST_OF_THICK_OBJ
int TYPE_CAST_OF_THICK_OBJ::errnum() const throw() { return 206; }
string TYPE_CAST_OF_THICK_OBJ::errtext() const throw() { return fkt+": TYPE_CAST_OF_THICK_OBJ"; }
TYPE_CAST_OF_THICK_OBJ::TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
TYPE_CAST_OF_THICK_OBJ::TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of class OP_WITH_WRONG_DIM
int OP_WITH_WRONG_DIM::errnum() const throw() { return 207; }
string OP_WITH_WRONG_DIM::errtext() const throw() { return fkt+": OP_WITH_WRONG_DIM"; }
OP_WITH_WRONG_DIM::OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
OP_WITH_WRONG_DIM::OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of class WRONG_STAGPREC
int WRONG_STAGPREC::errnum() const throw() { return 300; }
string WRONG_STAGPREC::errtext() const throw() { return fkt+": WRONG_STAGPREC"; }
WRONG_STAGPREC::WRONG_STAGPREC() throw() { fkt="<unknown function>"; }
WRONG_STAGPREC::WRONG_STAGPREC(const string &f) throw() { fkt=f; }

//member functions of class ELEMENT_NOT_IN_LONG
int ELEMENT_NOT_IN_LONG::errnum() const throw() { return 301; }
string ELEMENT_NOT_IN_LONG::errtext() const throw() { return fkt+": ELEMENT_NOT_IN_LONG"; }
ELEMENT_NOT_IN_LONG::ELEMENT_NOT_IN_LONG() throw() { fkt="<unknown function>"; }
ELEMENT_NOT_IN_LONG::ELEMENT_NOT_IN_LONG(const string &f) throw() { fkt=f; }

//member functions of class STD_FKT_OUT_OF_DEF
int STD_FKT_OUT_OF_DEF::errnum() const throw() { return 302; }
string STD_FKT_OUT_OF_DEF::errtext() const throw() { return fkt+": STD_FKT_OUT_OF_DEF"; }
STD_FKT_OUT_OF_DEF::STD_FKT_OUT_OF_DEF() throw() { fkt="<unknown function>"; }
STD_FKT_OUT_OF_DEF::STD_FKT_OUT_OF_DEF(const string &f) throw() { fkt=f; }

//member functions of class FAK_OVERFLOW
int FAK_OVERFLOW::errnum() const throw() { return 303; }
string FAK_OVERFLOW::errtext() const throw() { return fkt+": FAK_OVERFLOW"; }
FAK_OVERFLOW::FAK_OVERFLOW() throw() { fkt="<unknown function>"; }
FAK_OVERFLOW::FAK_OVERFLOW(const string &f) throw() { fkt=f; }

//member functions of class USE_OF_UNINITIALIZED_OBJ
int USE_OF_UNINITIALIZED_OBJ::errnum() const throw() { return 208; }
string USE_OF_UNINITIALIZED_OBJ::errtext() const throw() { return fkt+": USE_OF_UNINITIALIZED_OBJ"; }
USE_OF_UNINITIALIZED_OBJ::USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
USE_OF_UNINITIALIZED_OBJ::USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//------------- Occuring Errors --------------------------

//member functions of template class ERROR__ELEMENT_NOT_IN_VEC
template <class T>
int ERROR__ELEMENT_NOT_IN_VEC<T>::errnum() const throw() { return 7203; }
template <class T>
string ERROR__ELEMENT_NOT_IN_VEC<T>::errtext() const throw() { return fkt+": ERROR__ELEMENT_NOT_IN_VEC"; }
template <class T>
ERROR__ELEMENT_NOT_IN_VEC<T>::ERROR__ELEMENT_NOT_IN_VEC() throw() { fkt="<unknown function>"; }
template <class T>
ERROR__ELEMENT_NOT_IN_VEC<T>::ERROR__ELEMENT_NOT_IN_VEC(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ
template <class T>
int ERROR__TYPE_CAST_OF_THICK_OBJ<T>::errnum() const throw() { return 11206; }
template <class T>
string ERROR__TYPE_CAST_OF_THICK_OBJ<T>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
template <class T>
ERROR__TYPE_CAST_OF_THICK_OBJ<T>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
template <class T>
ERROR__TYPE_CAST_OF_THICK_OBJ<T>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ
template <class T>
int ERROR__USE_OF_UNINITIALIZED_OBJ<T>::errnum() const throw() { return 7208; }
template <class T>
string ERROR__USE_OF_UNINITIALIZED_OBJ<T>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
template <class T>
ERROR__USE_OF_UNINITIALIZED_OBJ<T>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
template <class T>
ERROR__USE_OF_UNINITIALIZED_OBJ<T>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM
template <class T>
int ERROR__OP_WITH_WRONG_DIM<T>::errnum() const throw() { return 11207; }
template <class T>
string ERROR__OP_WITH_WRONG_DIM<T>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
template <class T>
ERROR__OP_WITH_WRONG_DIM<T>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
template <class T>
ERROR__OP_WITH_WRONG_DIM<T>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES
template <class T>
int ERROR__WRONG_BOUNDARIES<T>::errnum() const throw() { return 11200; }
template <class T>
string ERROR__WRONG_BOUNDARIES<T>::errtext() const throw() { return fkt+": ERROR__WRONG_BOUNDARIES"; }
template <class T>
ERROR__WRONG_BOUNDARIES<T>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
template <class T>
ERROR__WRONG_BOUNDARIES<T>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY
template <class T>
int ERROR__NO_MORE_MEMORY<T>::errnum() const throw() { return 11002; }
template <class T>
string ERROR__NO_MORE_MEMORY<T>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
template <class T>
ERROR__NO_MORE_MEMORY<T>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
template <class T>
ERROR__NO_MORE_MEMORY<T>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG
template <class T>
int ERROR__SUB_ARRAY_TOO_BIG<T>::errnum() const throw() { return 11201; }
template <class T>
string ERROR__SUB_ARRAY_TOO_BIG<T>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
template <class T>
ERROR__SUB_ARRAY_TOO_BIG<T>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
template <class T>
ERROR__SUB_ARRAY_TOO_BIG<T>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ROW_OR_COL_NOT_IN_MAT
template <class T>
int ERROR__ROW_OR_COL_NOT_IN_MAT<T>::errnum() const throw() { return 11204; }
template <class T>
string ERROR__ROW_OR_COL_NOT_IN_MAT<T>::errtext() const throw() { return fkt+": ERROR__ROW_OR_COL_NOT_IN_MAT"; }
template <class T>
ERROR__ROW_OR_COL_NOT_IN_MAT<T>::ERROR__ROW_OR_COL_NOT_IN_MAT() throw() { fkt="<unknown function>"; }
template <class T>
ERROR__ROW_OR_COL_NOT_IN_MAT<T>::ERROR__ROW_OR_COL_NOT_IN_MAT(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_ROW_OR_COL
template <class T>
int ERROR__WRONG_ROW_OR_COL<T>::errnum() const throw() { return 11205; }
template <class T>
string ERROR__WRONG_ROW_OR_COL<T>::errtext() const throw() { return fkt+": ERROR__WRONG_ROW_OR_COL"; }
template <class T>
ERROR__WRONG_ROW_OR_COL<T>::ERROR__WRONG_ROW_OR_COL() throw() { fkt="<unknown function>"; }
template <class T>
ERROR__WRONG_ROW_OR_COL<T>::ERROR__WRONG_ROW_OR_COL(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<rvector>
int ERROR__WRONG_BOUNDARIES<rvector>::errnum() const throw() { return 7200; }
string ERROR__WRONG_BOUNDARIES<rvector>::errtext() const throw() { return fkt+": ERROR_RVECTOR_WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<rvector>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<rvector>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<rvector>
int ERROR__NO_MORE_MEMORY<rvector>::errnum() const throw() { return 7002; }
string ERROR__NO_MORE_MEMORY<rvector>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<rvector>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<rvector>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<rvector>
int ERROR__SUB_ARRAY_TOO_BIG<rvector>::errnum() const throw() { return 7201; }
string ERROR__SUB_ARRAY_TOO_BIG<rvector>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<rvector>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<rvector>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ELEMENT_NOT_IN_VEC<rvector>
int ERROR__ELEMENT_NOT_IN_VEC<rvector>::errnum() const throw() { return 7203; }
string ERROR__ELEMENT_NOT_IN_VEC<rvector>::errtext() const throw() { return fkt+": ERROR__ELEMENT_NOT_IN_VEC"; }
ERROR__ELEMENT_NOT_IN_VEC<rvector>::ERROR__ELEMENT_NOT_IN_VEC() throw() { fkt="<unknown function>"; }
ERROR__ELEMENT_NOT_IN_VEC<rvector>::ERROR__ELEMENT_NOT_IN_VEC(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<rvector>
int ERROR__TYPE_CAST_OF_THICK_OBJ<rvector>::errnum() const throw() { return 7206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<rvector>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<rvector>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<rvector>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<rmatrix>
int ERROR__TYPE_CAST_OF_THICK_OBJ<rmatrix>::errnum() const throw() { return 11206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<rmatrix>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<rmatrix>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<rmatrix>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<rvector>
int ERROR__USE_OF_UNINITIALIZED_OBJ<rvector>::errnum() const throw() { return 7208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<rvector>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<rvector>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<rvector>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<rmatrix>
int ERROR__USE_OF_UNINITIALIZED_OBJ<rmatrix>::errnum() const throw() { return 11208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<rmatrix>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<rmatrix>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<rmatrix>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<rvector>
int ERROR__OP_WITH_WRONG_DIM<rvector>::errnum() const throw() { return 7207; }
string ERROR__OP_WITH_WRONG_DIM<rvector>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<rvector>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<rvector>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<rmatrix>
int ERROR__OP_WITH_WRONG_DIM<rmatrix>::errnum() const throw() { return 11207; }
string ERROR__OP_WITH_WRONG_DIM<rmatrix>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<rmatrix>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<rmatrix>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<rmatrix>
int ERROR__WRONG_BOUNDARIES<rmatrix>::errnum() const throw() { return 11200; }
string ERROR__WRONG_BOUNDARIES<rmatrix>::errtext() const throw() { return fkt+": ERROR__WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<rmatrix>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<rmatrix>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<rmatrix>
int ERROR__NO_MORE_MEMORY<rmatrix>::errnum() const throw() { return 11002; }
string ERROR__NO_MORE_MEMORY<rmatrix>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<rmatrix>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<rmatrix>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<rmatrix>
int ERROR__SUB_ARRAY_TOO_BIG<rmatrix>::errnum() const throw() { return 11201; }
string ERROR__SUB_ARRAY_TOO_BIG<rmatrix>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<rmatrix>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<rmatrix>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ROW_OR_COL_NOT_IN_MAT<rmatrix>
int ERROR__ROW_OR_COL_NOT_IN_MAT<rmatrix>::errnum() const throw() { return 11204; }
string ERROR__ROW_OR_COL_NOT_IN_MAT<rmatrix>::errtext() const throw() { return fkt+": ERROR__ROW_OR_COL_NOT_IN_MAT"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<rmatrix>::ERROR__ROW_OR_COL_NOT_IN_MAT() throw() { fkt="<unknown function>"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<rmatrix>::ERROR__ROW_OR_COL_NOT_IN_MAT(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_ROW_OR_COL<rmatrix>
int ERROR__WRONG_ROW_OR_COL<rmatrix>::errnum() const throw() { return 11205; }
string ERROR__WRONG_ROW_OR_COL<rmatrix>::errtext() const throw() { return fkt+": ERROR__WRONG_ROW_OR_COL"; }
ERROR__WRONG_ROW_OR_COL<rmatrix>::ERROR__WRONG_ROW_OR_COL() throw() { fkt="<unknown function>"; }
ERROR__WRONG_ROW_OR_COL<rmatrix>::ERROR__WRONG_ROW_OR_COL(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<ivector>
int ERROR__WRONG_BOUNDARIES<ivector>::errnum() const throw() { return 8200; }
string ERROR__WRONG_BOUNDARIES<ivector>::errtext() const throw() { return fkt+": ERROR__WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<ivector>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<ivector>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<ivector>
int ERROR__NO_MORE_MEMORY<ivector>::errnum() const throw() { return 8002; }
string ERROR__NO_MORE_MEMORY<ivector>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<ivector>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<ivector>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<ivector>
int ERROR__SUB_ARRAY_TOO_BIG<ivector>::errnum() const throw() { return 8201; }
string ERROR__SUB_ARRAY_TOO_BIG<ivector>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<ivector>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<ivector>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ELEMENT_NOT_IN_VEC<ivector>
int ERROR__ELEMENT_NOT_IN_VEC<ivector>::errnum() const throw() { return 8203; }
string ERROR__ELEMENT_NOT_IN_VEC<ivector>::errtext() const throw() { return fkt+": ERROR__ELEMENT_NOT_IN_VEC"; }
ERROR__ELEMENT_NOT_IN_VEC<ivector>::ERROR__ELEMENT_NOT_IN_VEC() throw() { fkt="<unknown function>"; }
ERROR__ELEMENT_NOT_IN_VEC<ivector>::ERROR__ELEMENT_NOT_IN_VEC(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<ivector>
int ERROR__TYPE_CAST_OF_THICK_OBJ<ivector>::errnum() const throw() { return 8206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<ivector>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<ivector>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<ivector>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<imatrix>
int ERROR__TYPE_CAST_OF_THICK_OBJ<imatrix>::errnum() const throw() { return 12206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<imatrix>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<imatrix>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<imatrix>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<ivector>
int ERROR__USE_OF_UNINITIALIZED_OBJ<ivector>::errnum() const throw() { return 8208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<ivector>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<ivector>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<ivector>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<imatrix>
int ERROR__USE_OF_UNINITIALIZED_OBJ<imatrix>::errnum() const throw() { return 12208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<imatrix>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<imatrix>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<imatrix>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<ivector>
int ERROR__OP_WITH_WRONG_DIM<ivector>::errnum() const throw() { return 8207; }
string ERROR__OP_WITH_WRONG_DIM<ivector>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<ivector>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<ivector>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<imatrix>
int ERROR__OP_WITH_WRONG_DIM<imatrix>::errnum() const throw() { return 12207; }
string ERROR__OP_WITH_WRONG_DIM<imatrix>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<imatrix>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<imatrix>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<imatrix>
int ERROR__WRONG_BOUNDARIES<imatrix>::errnum() const throw() { return 12200; }
string ERROR__WRONG_BOUNDARIES<imatrix>::errtext() const throw() { return fkt+": ERROR__WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<imatrix>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<imatrix>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<imatrix>
int ERROR__NO_MORE_MEMORY<imatrix>::errnum() const throw() { return 12002; }
string ERROR__NO_MORE_MEMORY<imatrix>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<imatrix>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<imatrix>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<imatrix>
int ERROR__SUB_ARRAY_TOO_BIG<imatrix>::errnum() const throw() { return 12201; }
string ERROR__SUB_ARRAY_TOO_BIG<imatrix>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<imatrix>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<imatrix>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ROW_OR_COL_NOT_IN_MAT<imatrix>
int ERROR__ROW_OR_COL_NOT_IN_MAT<imatrix>::errnum() const throw() { return 12204; }
string ERROR__ROW_OR_COL_NOT_IN_MAT<imatrix>::errtext() const throw() { return fkt+": ERROR__ROW_OR_COL_NOT_IN_MAT"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<imatrix>::ERROR__ROW_OR_COL_NOT_IN_MAT() throw() { fkt="<unknown function>"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<imatrix>::ERROR__ROW_OR_COL_NOT_IN_MAT(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_ROW_OR_COL<imatrix>
int ERROR__WRONG_ROW_OR_COL<imatrix>::errnum() const throw() { return 12205; }
string ERROR__WRONG_ROW_OR_COL<imatrix>::errtext() const throw() { return fkt+": ERROR__WRONG_ROW_OR_COL"; }
ERROR__WRONG_ROW_OR_COL<imatrix>::ERROR__WRONG_ROW_OR_COL() throw() { fkt="<unknown function>"; }
ERROR__WRONG_ROW_OR_COL<imatrix>::ERROR__WRONG_ROW_OR_COL(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<rvector_slice>
int ERROR__WRONG_BOUNDARIES<rvector_slice>::errnum() const throw() { return 7200; }
string ERROR__WRONG_BOUNDARIES<rvector_slice>::errtext() const throw() { return fkt+": ERROR_RVECTOR_WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<rvector_slice>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<rvector_slice>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<rvector_slice>
int ERROR__NO_MORE_MEMORY<rvector_slice>::errnum() const throw() { return 7002; }
string ERROR__NO_MORE_MEMORY<rvector_slice>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<rvector_slice>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<rvector_slice>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<rvector_slice>
int ERROR__SUB_ARRAY_TOO_BIG<rvector_slice>::errnum() const throw() { return 7201; }
string ERROR__SUB_ARRAY_TOO_BIG<rvector_slice>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<rvector_slice>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<rvector_slice>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ELEMENT_NOT_IN_VEC<rvector_slice>
int ERROR__ELEMENT_NOT_IN_VEC<rvector_slice>::errnum() const throw() { return 7203; }
string ERROR__ELEMENT_NOT_IN_VEC<rvector_slice>::errtext() const throw() { return fkt+": ERROR__ELEMENT_NOT_IN_VEC"; }
ERROR__ELEMENT_NOT_IN_VEC<rvector_slice>::ERROR__ELEMENT_NOT_IN_VEC() throw() { fkt="<unknown function>"; }
ERROR__ELEMENT_NOT_IN_VEC<rvector_slice>::ERROR__ELEMENT_NOT_IN_VEC(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<rvector_slice>
int ERROR__TYPE_CAST_OF_THICK_OBJ<rvector_slice>::errnum() const throw() { return 7206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<rvector_slice>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<rvector_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<rvector_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<rmatrix_slice>
int ERROR__TYPE_CAST_OF_THICK_OBJ<rmatrix_slice>::errnum() const throw() { return 11206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<rmatrix_slice>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<rmatrix_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<rmatrix_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<rvector_slice>
int ERROR__USE_OF_UNINITIALIZED_OBJ<rvector_slice>::errnum() const throw() { return 7208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<rvector_slice>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<rvector_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<rvector_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<rmatrix_slice>
int ERROR__USE_OF_UNINITIALIZED_OBJ<rmatrix_slice>::errnum() const throw() { return 11208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<rmatrix_slice>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<rmatrix_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<rmatrix_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<rvector_slice>
int ERROR__OP_WITH_WRONG_DIM<rvector_slice>::errnum() const throw() { return 7207; }
string ERROR__OP_WITH_WRONG_DIM<rvector_slice>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<rvector_slice>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<rvector_slice>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<rmatrix_slice>
int ERROR__OP_WITH_WRONG_DIM<rmatrix_slice>::errnum() const throw() { return 11207; }
string ERROR__OP_WITH_WRONG_DIM<rmatrix_slice>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<rmatrix_slice>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<rmatrix_slice>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<rmatrix_slice>
int ERROR__WRONG_BOUNDARIES<rmatrix_slice>::errnum() const throw() { return 11200; }
string ERROR__WRONG_BOUNDARIES<rmatrix_slice>::errtext() const throw() { return fkt+": ERROR__WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<rmatrix_slice>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<rmatrix_slice>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<rmatrix_slice>
int ERROR__NO_MORE_MEMORY<rmatrix_slice>::errnum() const throw() { return 11002; }
string ERROR__NO_MORE_MEMORY<rmatrix_slice>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<rmatrix_slice>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<rmatrix_slice>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<rmatrix_slice>
int ERROR__SUB_ARRAY_TOO_BIG<rmatrix_slice>::errnum() const throw() { return 11201; }
string ERROR__SUB_ARRAY_TOO_BIG<rmatrix_slice>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<rmatrix_slice>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<rmatrix_slice>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ROW_OR_COL_NOT_IN_MAT<rmatrix_slice>
int ERROR__ROW_OR_COL_NOT_IN_MAT<rmatrix_slice>::errnum() const throw() { return 11204; }
string ERROR__ROW_OR_COL_NOT_IN_MAT<rmatrix_slice>::errtext() const throw() { return fkt+": ERROR__ROW_OR_COL_NOT_IN_MAT"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<rmatrix_slice>::ERROR__ROW_OR_COL_NOT_IN_MAT() throw() { fkt="<unknown function>"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<rmatrix_slice>::ERROR__ROW_OR_COL_NOT_IN_MAT(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_ROW_OR_COL<rmatrix_slice>
int ERROR__WRONG_ROW_OR_COL<rmatrix_slice>::errnum() const throw() { return 11205; }
string ERROR__WRONG_ROW_OR_COL<rmatrix_slice>::errtext() const throw() { return fkt+": ERROR__WRONG_ROW_OR_COL"; }
ERROR__WRONG_ROW_OR_COL<rmatrix_slice>::ERROR__WRONG_ROW_OR_COL() throw() { fkt="<unknown function>"; }
ERROR__WRONG_ROW_OR_COL<rmatrix_slice>::ERROR__WRONG_ROW_OR_COL(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<ivector_slice>
int ERROR__WRONG_BOUNDARIES<ivector_slice>::errnum() const throw() { return 8200; }
string ERROR__WRONG_BOUNDARIES<ivector_slice>::errtext() const throw() { return fkt+": ERROR__WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<ivector_slice>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<ivector_slice>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<ivector_slice>
int ERROR__NO_MORE_MEMORY<ivector_slice>::errnum() const throw() { return 8002; }
string ERROR__NO_MORE_MEMORY<ivector_slice>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<ivector_slice>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<ivector_slice>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<ivector_slice>
int ERROR__SUB_ARRAY_TOO_BIG<ivector_slice>::errnum() const throw() { return 8201; }
string ERROR__SUB_ARRAY_TOO_BIG<ivector_slice>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<ivector_slice>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<ivector_slice>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ELEMENT_NOT_IN_VEC<ivector_slice>
int ERROR__ELEMENT_NOT_IN_VEC<ivector_slice>::errnum() const throw() { return 8203; }
string ERROR__ELEMENT_NOT_IN_VEC<ivector_slice>::errtext() const throw() { return fkt+": ERROR__ELEMENT_NOT_IN_VEC"; }
ERROR__ELEMENT_NOT_IN_VEC<ivector_slice>::ERROR__ELEMENT_NOT_IN_VEC() throw() { fkt="<unknown function>"; }
ERROR__ELEMENT_NOT_IN_VEC<ivector_slice>::ERROR__ELEMENT_NOT_IN_VEC(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<ivector_slice>
int ERROR__TYPE_CAST_OF_THICK_OBJ<ivector_slice>::errnum() const throw() { return 8206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<ivector_slice>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<ivector_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<ivector_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<imatrix_slice>
int ERROR__TYPE_CAST_OF_THICK_OBJ<imatrix_slice>::errnum() const throw() { return 12206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<imatrix_slice>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<imatrix_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<imatrix_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<ivector_slice>
int ERROR__USE_OF_UNINITIALIZED_OBJ<ivector_slice>::errnum() const throw() { return 8208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<ivector_slice>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<ivector_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<ivector_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<imatrix_slice>
int ERROR__USE_OF_UNINITIALIZED_OBJ<imatrix_slice>::errnum() const throw() { return 12208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<imatrix_slice>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<imatrix_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<imatrix_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<ivector_slice>
int ERROR__OP_WITH_WRONG_DIM<ivector_slice>::errnum() const throw() { return 8207; }
string ERROR__OP_WITH_WRONG_DIM<ivector_slice>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<ivector_slice>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<ivector_slice>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<imatrix_slice>
int ERROR__OP_WITH_WRONG_DIM<imatrix_slice>::errnum() const throw() { return 12207; }
string ERROR__OP_WITH_WRONG_DIM<imatrix_slice>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<imatrix_slice>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<imatrix_slice>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<imatrix_slice>
int ERROR__WRONG_BOUNDARIES<imatrix_slice>::errnum() const throw() { return 12200; }
string ERROR__WRONG_BOUNDARIES<imatrix_slice>::errtext() const throw() { return fkt+": ERROR__WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<imatrix_slice>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<imatrix_slice>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<imatrix_slice>
int ERROR__NO_MORE_MEMORY<imatrix_slice>::errnum() const throw() { return 12002; }
string ERROR__NO_MORE_MEMORY<imatrix_slice>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<imatrix_slice>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<imatrix_slice>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<imatrix_slice>
int ERROR__SUB_ARRAY_TOO_BIG<imatrix_slice>::errnum() const throw() { return 12201; }
string ERROR__SUB_ARRAY_TOO_BIG<imatrix_slice>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<imatrix_slice>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<imatrix_slice>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ROW_OR_COL_NOT_IN_MAT<imatrix_slice>
int ERROR__ROW_OR_COL_NOT_IN_MAT<imatrix_slice>::errnum() const throw() { return 12204; }
string ERROR__ROW_OR_COL_NOT_IN_MAT<imatrix_slice>::errtext() const throw() { return fkt+": ERROR__ROW_OR_COL_NOT_IN_MAT"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<imatrix_slice>::ERROR__ROW_OR_COL_NOT_IN_MAT() throw() { fkt="<unknown function>"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<imatrix_slice>::ERROR__ROW_OR_COL_NOT_IN_MAT(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_ROW_OR_COL<imatrix_slice>
int ERROR__WRONG_ROW_OR_COL<imatrix_slice>::errnum() const throw() { return 12205; }
string ERROR__WRONG_ROW_OR_COL<imatrix_slice>::errtext() const throw() { return fkt+": ERROR__WRONG_ROW_OR_COL"; }
ERROR__WRONG_ROW_OR_COL<imatrix_slice>::ERROR__WRONG_ROW_OR_COL() throw() { fkt="<unknown function>"; }
ERROR__WRONG_ROW_OR_COL<imatrix_slice>::ERROR__WRONG_ROW_OR_COL(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<rmatrix_subv>
int ERROR__TYPE_CAST_OF_THICK_OBJ<rmatrix_subv>::errnum() const throw() { return 11206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<rmatrix_subv>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<rmatrix_subv>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<rmatrix_subv>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<rmatrix_subv>
int ERROR__USE_OF_UNINITIALIZED_OBJ<rmatrix_subv>::errnum() const throw() { return 11208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<rmatrix_subv>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<rmatrix_subv>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<rmatrix_subv>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<rmatrix_subv>
int ERROR__OP_WITH_WRONG_DIM<rmatrix_subv>::errnum() const throw() { return 11207; }
string ERROR__OP_WITH_WRONG_DIM<rmatrix_subv>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<rmatrix_subv>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<rmatrix_subv>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<rmatrix_subv>
int ERROR__WRONG_BOUNDARIES<rmatrix_subv>::errnum() const throw() { return 11200; }
string ERROR__WRONG_BOUNDARIES<rmatrix_subv>::errtext() const throw() { return fkt+": ERROR__WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<rmatrix_subv>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<rmatrix_subv>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<rmatrix_subv>
int ERROR__NO_MORE_MEMORY<rmatrix_subv>::errnum() const throw() { return 11002; }
string ERROR__NO_MORE_MEMORY<rmatrix_subv>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<rmatrix_subv>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<rmatrix_subv>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<rmatrix_subv>
int ERROR__SUB_ARRAY_TOO_BIG<rmatrix_subv>::errnum() const throw() { return 11201; }
string ERROR__SUB_ARRAY_TOO_BIG<rmatrix_subv>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<rmatrix_subv>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<rmatrix_subv>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ROW_OR_COL_NOT_IN_MAT<rmatrix_subv>
int ERROR__ROW_OR_COL_NOT_IN_MAT<rmatrix_subv>::errnum() const throw() { return 11204; }
string ERROR__ROW_OR_COL_NOT_IN_MAT<rmatrix_subv>::errtext() const throw() { return fkt+": ERROR__ROW_OR_COL_NOT_IN_MAT"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<rmatrix_subv>::ERROR__ROW_OR_COL_NOT_IN_MAT() throw() { fkt="<unknown function>"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<rmatrix_subv>::ERROR__ROW_OR_COL_NOT_IN_MAT(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_ROW_OR_COL<rmatrix_subv>
int ERROR__WRONG_ROW_OR_COL<rmatrix_subv>::errnum() const throw() { return 11205; }
string ERROR__WRONG_ROW_OR_COL<rmatrix_subv>::errtext() const throw() { return fkt+": ERROR__WRONG_ROW_OR_COL"; }
ERROR__WRONG_ROW_OR_COL<rmatrix_subv>::ERROR__WRONG_ROW_OR_COL() throw() { fkt="<unknown function>"; }
ERROR__WRONG_ROW_OR_COL<rmatrix_subv>::ERROR__WRONG_ROW_OR_COL(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<imatrix_subv>
int ERROR__TYPE_CAST_OF_THICK_OBJ<imatrix_subv>::errnum() const throw() { return 12206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<imatrix_subv>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<imatrix_subv>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<imatrix_subv>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<imatrix_subv>
int ERROR__USE_OF_UNINITIALIZED_OBJ<imatrix_subv>::errnum() const throw() { return 12208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<imatrix_subv>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<imatrix_subv>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<imatrix_subv>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<imatrix_subv>
int ERROR__OP_WITH_WRONG_DIM<imatrix_subv>::errnum() const throw() { return 12207; }
string ERROR__OP_WITH_WRONG_DIM<imatrix_subv>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<imatrix_subv>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<imatrix_subv>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<imatrix_subv>
int ERROR__WRONG_BOUNDARIES<imatrix_subv>::errnum() const throw() { return 12200; }
string ERROR__WRONG_BOUNDARIES<imatrix_subv>::errtext() const throw() { return fkt+": ERROR__WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<imatrix_subv>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<imatrix_subv>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<imatrix_subv>
int ERROR__NO_MORE_MEMORY<imatrix_subv>::errnum() const throw() { return 12002; }
string ERROR__NO_MORE_MEMORY<imatrix_subv>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<imatrix_subv>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<imatrix_subv>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<imatrix_subv>
int ERROR__SUB_ARRAY_TOO_BIG<imatrix_subv>::errnum() const throw() { return 12201; }
string ERROR__SUB_ARRAY_TOO_BIG<imatrix_subv>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<imatrix_subv>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<imatrix_subv>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ROW_OR_COL_NOT_IN_MAT<imatrix_subv>
int ERROR__ROW_OR_COL_NOT_IN_MAT<imatrix_subv>::errnum() const throw() { return 12204; }
string ERROR__ROW_OR_COL_NOT_IN_MAT<imatrix_subv>::errtext() const throw() { return fkt+": ERROR__ROW_OR_COL_NOT_IN_MAT"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<imatrix_subv>::ERROR__ROW_OR_COL_NOT_IN_MAT() throw() { fkt="<unknown function>"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<imatrix_subv>::ERROR__ROW_OR_COL_NOT_IN_MAT(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_ROW_OR_COL<imatrix_subv>
int ERROR__WRONG_ROW_OR_COL<imatrix_subv>::errnum() const throw() { return 12205; }
string ERROR__WRONG_ROW_OR_COL<imatrix_subv>::errtext() const throw() { return fkt+": ERROR__WRONG_ROW_OR_COL"; }
ERROR__WRONG_ROW_OR_COL<imatrix_subv>::ERROR__WRONG_ROW_OR_COL() throw() { fkt="<unknown function>"; }
ERROR__WRONG_ROW_OR_COL<imatrix_subv>::ERROR__WRONG_ROW_OR_COL(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<intvector>
int ERROR__WRONG_BOUNDARIES<intvector>::errnum() const throw() { return 7200; }
string ERROR__WRONG_BOUNDARIES<intvector>::errtext() const throw() { return fkt+": ERROR_INTVECTOR_WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<intvector>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<intvector>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<intvector>
int ERROR__NO_MORE_MEMORY<intvector>::errnum() const throw() { return 7002; }
string ERROR__NO_MORE_MEMORY<intvector>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<intvector>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<intvector>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<intvector>
int ERROR__SUB_ARRAY_TOO_BIG<intvector>::errnum() const throw() { return 7201; }
string ERROR__SUB_ARRAY_TOO_BIG<intvector>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<intvector>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<intvector>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ELEMENT_NOT_IN_VEC<intvector>
int ERROR__ELEMENT_NOT_IN_VEC<intvector>::errnum() const throw() { return 7203; }
string ERROR__ELEMENT_NOT_IN_VEC<intvector>::errtext() const throw() { return fkt+": ERROR__ELEMENT_NOT_IN_VEC"; }
ERROR__ELEMENT_NOT_IN_VEC<intvector>::ERROR__ELEMENT_NOT_IN_VEC() throw() { fkt="<unknown function>"; }
ERROR__ELEMENT_NOT_IN_VEC<intvector>::ERROR__ELEMENT_NOT_IN_VEC(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<intvector>
int ERROR__TYPE_CAST_OF_THICK_OBJ<intvector>::errnum() const throw() { return 7206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<intvector>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<intvector>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<intvector>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<intmatrix>
int ERROR__TYPE_CAST_OF_THICK_OBJ<intmatrix>::errnum() const throw() { return 11206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<intmatrix>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<intmatrix>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<intmatrix>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<intvector>
int ERROR__USE_OF_UNINITIALIZED_OBJ<intvector>::errnum() const throw() { return 7208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<intvector>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<intvector>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<intvector>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<intmatrix>
int ERROR__USE_OF_UNINITIALIZED_OBJ<intmatrix>::errnum() const throw() { return 11208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<intmatrix>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<intmatrix>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<intmatrix>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<intvector>
int ERROR__OP_WITH_WRONG_DIM<intvector>::errnum() const throw() { return 7207; }
string ERROR__OP_WITH_WRONG_DIM<intvector>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<intvector>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<intvector>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<intmatrix>
int ERROR__OP_WITH_WRONG_DIM<intmatrix>::errnum() const throw() { return 11207; }
string ERROR__OP_WITH_WRONG_DIM<intmatrix>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<intmatrix>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<intmatrix>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<intmatrix>
int ERROR__WRONG_BOUNDARIES<intmatrix>::errnum() const throw() { return 11200; }
string ERROR__WRONG_BOUNDARIES<intmatrix>::errtext() const throw() { return fkt+": ERROR__WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<intmatrix>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<intmatrix>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<intmatrix>
int ERROR__NO_MORE_MEMORY<intmatrix>::errnum() const throw() { return 11002; }
string ERROR__NO_MORE_MEMORY<intmatrix>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<intmatrix>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<intmatrix>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<intmatrix>
int ERROR__SUB_ARRAY_TOO_BIG<intmatrix>::errnum() const throw() { return 11201; }
string ERROR__SUB_ARRAY_TOO_BIG<intmatrix>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<intmatrix>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<intmatrix>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ROW_OR_COL_NOT_IN_MAT<intmatrix>
int ERROR__ROW_OR_COL_NOT_IN_MAT<intmatrix>::errnum() const throw() { return 11204; }
string ERROR__ROW_OR_COL_NOT_IN_MAT<intmatrix>::errtext() const throw() { return fkt+": ERROR__ROW_OR_COL_NOT_IN_MAT"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<intmatrix>::ERROR__ROW_OR_COL_NOT_IN_MAT() throw() { fkt="<unknown function>"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<intmatrix>::ERROR__ROW_OR_COL_NOT_IN_MAT(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_ROW_OR_COL<intmatrix>
int ERROR__WRONG_ROW_OR_COL<intmatrix>::errnum() const throw() { return 11205; }
string ERROR__WRONG_ROW_OR_COL<intmatrix>::errtext() const throw() { return fkt+": ERROR__WRONG_ROW_OR_COL"; }
ERROR__WRONG_ROW_OR_COL<intmatrix>::ERROR__WRONG_ROW_OR_COL() throw() { fkt="<unknown function>"; }
ERROR__WRONG_ROW_OR_COL<intmatrix>::ERROR__WRONG_ROW_OR_COL(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<intvector_slice>
int ERROR__WRONG_BOUNDARIES<intvector_slice>::errnum() const throw() { return 7200; }
string ERROR__WRONG_BOUNDARIES<intvector_slice>::errtext() const throw() { return fkt+": ERROR_INTVECTOR_WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<intvector_slice>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<intvector_slice>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<intvector_slice>
int ERROR__NO_MORE_MEMORY<intvector_slice>::errnum() const throw() { return 7002; }
string ERROR__NO_MORE_MEMORY<intvector_slice>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<intvector_slice>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<intvector_slice>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<intvector_slice>
int ERROR__SUB_ARRAY_TOO_BIG<intvector_slice>::errnum() const throw() { return 7201; }
string ERROR__SUB_ARRAY_TOO_BIG<intvector_slice>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<intvector_slice>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<intvector_slice>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ELEMENT_NOT_IN_VEC<intvector_slice>
int ERROR__ELEMENT_NOT_IN_VEC<intvector_slice>::errnum() const throw() { return 7203; }
string ERROR__ELEMENT_NOT_IN_VEC<intvector_slice>::errtext() const throw() { return fkt+": ERROR__ELEMENT_NOT_IN_VEC"; }
ERROR__ELEMENT_NOT_IN_VEC<intvector_slice>::ERROR__ELEMENT_NOT_IN_VEC() throw() { fkt="<unknown function>"; }
ERROR__ELEMENT_NOT_IN_VEC<intvector_slice>::ERROR__ELEMENT_NOT_IN_VEC(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<intvector_slice>
int ERROR__TYPE_CAST_OF_THICK_OBJ<intvector_slice>::errnum() const throw() { return 7206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<intvector_slice>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<intvector_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<intvector_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<intmatrix_slice>
int ERROR__TYPE_CAST_OF_THICK_OBJ<intmatrix_slice>::errnum() const throw() { return 11206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<intmatrix_slice>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<intmatrix_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<intmatrix_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<intvector_slice>
int ERROR__USE_OF_UNINITIALIZED_OBJ<intvector_slice>::errnum() const throw() { return 7208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<intvector_slice>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<intvector_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<intvector_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<intmatrix_slice>
int ERROR__USE_OF_UNINITIALIZED_OBJ<intmatrix_slice>::errnum() const throw() { return 11208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<intmatrix_slice>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<intmatrix_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<intmatrix_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<intvector_slice>
int ERROR__OP_WITH_WRONG_DIM<intvector_slice>::errnum() const throw() { return 7207; }
string ERROR__OP_WITH_WRONG_DIM<intvector_slice>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<intvector_slice>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<intvector_slice>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<intmatrix_slice>
int ERROR__OP_WITH_WRONG_DIM<intmatrix_slice>::errnum() const throw() { return 11207; }
string ERROR__OP_WITH_WRONG_DIM<intmatrix_slice>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<intmatrix_slice>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<intmatrix_slice>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<intmatrix_slice>
int ERROR__WRONG_BOUNDARIES<intmatrix_slice>::errnum() const throw() { return 11200; }
string ERROR__WRONG_BOUNDARIES<intmatrix_slice>::errtext() const throw() { return fkt+": ERROR__WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<intmatrix_slice>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<intmatrix_slice>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<intmatrix_slice>
int ERROR__NO_MORE_MEMORY<intmatrix_slice>::errnum() const throw() { return 11002; }
string ERROR__NO_MORE_MEMORY<intmatrix_slice>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<intmatrix_slice>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<intmatrix_slice>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<intmatrix_slice>
int ERROR__SUB_ARRAY_TOO_BIG<intmatrix_slice>::errnum() const throw() { return 11201; }
string ERROR__SUB_ARRAY_TOO_BIG<intmatrix_slice>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<intmatrix_slice>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<intmatrix_slice>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ROW_OR_COL_NOT_IN_MAT<intmatrix_slice>
int ERROR__ROW_OR_COL_NOT_IN_MAT<intmatrix_slice>::errnum() const throw() { return 11204; }
string ERROR__ROW_OR_COL_NOT_IN_MAT<intmatrix_slice>::errtext() const throw() { return fkt+": ERROR__ROW_OR_COL_NOT_IN_MAT"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<intmatrix_slice>::ERROR__ROW_OR_COL_NOT_IN_MAT() throw() { fkt="<unknown function>"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<intmatrix_slice>::ERROR__ROW_OR_COL_NOT_IN_MAT(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_ROW_OR_COL<intmatrix_slice>
int ERROR__WRONG_ROW_OR_COL<intmatrix_slice>::errnum() const throw() { return 11205; }
string ERROR__WRONG_ROW_OR_COL<intmatrix_slice>::errtext() const throw() { return fkt+": ERROR__WRONG_ROW_OR_COL"; }
ERROR__WRONG_ROW_OR_COL<intmatrix_slice>::ERROR__WRONG_ROW_OR_COL() throw() { fkt="<unknown function>"; }
ERROR__WRONG_ROW_OR_COL<intmatrix_slice>::ERROR__WRONG_ROW_OR_COL(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<intmatrix_subv>
int ERROR__TYPE_CAST_OF_THICK_OBJ<intmatrix_subv>::errnum() const throw() { return 11206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<intmatrix_subv>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<intmatrix_subv>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<intmatrix_subv>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<intmatrix_subv>
int ERROR__USE_OF_UNINITIALIZED_OBJ<intmatrix_subv>::errnum() const throw() { return 11208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<intmatrix_subv>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<intmatrix_subv>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<intmatrix_subv>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<intmatrix_subv>
int ERROR__OP_WITH_WRONG_DIM<intmatrix_subv>::errnum() const throw() { return 11207; }
string ERROR__OP_WITH_WRONG_DIM<intmatrix_subv>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<intmatrix_subv>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<intmatrix_subv>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<intmatrix_subv>
int ERROR__WRONG_BOUNDARIES<intmatrix_subv>::errnum() const throw() { return 11200; }
string ERROR__WRONG_BOUNDARIES<intmatrix_subv>::errtext() const throw() { return fkt+": ERROR__WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<intmatrix_subv>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<intmatrix_subv>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<intmatrix_subv>
int ERROR__NO_MORE_MEMORY<intmatrix_subv>::errnum() const throw() { return 11002; }
string ERROR__NO_MORE_MEMORY<intmatrix_subv>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<intmatrix_subv>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<intmatrix_subv>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<intmatrix_subv>
int ERROR__SUB_ARRAY_TOO_BIG<intmatrix_subv>::errnum() const throw() { return 11201; }
string ERROR__SUB_ARRAY_TOO_BIG<intmatrix_subv>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<intmatrix_subv>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<intmatrix_subv>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ROW_OR_COL_NOT_IN_MAT<intmatrix_subv>
int ERROR__ROW_OR_COL_NOT_IN_MAT<intmatrix_subv>::errnum() const throw() { return 11204; }
string ERROR__ROW_OR_COL_NOT_IN_MAT<intmatrix_subv>::errtext() const throw() { return fkt+": ERROR__ROW_OR_COL_NOT_IN_MAT"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<intmatrix_subv>::ERROR__ROW_OR_COL_NOT_IN_MAT() throw() { fkt="<unknown function>"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<intmatrix_subv>::ERROR__ROW_OR_COL_NOT_IN_MAT(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_ROW_OR_COL<intmatrix_subv>
int ERROR__WRONG_ROW_OR_COL<intmatrix_subv>::errnum() const throw() { return 11205; }
string ERROR__WRONG_ROW_OR_COL<intmatrix_subv>::errtext() const throw() { return fkt+": ERROR__WRONG_ROW_OR_COL"; }
ERROR__WRONG_ROW_OR_COL<intmatrix_subv>::ERROR__WRONG_ROW_OR_COL() throw() { fkt="<unknown function>"; }
ERROR__WRONG_ROW_OR_COL<intmatrix_subv>::ERROR__WRONG_ROW_OR_COL(const string &f) throw() { fkt=f; }

// ===

//member functions of template class ERROR__WRONG_BOUNDARIES<cvector>
int ERROR__WRONG_BOUNDARIES<cvector>::errnum() const throw() { return 7200; }
string ERROR__WRONG_BOUNDARIES<cvector>::errtext() const throw() { return fkt+": ERROR_CVECTOR_WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<cvector>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<cvector>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<cvector>
int ERROR__NO_MORE_MEMORY<cvector>::errnum() const throw() { return 7002; }
string ERROR__NO_MORE_MEMORY<cvector>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<cvector>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<cvector>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<cvector>
int ERROR__SUB_ARRAY_TOO_BIG<cvector>::errnum() const throw() { return 7201; }
string ERROR__SUB_ARRAY_TOO_BIG<cvector>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<cvector>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<cvector>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ELEMENT_NOT_IN_VEC<cvector>
int ERROR__ELEMENT_NOT_IN_VEC<cvector>::errnum() const throw() { return 7203; }
string ERROR__ELEMENT_NOT_IN_VEC<cvector>::errtext() const throw() { return fkt+": ERROR__ELEMENT_NOT_IN_VEC"; }
ERROR__ELEMENT_NOT_IN_VEC<cvector>::ERROR__ELEMENT_NOT_IN_VEC() throw() { fkt="<unknown function>"; }
ERROR__ELEMENT_NOT_IN_VEC<cvector>::ERROR__ELEMENT_NOT_IN_VEC(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<cvector>
int ERROR__TYPE_CAST_OF_THICK_OBJ<cvector>::errnum() const throw() { return 7206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<cvector>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<cvector>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<cvector>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<cmatrix>
int ERROR__TYPE_CAST_OF_THICK_OBJ<cmatrix>::errnum() const throw() { return 11206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<cmatrix>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<cmatrix>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<cmatrix>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<cvector>
int ERROR__USE_OF_UNINITIALIZED_OBJ<cvector>::errnum() const throw() { return 7208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<cvector>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<cvector>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<cvector>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<cmatrix>
int ERROR__USE_OF_UNINITIALIZED_OBJ<cmatrix>::errnum() const throw() { return 11208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<cmatrix>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<cmatrix>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<cmatrix>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<cvector>
int ERROR__OP_WITH_WRONG_DIM<cvector>::errnum() const throw() { return 7207; }
string ERROR__OP_WITH_WRONG_DIM<cvector>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<cvector>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<cvector>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<cmatrix>
int ERROR__OP_WITH_WRONG_DIM<cmatrix>::errnum() const throw() { return 11207; }
string ERROR__OP_WITH_WRONG_DIM<cmatrix>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<cmatrix>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<cmatrix>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<cmatrix>
int ERROR__WRONG_BOUNDARIES<cmatrix>::errnum() const throw() { return 11200; }
string ERROR__WRONG_BOUNDARIES<cmatrix>::errtext() const throw() { return fkt+": ERROR__WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<cmatrix>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<cmatrix>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<cmatrix>
int ERROR__NO_MORE_MEMORY<cmatrix>::errnum() const throw() { return 11002; }
string ERROR__NO_MORE_MEMORY<cmatrix>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<cmatrix>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<cmatrix>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<cmatrix>
int ERROR__SUB_ARRAY_TOO_BIG<cmatrix>::errnum() const throw() { return 11201; }
string ERROR__SUB_ARRAY_TOO_BIG<cmatrix>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<cmatrix>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<cmatrix>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ROW_OR_COL_NOT_IN_MAT<cmatrix>
int ERROR__ROW_OR_COL_NOT_IN_MAT<cmatrix>::errnum() const throw() { return 11204; }
string ERROR__ROW_OR_COL_NOT_IN_MAT<cmatrix>::errtext() const throw() { return fkt+": ERROR__ROW_OR_COL_NOT_IN_MAT"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<cmatrix>::ERROR__ROW_OR_COL_NOT_IN_MAT() throw() { fkt="<unknown function>"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<cmatrix>::ERROR__ROW_OR_COL_NOT_IN_MAT(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_ROW_OR_COL<cmatrix>
int ERROR__WRONG_ROW_OR_COL<cmatrix>::errnum() const throw() { return 11205; }
string ERROR__WRONG_ROW_OR_COL<cmatrix>::errtext() const throw() { return fkt+": ERROR__WRONG_ROW_OR_COL"; }
ERROR__WRONG_ROW_OR_COL<cmatrix>::ERROR__WRONG_ROW_OR_COL() throw() { fkt="<unknown function>"; }
ERROR__WRONG_ROW_OR_COL<cmatrix>::ERROR__WRONG_ROW_OR_COL(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<cvector_slice>
int ERROR__WRONG_BOUNDARIES<cvector_slice>::errnum() const throw() { return 7200; }
string ERROR__WRONG_BOUNDARIES<cvector_slice>::errtext() const throw() { return fkt+": ERROR_CVECTOR_WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<cvector_slice>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<cvector_slice>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<cvector_slice>
int ERROR__NO_MORE_MEMORY<cvector_slice>::errnum() const throw() { return 7002; }
string ERROR__NO_MORE_MEMORY<cvector_slice>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<cvector_slice>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<cvector_slice>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<cvector_slice>
int ERROR__SUB_ARRAY_TOO_BIG<cvector_slice>::errnum() const throw() { return 7201; }
string ERROR__SUB_ARRAY_TOO_BIG<cvector_slice>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<cvector_slice>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<cvector_slice>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ELEMENT_NOT_IN_VEC<cvector_slice>
int ERROR__ELEMENT_NOT_IN_VEC<cvector_slice>::errnum() const throw() { return 7203; }
string ERROR__ELEMENT_NOT_IN_VEC<cvector_slice>::errtext() const throw() { return fkt+": ERROR__ELEMENT_NOT_IN_VEC"; }
ERROR__ELEMENT_NOT_IN_VEC<cvector_slice>::ERROR__ELEMENT_NOT_IN_VEC() throw() { fkt="<unknown function>"; }
ERROR__ELEMENT_NOT_IN_VEC<cvector_slice>::ERROR__ELEMENT_NOT_IN_VEC(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<cvector_slice>
int ERROR__TYPE_CAST_OF_THICK_OBJ<cvector_slice>::errnum() const throw() { return 7206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<cvector_slice>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<cvector_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<cvector_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<cmatrix_slice>
int ERROR__TYPE_CAST_OF_THICK_OBJ<cmatrix_slice>::errnum() const throw() { return 11206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<cmatrix_slice>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<cmatrix_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<cmatrix_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<cvector_slice>
int ERROR__USE_OF_UNINITIALIZED_OBJ<cvector_slice>::errnum() const throw() { return 7208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<cvector_slice>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<cvector_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<cvector_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<cmatrix_slice>
int ERROR__USE_OF_UNINITIALIZED_OBJ<cmatrix_slice>::errnum() const throw() { return 11208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<cmatrix_slice>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<cmatrix_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<cmatrix_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<cvector_slice>
int ERROR__OP_WITH_WRONG_DIM<cvector_slice>::errnum() const throw() { return 7207; }
string ERROR__OP_WITH_WRONG_DIM<cvector_slice>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<cvector_slice>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<cvector_slice>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<cmatrix_slice>
int ERROR__OP_WITH_WRONG_DIM<cmatrix_slice>::errnum() const throw() { return 11207; }
string ERROR__OP_WITH_WRONG_DIM<cmatrix_slice>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<cmatrix_slice>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<cmatrix_slice>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<cmatrix_slice>
int ERROR__WRONG_BOUNDARIES<cmatrix_slice>::errnum() const throw() { return 11200; }
string ERROR__WRONG_BOUNDARIES<cmatrix_slice>::errtext() const throw() { return fkt+": ERROR__WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<cmatrix_slice>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<cmatrix_slice>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<cmatrix_slice>
int ERROR__NO_MORE_MEMORY<cmatrix_slice>::errnum() const throw() { return 11002; }
string ERROR__NO_MORE_MEMORY<cmatrix_slice>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<cmatrix_slice>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<cmatrix_slice>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<cmatrix_slice>
int ERROR__SUB_ARRAY_TOO_BIG<cmatrix_slice>::errnum() const throw() { return 11201; }
string ERROR__SUB_ARRAY_TOO_BIG<cmatrix_slice>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<cmatrix_slice>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<cmatrix_slice>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ROW_OR_COL_NOT_IN_MAT<cmatrix_slice>
int ERROR__ROW_OR_COL_NOT_IN_MAT<cmatrix_slice>::errnum() const throw() { return 11204; }
string ERROR__ROW_OR_COL_NOT_IN_MAT<cmatrix_slice>::errtext() const throw() { return fkt+": ERROR__ROW_OR_COL_NOT_IN_MAT"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<cmatrix_slice>::ERROR__ROW_OR_COL_NOT_IN_MAT() throw() { fkt="<unknown function>"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<cmatrix_slice>::ERROR__ROW_OR_COL_NOT_IN_MAT(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_ROW_OR_COL<cmatrix_slice>
int ERROR__WRONG_ROW_OR_COL<cmatrix_slice>::errnum() const throw() { return 11205; }
string ERROR__WRONG_ROW_OR_COL<cmatrix_slice>::errtext() const throw() { return fkt+": ERROR__WRONG_ROW_OR_COL"; }
ERROR__WRONG_ROW_OR_COL<cmatrix_slice>::ERROR__WRONG_ROW_OR_COL() throw() { fkt="<unknown function>"; }
ERROR__WRONG_ROW_OR_COL<cmatrix_slice>::ERROR__WRONG_ROW_OR_COL(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<cmatrix_subv>
int ERROR__TYPE_CAST_OF_THICK_OBJ<cmatrix_subv>::errnum() const throw() { return 11206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<cmatrix_subv>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<cmatrix_subv>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<cmatrix_subv>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<cmatrix_subv>
int ERROR__USE_OF_UNINITIALIZED_OBJ<cmatrix_subv>::errnum() const throw() { return 11208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<cmatrix_subv>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<cmatrix_subv>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<cmatrix_subv>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<cmatrix_subv>
int ERROR__OP_WITH_WRONG_DIM<cmatrix_subv>::errnum() const throw() { return 11207; }
string ERROR__OP_WITH_WRONG_DIM<cmatrix_subv>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<cmatrix_subv>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<cmatrix_subv>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<cmatrix_subv>
int ERROR__WRONG_BOUNDARIES<cmatrix_subv>::errnum() const throw() { return 11200; }
string ERROR__WRONG_BOUNDARIES<cmatrix_subv>::errtext() const throw() { return fkt+": ERROR__WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<cmatrix_subv>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<cmatrix_subv>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<cmatrix_subv>
int ERROR__NO_MORE_MEMORY<cmatrix_subv>::errnum() const throw() { return 11002; }
string ERROR__NO_MORE_MEMORY<cmatrix_subv>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<cmatrix_subv>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<cmatrix_subv>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<cmatrix_subv>
int ERROR__SUB_ARRAY_TOO_BIG<cmatrix_subv>::errnum() const throw() { return 11201; }
string ERROR__SUB_ARRAY_TOO_BIG<cmatrix_subv>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<cmatrix_subv>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<cmatrix_subv>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ROW_OR_COL_NOT_IN_MAT<cmatrix_subv>
int ERROR__ROW_OR_COL_NOT_IN_MAT<cmatrix_subv>::errnum() const throw() { return 11204; }
string ERROR__ROW_OR_COL_NOT_IN_MAT<cmatrix_subv>::errtext() const throw() { return fkt+": ERROR__ROW_OR_COL_NOT_IN_MAT"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<cmatrix_subv>::ERROR__ROW_OR_COL_NOT_IN_MAT() throw() { fkt="<unknown function>"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<cmatrix_subv>::ERROR__ROW_OR_COL_NOT_IN_MAT(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_ROW_OR_COL<cmatrix_subv>
int ERROR__WRONG_ROW_OR_COL<cmatrix_subv>::errnum() const throw() { return 11205; }
string ERROR__WRONG_ROW_OR_COL<cmatrix_subv>::errtext() const throw() { return fkt+": ERROR__WRONG_ROW_OR_COL"; }
ERROR__WRONG_ROW_OR_COL<cmatrix_subv>::ERROR__WRONG_ROW_OR_COL() throw() { fkt="<unknown function>"; }
ERROR__WRONG_ROW_OR_COL<cmatrix_subv>::ERROR__WRONG_ROW_OR_COL(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<civector>
int ERROR__WRONG_BOUNDARIES<civector>::errnum() const throw() { return 7200; }
string ERROR__WRONG_BOUNDARIES<civector>::errtext() const throw() { return fkt+": ERROR_CIVECTOR_WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<civector>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<civector>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<civector>
int ERROR__NO_MORE_MEMORY<civector>::errnum() const throw() { return 7002; }
string ERROR__NO_MORE_MEMORY<civector>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<civector>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<civector>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<civector>
int ERROR__SUB_ARRAY_TOO_BIG<civector>::errnum() const throw() { return 7201; }
string ERROR__SUB_ARRAY_TOO_BIG<civector>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<civector>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<civector>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ELEMENT_NOT_IN_VEC<civector>
int ERROR__ELEMENT_NOT_IN_VEC<civector>::errnum() const throw() { return 7203; }
string ERROR__ELEMENT_NOT_IN_VEC<civector>::errtext() const throw() { return fkt+": ERROR__ELEMENT_NOT_IN_VEC"; }
ERROR__ELEMENT_NOT_IN_VEC<civector>::ERROR__ELEMENT_NOT_IN_VEC() throw() { fkt="<unknown function>"; }
ERROR__ELEMENT_NOT_IN_VEC<civector>::ERROR__ELEMENT_NOT_IN_VEC(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<civector>
int ERROR__TYPE_CAST_OF_THICK_OBJ<civector>::errnum() const throw() { return 7206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<civector>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<civector>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<civector>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<cimatrix>
int ERROR__TYPE_CAST_OF_THICK_OBJ<cimatrix>::errnum() const throw() { return 11206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<cimatrix>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<cimatrix>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<cimatrix>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<civector>
int ERROR__USE_OF_UNINITIALIZED_OBJ<civector>::errnum() const throw() { return 7208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<civector>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<civector>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<civector>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<cimatrix>
int ERROR__USE_OF_UNINITIALIZED_OBJ<cimatrix>::errnum() const throw() { return 11208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<cimatrix>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<cimatrix>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<cimatrix>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<civector>
int ERROR__OP_WITH_WRONG_DIM<civector>::errnum() const throw() { return 7207; }
string ERROR__OP_WITH_WRONG_DIM<civector>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<civector>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<civector>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<cimatrix>
int ERROR__OP_WITH_WRONG_DIM<cimatrix>::errnum() const throw() { return 11207; }
string ERROR__OP_WITH_WRONG_DIM<cimatrix>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<cimatrix>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<cimatrix>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<cimatrix>
int ERROR__WRONG_BOUNDARIES<cimatrix>::errnum() const throw() { return 11200; }
string ERROR__WRONG_BOUNDARIES<cimatrix>::errtext() const throw() { return fkt+": ERROR__WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<cimatrix>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<cimatrix>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<cimatrix>
int ERROR__NO_MORE_MEMORY<cimatrix>::errnum() const throw() { return 11002; }
string ERROR__NO_MORE_MEMORY<cimatrix>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<cimatrix>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<cimatrix>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<cimatrix>
int ERROR__SUB_ARRAY_TOO_BIG<cimatrix>::errnum() const throw() { return 11201; }
string ERROR__SUB_ARRAY_TOO_BIG<cimatrix>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<cimatrix>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<cimatrix>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ROW_OR_COL_NOT_IN_MAT<cimatrix>
int ERROR__ROW_OR_COL_NOT_IN_MAT<cimatrix>::errnum() const throw() { return 11204; }
string ERROR__ROW_OR_COL_NOT_IN_MAT<cimatrix>::errtext() const throw() { return fkt+": ERROR__ROW_OR_COL_NOT_IN_MAT"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<cimatrix>::ERROR__ROW_OR_COL_NOT_IN_MAT() throw() { fkt="<unknown function>"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<cimatrix>::ERROR__ROW_OR_COL_NOT_IN_MAT(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_ROW_OR_COL<cimatrix>
int ERROR__WRONG_ROW_OR_COL<cimatrix>::errnum() const throw() { return 11205; }
string ERROR__WRONG_ROW_OR_COL<cimatrix>::errtext() const throw() { return fkt+": ERROR__WRONG_ROW_OR_COL"; }
ERROR__WRONG_ROW_OR_COL<cimatrix>::ERROR__WRONG_ROW_OR_COL() throw() { fkt="<unknown function>"; }
ERROR__WRONG_ROW_OR_COL<cimatrix>::ERROR__WRONG_ROW_OR_COL(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<civector_slice>
int ERROR__WRONG_BOUNDARIES<civector_slice>::errnum() const throw() { return 7200; }
string ERROR__WRONG_BOUNDARIES<civector_slice>::errtext() const throw() { return fkt+": ERROR_CIVECTOR_WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<civector_slice>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<civector_slice>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<civector_slice>
int ERROR__NO_MORE_MEMORY<civector_slice>::errnum() const throw() { return 7002; }
string ERROR__NO_MORE_MEMORY<civector_slice>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<civector_slice>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<civector_slice>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<civector_slice>
int ERROR__SUB_ARRAY_TOO_BIG<civector_slice>::errnum() const throw() { return 7201; }
string ERROR__SUB_ARRAY_TOO_BIG<civector_slice>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<civector_slice>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<civector_slice>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ELEMENT_NOT_IN_VEC<civector_slice>
int ERROR__ELEMENT_NOT_IN_VEC<civector_slice>::errnum() const throw() { return 7203; }
string ERROR__ELEMENT_NOT_IN_VEC<civector_slice>::errtext() const throw() { return fkt+": ERROR__ELEMENT_NOT_IN_VEC"; }
ERROR__ELEMENT_NOT_IN_VEC<civector_slice>::ERROR__ELEMENT_NOT_IN_VEC() throw() { fkt="<unknown function>"; }
ERROR__ELEMENT_NOT_IN_VEC<civector_slice>::ERROR__ELEMENT_NOT_IN_VEC(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<civector_slice>
int ERROR__TYPE_CAST_OF_THICK_OBJ<civector_slice>::errnum() const throw() { return 7206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<civector_slice>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<civector_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<civector_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<cimatrix_slice>
int ERROR__TYPE_CAST_OF_THICK_OBJ<cimatrix_slice>::errnum() const throw() { return 11206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<cimatrix_slice>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<cimatrix_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<cimatrix_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<civector_slice>
int ERROR__USE_OF_UNINITIALIZED_OBJ<civector_slice>::errnum() const throw() { return 7208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<civector_slice>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<civector_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<civector_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<cimatrix_slice>
int ERROR__USE_OF_UNINITIALIZED_OBJ<cimatrix_slice>::errnum() const throw() { return 11208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<cimatrix_slice>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<cimatrix_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<cimatrix_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<civector_slice>
int ERROR__OP_WITH_WRONG_DIM<civector_slice>::errnum() const throw() { return 7207; }
string ERROR__OP_WITH_WRONG_DIM<civector_slice>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<civector_slice>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<civector_slice>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<cimatrix_slice>
int ERROR__OP_WITH_WRONG_DIM<cimatrix_slice>::errnum() const throw() { return 11207; }
string ERROR__OP_WITH_WRONG_DIM<cimatrix_slice>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<cimatrix_slice>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<cimatrix_slice>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<cimatrix_slice>
int ERROR__WRONG_BOUNDARIES<cimatrix_slice>::errnum() const throw() { return 11200; }
string ERROR__WRONG_BOUNDARIES<cimatrix_slice>::errtext() const throw() { return fkt+": ERROR__WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<cimatrix_slice>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<cimatrix_slice>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<cimatrix_slice>
int ERROR__NO_MORE_MEMORY<cimatrix_slice>::errnum() const throw() { return 11002; }
string ERROR__NO_MORE_MEMORY<cimatrix_slice>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<cimatrix_slice>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<cimatrix_slice>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<cimatrix_slice>
int ERROR__SUB_ARRAY_TOO_BIG<cimatrix_slice>::errnum() const throw() { return 11201; }
string ERROR__SUB_ARRAY_TOO_BIG<cimatrix_slice>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<cimatrix_slice>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<cimatrix_slice>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ROW_OR_COL_NOT_IN_MAT<cimatrix_slice>
int ERROR__ROW_OR_COL_NOT_IN_MAT<cimatrix_slice>::errnum() const throw() { return 11204; }
string ERROR__ROW_OR_COL_NOT_IN_MAT<cimatrix_slice>::errtext() const throw() { return fkt+": ERROR__ROW_OR_COL_NOT_IN_MAT"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<cimatrix_slice>::ERROR__ROW_OR_COL_NOT_IN_MAT() throw() { fkt="<unknown function>"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<cimatrix_slice>::ERROR__ROW_OR_COL_NOT_IN_MAT(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_ROW_OR_COL<cimatrix_slice>
int ERROR__WRONG_ROW_OR_COL<cimatrix_slice>::errnum() const throw() { return 11205; }
string ERROR__WRONG_ROW_OR_COL<cimatrix_slice>::errtext() const throw() { return fkt+": ERROR__WRONG_ROW_OR_COL"; }
ERROR__WRONG_ROW_OR_COL<cimatrix_slice>::ERROR__WRONG_ROW_OR_COL() throw() { fkt="<unknown function>"; }
ERROR__WRONG_ROW_OR_COL<cimatrix_slice>::ERROR__WRONG_ROW_OR_COL(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<cimatrix_subv>
int ERROR__TYPE_CAST_OF_THICK_OBJ<cimatrix_subv>::errnum() const throw() { return 11206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<cimatrix_subv>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<cimatrix_subv>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<cimatrix_subv>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<cimatrix_subv>
int ERROR__USE_OF_UNINITIALIZED_OBJ<cimatrix_subv>::errnum() const throw() { return 11208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<cimatrix_subv>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<cimatrix_subv>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<cimatrix_subv>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<cimatrix_subv>
int ERROR__OP_WITH_WRONG_DIM<cimatrix_subv>::errnum() const throw() { return 11207; }
string ERROR__OP_WITH_WRONG_DIM<cimatrix_subv>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<cimatrix_subv>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<cimatrix_subv>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<cimatrix_subv>
int ERROR__WRONG_BOUNDARIES<cimatrix_subv>::errnum() const throw() { return 11200; }
string ERROR__WRONG_BOUNDARIES<cimatrix_subv>::errtext() const throw() { return fkt+": ERROR__WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<cimatrix_subv>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<cimatrix_subv>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<cimatrix_subv>
int ERROR__NO_MORE_MEMORY<cimatrix_subv>::errnum() const throw() { return 11002; }
string ERROR__NO_MORE_MEMORY<cimatrix_subv>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<cimatrix_subv>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<cimatrix_subv>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<cimatrix_subv>
int ERROR__SUB_ARRAY_TOO_BIG<cimatrix_subv>::errnum() const throw() { return 11201; }
string ERROR__SUB_ARRAY_TOO_BIG<cimatrix_subv>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<cimatrix_subv>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<cimatrix_subv>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ROW_OR_COL_NOT_IN_MAT<cimatrix_subv>
int ERROR__ROW_OR_COL_NOT_IN_MAT<cimatrix_subv>::errnum() const throw() { return 11204; }
string ERROR__ROW_OR_COL_NOT_IN_MAT<cimatrix_subv>::errtext() const throw() { return fkt+": ERROR__ROW_OR_COL_NOT_IN_MAT"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<cimatrix_subv>::ERROR__ROW_OR_COL_NOT_IN_MAT() throw() { fkt="<unknown function>"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<cimatrix_subv>::ERROR__ROW_OR_COL_NOT_IN_MAT(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_ROW_OR_COL<cimatrix_subv>
int ERROR__WRONG_ROW_OR_COL<cimatrix_subv>::errnum() const throw() { return 11205; }
string ERROR__WRONG_ROW_OR_COL<cimatrix_subv>::errtext() const throw() { return fkt+": ERROR__WRONG_ROW_OR_COL"; }
ERROR__WRONG_ROW_OR_COL<cimatrix_subv>::ERROR__WRONG_ROW_OR_COL() throw() { fkt="<unknown function>"; }
ERROR__WRONG_ROW_OR_COL<cimatrix_subv>::ERROR__WRONG_ROW_OR_COL(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<l_rvector>
int ERROR__WRONG_BOUNDARIES<l_rvector>::errnum() const throw() { return 7200; }
string ERROR__WRONG_BOUNDARIES<l_rvector>::errtext() const throw() { return fkt+": ERROR_LRVECTOR_WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<l_rvector>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<l_rvector>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<l_rvector>
int ERROR__NO_MORE_MEMORY<l_rvector>::errnum() const throw() { return 7002; }
string ERROR__NO_MORE_MEMORY<l_rvector>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<l_rvector>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<l_rvector>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<l_rvector>
int ERROR__SUB_ARRAY_TOO_BIG<l_rvector>::errnum() const throw() { return 7201; }
string ERROR__SUB_ARRAY_TOO_BIG<l_rvector>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<l_rvector>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<l_rvector>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ELEMENT_NOT_IN_VEC<l_rvector>
int ERROR__ELEMENT_NOT_IN_VEC<l_rvector>::errnum() const throw() { return 7203; }
string ERROR__ELEMENT_NOT_IN_VEC<l_rvector>::errtext() const throw() { return fkt+": ERROR__ELEMENT_NOT_IN_VEC"; }
ERROR__ELEMENT_NOT_IN_VEC<l_rvector>::ERROR__ELEMENT_NOT_IN_VEC() throw() { fkt="<unknown function>"; }
ERROR__ELEMENT_NOT_IN_VEC<l_rvector>::ERROR__ELEMENT_NOT_IN_VEC(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<l_rvector>
int ERROR__TYPE_CAST_OF_THICK_OBJ<l_rvector>::errnum() const throw() { return 7206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<l_rvector>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<l_rvector>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<l_rvector>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<l_rmatrix>
int ERROR__TYPE_CAST_OF_THICK_OBJ<l_rmatrix>::errnum() const throw() { return 11206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<l_rmatrix>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<l_rmatrix>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<l_rmatrix>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<l_rvector>
int ERROR__USE_OF_UNINITIALIZED_OBJ<l_rvector>::errnum() const throw() { return 7208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<l_rvector>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<l_rvector>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<l_rvector>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<l_rmatrix>
int ERROR__USE_OF_UNINITIALIZED_OBJ<l_rmatrix>::errnum() const throw() { return 11208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<l_rmatrix>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<l_rmatrix>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<l_rmatrix>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<l_rvector>
int ERROR__OP_WITH_WRONG_DIM<l_rvector>::errnum() const throw() { return 7207; }
string ERROR__OP_WITH_WRONG_DIM<l_rvector>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<l_rvector>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<l_rvector>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<l_rmatrix>
int ERROR__OP_WITH_WRONG_DIM<l_rmatrix>::errnum() const throw() { return 11207; }
string ERROR__OP_WITH_WRONG_DIM<l_rmatrix>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<l_rmatrix>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<l_rmatrix>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<l_rmatrix>
int ERROR__WRONG_BOUNDARIES<l_rmatrix>::errnum() const throw() { return 11200; }
string ERROR__WRONG_BOUNDARIES<l_rmatrix>::errtext() const throw() { return fkt+": ERROR__WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<l_rmatrix>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<l_rmatrix>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<l_rmatrix>
int ERROR__NO_MORE_MEMORY<l_rmatrix>::errnum() const throw() { return 11002; }
string ERROR__NO_MORE_MEMORY<l_rmatrix>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<l_rmatrix>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<l_rmatrix>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<l_rmatrix>
int ERROR__SUB_ARRAY_TOO_BIG<l_rmatrix>::errnum() const throw() { return 11201; }
string ERROR__SUB_ARRAY_TOO_BIG<l_rmatrix>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<l_rmatrix>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<l_rmatrix>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ROW_OR_COL_NOT_IN_MAT<l_rmatrix>
int ERROR__ROW_OR_COL_NOT_IN_MAT<l_rmatrix>::errnum() const throw() { return 11204; }
string ERROR__ROW_OR_COL_NOT_IN_MAT<l_rmatrix>::errtext() const throw() { return fkt+": ERROR__ROW_OR_COL_NOT_IN_MAT"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<l_rmatrix>::ERROR__ROW_OR_COL_NOT_IN_MAT() throw() { fkt="<unknown function>"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<l_rmatrix>::ERROR__ROW_OR_COL_NOT_IN_MAT(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_ROW_OR_COL<l_rmatrix>
int ERROR__WRONG_ROW_OR_COL<l_rmatrix>::errnum() const throw() { return 11205; }
string ERROR__WRONG_ROW_OR_COL<l_rmatrix>::errtext() const throw() { return fkt+": ERROR__WRONG_ROW_OR_COL"; }
ERROR__WRONG_ROW_OR_COL<l_rmatrix>::ERROR__WRONG_ROW_OR_COL() throw() { fkt="<unknown function>"; }
ERROR__WRONG_ROW_OR_COL<l_rmatrix>::ERROR__WRONG_ROW_OR_COL(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<l_rvector_slice>
int ERROR__WRONG_BOUNDARIES<l_rvector_slice>::errnum() const throw() { return 7200; }
string ERROR__WRONG_BOUNDARIES<l_rvector_slice>::errtext() const throw() { return fkt+": ERROR_LRVECTOR_WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<l_rvector_slice>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<l_rvector_slice>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<l_rvector_slice>
int ERROR__NO_MORE_MEMORY<l_rvector_slice>::errnum() const throw() { return 7002; }
string ERROR__NO_MORE_MEMORY<l_rvector_slice>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<l_rvector_slice>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<l_rvector_slice>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<l_rvector_slice>
int ERROR__SUB_ARRAY_TOO_BIG<l_rvector_slice>::errnum() const throw() { return 7201; }
string ERROR__SUB_ARRAY_TOO_BIG<l_rvector_slice>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<l_rvector_slice>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<l_rvector_slice>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ELEMENT_NOT_IN_VEC<l_rvector_slice>
int ERROR__ELEMENT_NOT_IN_VEC<l_rvector_slice>::errnum() const throw() { return 7203; }
string ERROR__ELEMENT_NOT_IN_VEC<l_rvector_slice>::errtext() const throw() { return fkt+": ERROR__ELEMENT_NOT_IN_VEC"; }
ERROR__ELEMENT_NOT_IN_VEC<l_rvector_slice>::ERROR__ELEMENT_NOT_IN_VEC() throw() { fkt="<unknown function>"; }
ERROR__ELEMENT_NOT_IN_VEC<l_rvector_slice>::ERROR__ELEMENT_NOT_IN_VEC(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<l_rvector_slice>
int ERROR__TYPE_CAST_OF_THICK_OBJ<l_rvector_slice>::errnum() const throw() { return 7206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<l_rvector_slice>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<l_rvector_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<l_rvector_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<l_rmatrix_slice>
int ERROR__TYPE_CAST_OF_THICK_OBJ<l_rmatrix_slice>::errnum() const throw() { return 11206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<l_rmatrix_slice>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<l_rmatrix_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<l_rmatrix_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<l_rvector_slice>
int ERROR__USE_OF_UNINITIALIZED_OBJ<l_rvector_slice>::errnum() const throw() { return 7208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<l_rvector_slice>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<l_rvector_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<l_rvector_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<l_rmatrix_slice>
int ERROR__USE_OF_UNINITIALIZED_OBJ<l_rmatrix_slice>::errnum() const throw() { return 11208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<l_rmatrix_slice>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<l_rmatrix_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<l_rmatrix_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<l_rvector_slice>
int ERROR__OP_WITH_WRONG_DIM<l_rvector_slice>::errnum() const throw() { return 7207; }
string ERROR__OP_WITH_WRONG_DIM<l_rvector_slice>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<l_rvector_slice>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<l_rvector_slice>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<l_rmatrix_slice>
int ERROR__OP_WITH_WRONG_DIM<l_rmatrix_slice>::errnum() const throw() { return 11207; }
string ERROR__OP_WITH_WRONG_DIM<l_rmatrix_slice>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<l_rmatrix_slice>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<l_rmatrix_slice>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<l_rmatrix_slice>
int ERROR__WRONG_BOUNDARIES<l_rmatrix_slice>::errnum() const throw() { return 11200; }
string ERROR__WRONG_BOUNDARIES<l_rmatrix_slice>::errtext() const throw() { return fkt+": ERROR__WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<l_rmatrix_slice>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<l_rmatrix_slice>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<l_rmatrix_slice>
int ERROR__NO_MORE_MEMORY<l_rmatrix_slice>::errnum() const throw() { return 11002; }
string ERROR__NO_MORE_MEMORY<l_rmatrix_slice>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<l_rmatrix_slice>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<l_rmatrix_slice>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<l_rmatrix_slice>
int ERROR__SUB_ARRAY_TOO_BIG<l_rmatrix_slice>::errnum() const throw() { return 11201; }
string ERROR__SUB_ARRAY_TOO_BIG<l_rmatrix_slice>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<l_rmatrix_slice>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<l_rmatrix_slice>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ROW_OR_COL_NOT_IN_MAT<l_rmatrix_slice>
int ERROR__ROW_OR_COL_NOT_IN_MAT<l_rmatrix_slice>::errnum() const throw() { return 11204; }
string ERROR__ROW_OR_COL_NOT_IN_MAT<l_rmatrix_slice>::errtext() const throw() { return fkt+": ERROR__ROW_OR_COL_NOT_IN_MAT"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<l_rmatrix_slice>::ERROR__ROW_OR_COL_NOT_IN_MAT() throw() { fkt="<unknown function>"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<l_rmatrix_slice>::ERROR__ROW_OR_COL_NOT_IN_MAT(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_ROW_OR_COL<l_rmatrix_slice>
int ERROR__WRONG_ROW_OR_COL<l_rmatrix_slice>::errnum() const throw() { return 11205; }
string ERROR__WRONG_ROW_OR_COL<l_rmatrix_slice>::errtext() const throw() { return fkt+": ERROR__WRONG_ROW_OR_COL"; }
ERROR__WRONG_ROW_OR_COL<l_rmatrix_slice>::ERROR__WRONG_ROW_OR_COL() throw() { fkt="<unknown function>"; }
ERROR__WRONG_ROW_OR_COL<l_rmatrix_slice>::ERROR__WRONG_ROW_OR_COL(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<l_rmatrix_subv>
int ERROR__TYPE_CAST_OF_THICK_OBJ<l_rmatrix_subv>::errnum() const throw() { return 11206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<l_rmatrix_subv>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<l_rmatrix_subv>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<l_rmatrix_subv>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<l_rmatrix_subv>
int ERROR__USE_OF_UNINITIALIZED_OBJ<l_rmatrix_subv>::errnum() const throw() { return 11208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<l_rmatrix_subv>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<l_rmatrix_subv>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<l_rmatrix_subv>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<l_rmatrix_subv>
int ERROR__OP_WITH_WRONG_DIM<l_rmatrix_subv>::errnum() const throw() { return 11207; }
string ERROR__OP_WITH_WRONG_DIM<l_rmatrix_subv>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<l_rmatrix_subv>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<l_rmatrix_subv>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<l_rmatrix_subv>
int ERROR__WRONG_BOUNDARIES<l_rmatrix_subv>::errnum() const throw() { return 11200; }
string ERROR__WRONG_BOUNDARIES<l_rmatrix_subv>::errtext() const throw() { return fkt+": ERROR__WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<l_rmatrix_subv>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<l_rmatrix_subv>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<l_rmatrix_subv>
int ERROR__NO_MORE_MEMORY<l_rmatrix_subv>::errnum() const throw() { return 11002; }
string ERROR__NO_MORE_MEMORY<l_rmatrix_subv>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<l_rmatrix_subv>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<l_rmatrix_subv>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<l_rmatrix_subv>
int ERROR__SUB_ARRAY_TOO_BIG<l_rmatrix_subv>::errnum() const throw() { return 11201; }
string ERROR__SUB_ARRAY_TOO_BIG<l_rmatrix_subv>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<l_rmatrix_subv>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<l_rmatrix_subv>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ROW_OR_COL_NOT_IN_MAT<l_rmatrix_subv>
int ERROR__ROW_OR_COL_NOT_IN_MAT<l_rmatrix_subv>::errnum() const throw() { return 11204; }
string ERROR__ROW_OR_COL_NOT_IN_MAT<l_rmatrix_subv>::errtext() const throw() { return fkt+": ERROR__ROW_OR_COL_NOT_IN_MAT"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<l_rmatrix_subv>::ERROR__ROW_OR_COL_NOT_IN_MAT() throw() { fkt="<unknown function>"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<l_rmatrix_subv>::ERROR__ROW_OR_COL_NOT_IN_MAT(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_ROW_OR_COL<l_rmatrix_subv>
int ERROR__WRONG_ROW_OR_COL<l_rmatrix_subv>::errnum() const throw() { return 11205; }
string ERROR__WRONG_ROW_OR_COL<l_rmatrix_subv>::errtext() const throw() { return fkt+": ERROR__WRONG_ROW_OR_COL"; }
ERROR__WRONG_ROW_OR_COL<l_rmatrix_subv>::ERROR__WRONG_ROW_OR_COL() throw() { fkt="<unknown function>"; }
ERROR__WRONG_ROW_OR_COL<l_rmatrix_subv>::ERROR__WRONG_ROW_OR_COL(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<l_ivector>
int ERROR__WRONG_BOUNDARIES<l_ivector>::errnum() const throw() { return 7200; }
string ERROR__WRONG_BOUNDARIES<l_ivector>::errtext() const throw() { return fkt+": ERROR_LIVECTOR_WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<l_ivector>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<l_ivector>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<l_ivector>
int ERROR__NO_MORE_MEMORY<l_ivector>::errnum() const throw() { return 7002; }
string ERROR__NO_MORE_MEMORY<l_ivector>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<l_ivector>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<l_ivector>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<l_ivector>
int ERROR__SUB_ARRAY_TOO_BIG<l_ivector>::errnum() const throw() { return 7201; }
string ERROR__SUB_ARRAY_TOO_BIG<l_ivector>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<l_ivector>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<l_ivector>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ELEMENT_NOT_IN_VEC<l_ivector>
int ERROR__ELEMENT_NOT_IN_VEC<l_ivector>::errnum() const throw() { return 7203; }
string ERROR__ELEMENT_NOT_IN_VEC<l_ivector>::errtext() const throw() { return fkt+": ERROR__ELEMENT_NOT_IN_VEC"; }
ERROR__ELEMENT_NOT_IN_VEC<l_ivector>::ERROR__ELEMENT_NOT_IN_VEC() throw() { fkt="<unknown function>"; }
ERROR__ELEMENT_NOT_IN_VEC<l_ivector>::ERROR__ELEMENT_NOT_IN_VEC(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<l_ivector>
int ERROR__TYPE_CAST_OF_THICK_OBJ<l_ivector>::errnum() const throw() { return 7206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<l_ivector>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<l_ivector>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<l_ivector>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<l_imatrix>
int ERROR__TYPE_CAST_OF_THICK_OBJ<l_imatrix>::errnum() const throw() { return 11206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<l_imatrix>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<l_imatrix>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<l_imatrix>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<l_ivector>
int ERROR__USE_OF_UNINITIALIZED_OBJ<l_ivector>::errnum() const throw() { return 7208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<l_ivector>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<l_ivector>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<l_ivector>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<l_imatrix>
int ERROR__USE_OF_UNINITIALIZED_OBJ<l_imatrix>::errnum() const throw() { return 11208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<l_imatrix>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<l_imatrix>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<l_imatrix>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<l_ivector>
int ERROR__OP_WITH_WRONG_DIM<l_ivector>::errnum() const throw() { return 7207; }
string ERROR__OP_WITH_WRONG_DIM<l_ivector>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<l_ivector>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<l_ivector>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<l_imatrix>
int ERROR__OP_WITH_WRONG_DIM<l_imatrix>::errnum() const throw() { return 11207; }
string ERROR__OP_WITH_WRONG_DIM<l_imatrix>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<l_imatrix>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<l_imatrix>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<l_imatrix>
int ERROR__WRONG_BOUNDARIES<l_imatrix>::errnum() const throw() { return 11200; }
string ERROR__WRONG_BOUNDARIES<l_imatrix>::errtext() const throw() { return fkt+": ERROR__WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<l_imatrix>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<l_imatrix>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<l_imatrix>
int ERROR__NO_MORE_MEMORY<l_imatrix>::errnum() const throw() { return 11002; }
string ERROR__NO_MORE_MEMORY<l_imatrix>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<l_imatrix>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<l_imatrix>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<l_imatrix>
int ERROR__SUB_ARRAY_TOO_BIG<l_imatrix>::errnum() const throw() { return 11201; }
string ERROR__SUB_ARRAY_TOO_BIG<l_imatrix>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<l_imatrix>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<l_imatrix>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ROW_OR_COL_NOT_IN_MAT<l_imatrix>
int ERROR__ROW_OR_COL_NOT_IN_MAT<l_imatrix>::errnum() const throw() { return 11204; }
string ERROR__ROW_OR_COL_NOT_IN_MAT<l_imatrix>::errtext() const throw() { return fkt+": ERROR__ROW_OR_COL_NOT_IN_MAT"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<l_imatrix>::ERROR__ROW_OR_COL_NOT_IN_MAT() throw() { fkt="<unknown function>"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<l_imatrix>::ERROR__ROW_OR_COL_NOT_IN_MAT(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_ROW_OR_COL<l_imatrix>
int ERROR__WRONG_ROW_OR_COL<l_imatrix>::errnum() const throw() { return 11205; }
string ERROR__WRONG_ROW_OR_COL<l_imatrix>::errtext() const throw() { return fkt+": ERROR__WRONG_ROW_OR_COL"; }
ERROR__WRONG_ROW_OR_COL<l_imatrix>::ERROR__WRONG_ROW_OR_COL() throw() { fkt="<unknown function>"; }
ERROR__WRONG_ROW_OR_COL<l_imatrix>::ERROR__WRONG_ROW_OR_COL(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<l_ivector_slice>
int ERROR__WRONG_BOUNDARIES<l_ivector_slice>::errnum() const throw() { return 7200; }
string ERROR__WRONG_BOUNDARIES<l_ivector_slice>::errtext() const throw() { return fkt+": ERROR_LIVECTOR_WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<l_ivector_slice>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<l_ivector_slice>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<l_ivector_slice>
int ERROR__NO_MORE_MEMORY<l_ivector_slice>::errnum() const throw() { return 7002; }
string ERROR__NO_MORE_MEMORY<l_ivector_slice>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<l_ivector_slice>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<l_ivector_slice>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<l_ivector_slice>
int ERROR__SUB_ARRAY_TOO_BIG<l_ivector_slice>::errnum() const throw() { return 7201; }
string ERROR__SUB_ARRAY_TOO_BIG<l_ivector_slice>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<l_ivector_slice>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<l_ivector_slice>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ELEMENT_NOT_IN_VEC<l_ivector_slice>
int ERROR__ELEMENT_NOT_IN_VEC<l_ivector_slice>::errnum() const throw() { return 7203; }
string ERROR__ELEMENT_NOT_IN_VEC<l_ivector_slice>::errtext() const throw() { return fkt+": ERROR__ELEMENT_NOT_IN_VEC"; }
ERROR__ELEMENT_NOT_IN_VEC<l_ivector_slice>::ERROR__ELEMENT_NOT_IN_VEC() throw() { fkt="<unknown function>"; }
ERROR__ELEMENT_NOT_IN_VEC<l_ivector_slice>::ERROR__ELEMENT_NOT_IN_VEC(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<l_ivector_slice>
int ERROR__TYPE_CAST_OF_THICK_OBJ<l_ivector_slice>::errnum() const throw() { return 7206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<l_ivector_slice>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<l_ivector_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<l_ivector_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<l_imatrix_slice>
int ERROR__TYPE_CAST_OF_THICK_OBJ<l_imatrix_slice>::errnum() const throw() { return 11206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<l_imatrix_slice>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<l_imatrix_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<l_imatrix_slice>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<l_ivector_slice>
int ERROR__USE_OF_UNINITIALIZED_OBJ<l_ivector_slice>::errnum() const throw() { return 7208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<l_ivector_slice>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<l_ivector_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<l_ivector_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<l_imatrix_slice>
int ERROR__USE_OF_UNINITIALIZED_OBJ<l_imatrix_slice>::errnum() const throw() { return 11208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<l_imatrix_slice>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<l_imatrix_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<l_imatrix_slice>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<l_ivector_slice>
int ERROR__OP_WITH_WRONG_DIM<l_ivector_slice>::errnum() const throw() { return 7207; }
string ERROR__OP_WITH_WRONG_DIM<l_ivector_slice>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<l_ivector_slice>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<l_ivector_slice>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<l_imatrix_slice>
int ERROR__OP_WITH_WRONG_DIM<l_imatrix_slice>::errnum() const throw() { return 11207; }
string ERROR__OP_WITH_WRONG_DIM<l_imatrix_slice>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<l_imatrix_slice>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<l_imatrix_slice>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<l_imatrix_slice>
int ERROR__WRONG_BOUNDARIES<l_imatrix_slice>::errnum() const throw() { return 11200; }
string ERROR__WRONG_BOUNDARIES<l_imatrix_slice>::errtext() const throw() { return fkt+": ERROR__WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<l_imatrix_slice>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<l_imatrix_slice>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<l_imatrix_slice>
int ERROR__NO_MORE_MEMORY<l_imatrix_slice>::errnum() const throw() { return 11002; }
string ERROR__NO_MORE_MEMORY<l_imatrix_slice>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<l_imatrix_slice>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<l_imatrix_slice>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<l_imatrix_slice>
int ERROR__SUB_ARRAY_TOO_BIG<l_imatrix_slice>::errnum() const throw() { return 11201; }
string ERROR__SUB_ARRAY_TOO_BIG<l_imatrix_slice>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<l_imatrix_slice>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<l_imatrix_slice>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ROW_OR_COL_NOT_IN_MAT<l_imatrix_slice>
int ERROR__ROW_OR_COL_NOT_IN_MAT<l_imatrix_slice>::errnum() const throw() { return 11204; }
string ERROR__ROW_OR_COL_NOT_IN_MAT<l_imatrix_slice>::errtext() const throw() { return fkt+": ERROR__ROW_OR_COL_NOT_IN_MAT"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<l_imatrix_slice>::ERROR__ROW_OR_COL_NOT_IN_MAT() throw() { fkt="<unknown function>"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<l_imatrix_slice>::ERROR__ROW_OR_COL_NOT_IN_MAT(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_ROW_OR_COL<l_imatrix_slice>
int ERROR__WRONG_ROW_OR_COL<l_imatrix_slice>::errnum() const throw() { return 11205; }
string ERROR__WRONG_ROW_OR_COL<l_imatrix_slice>::errtext() const throw() { return fkt+": ERROR__WRONG_ROW_OR_COL"; }
ERROR__WRONG_ROW_OR_COL<l_imatrix_slice>::ERROR__WRONG_ROW_OR_COL() throw() { fkt="<unknown function>"; }
ERROR__WRONG_ROW_OR_COL<l_imatrix_slice>::ERROR__WRONG_ROW_OR_COL(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__TYPE_CAST_OF_THICK_OBJ<l_imatrix_subv>
int ERROR__TYPE_CAST_OF_THICK_OBJ<l_imatrix_subv>::errnum() const throw() { return 11206; }
string ERROR__TYPE_CAST_OF_THICK_OBJ<l_imatrix_subv>::errtext() const throw() { return fkt+": ERROR__TYPE_CAST_OF_THICK_OBJ"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<l_imatrix_subv>::ERROR__TYPE_CAST_OF_THICK_OBJ() throw() { fkt="<unknown function>"; }
ERROR__TYPE_CAST_OF_THICK_OBJ<l_imatrix_subv>::ERROR__TYPE_CAST_OF_THICK_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__USE_OF_UNINITIALIZED_OBJ<l_imatrix_subv>
int ERROR__USE_OF_UNINITIALIZED_OBJ<l_imatrix_subv>::errnum() const throw() { return 11208; }
string ERROR__USE_OF_UNINITIALIZED_OBJ<l_imatrix_subv>::errtext() const throw() { return fkt+": ERROR__USE_OF_UNINITIALIZED_OBJ"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<l_imatrix_subv>::ERROR__USE_OF_UNINITIALIZED_OBJ() throw() { fkt="<unknown function>"; }
ERROR__USE_OF_UNINITIALIZED_OBJ<l_imatrix_subv>::ERROR__USE_OF_UNINITIALIZED_OBJ(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__OP_WITH_WRONG_DIM<l_imatrix_subv>
int ERROR__OP_WITH_WRONG_DIM<l_imatrix_subv>::errnum() const throw() { return 11207; }
string ERROR__OP_WITH_WRONG_DIM<l_imatrix_subv>::errtext() const throw() { return fkt+": ERROR__OP_WITH_WRONG_DIM"; }
ERROR__OP_WITH_WRONG_DIM<l_imatrix_subv>::ERROR__OP_WITH_WRONG_DIM() throw() { fkt="<unknown function>"; }
ERROR__OP_WITH_WRONG_DIM<l_imatrix_subv>::ERROR__OP_WITH_WRONG_DIM(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_BOUNDARIES<l_imatrix_subv>
int ERROR__WRONG_BOUNDARIES<l_imatrix_subv>::errnum() const throw() { return 11200; }
string ERROR__WRONG_BOUNDARIES<l_imatrix_subv>::errtext() const throw() { return fkt+": ERROR__WRONG_BOUNDARIES"; }
ERROR__WRONG_BOUNDARIES<l_imatrix_subv>::ERROR__WRONG_BOUNDARIES() throw() { fkt="<unknown function>"; }
ERROR__WRONG_BOUNDARIES<l_imatrix_subv>::ERROR__WRONG_BOUNDARIES(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__NO_MORE_MEMORY<l_imatrix_subv>
int ERROR__NO_MORE_MEMORY<l_imatrix_subv>::errnum() const throw() { return 11002; }
string ERROR__NO_MORE_MEMORY<l_imatrix_subv>::errtext() const throw() { return fkt+": ERROR__NO_MORE_MEMORY"; }
ERROR__NO_MORE_MEMORY<l_imatrix_subv>::ERROR__NO_MORE_MEMORY() throw() { fkt="<unknown function>"; }
ERROR__NO_MORE_MEMORY<l_imatrix_subv>::ERROR__NO_MORE_MEMORY(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__SUB_ARRAY_TOO_BIG<l_imatrix_subv>
int ERROR__SUB_ARRAY_TOO_BIG<l_imatrix_subv>::errnum() const throw() { return 11201; }
string ERROR__SUB_ARRAY_TOO_BIG<l_imatrix_subv>::errtext() const throw() { return fkt+": ERROR__SUB_ARRAY_TOO_BIG"; }
ERROR__SUB_ARRAY_TOO_BIG<l_imatrix_subv>::ERROR__SUB_ARRAY_TOO_BIG() throw() { fkt="<unknown function>"; }
ERROR__SUB_ARRAY_TOO_BIG<l_imatrix_subv>::ERROR__SUB_ARRAY_TOO_BIG(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__ROW_OR_COL_NOT_IN_MAT<l_imatrix_subv>
int ERROR__ROW_OR_COL_NOT_IN_MAT<l_imatrix_subv>::errnum() const throw() { return 11204; }
string ERROR__ROW_OR_COL_NOT_IN_MAT<l_imatrix_subv>::errtext() const throw() { return fkt+": ERROR__ROW_OR_COL_NOT_IN_MAT"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<l_imatrix_subv>::ERROR__ROW_OR_COL_NOT_IN_MAT() throw() { fkt="<unknown function>"; }
ERROR__ROW_OR_COL_NOT_IN_MAT<l_imatrix_subv>::ERROR__ROW_OR_COL_NOT_IN_MAT(const string &f) throw() { fkt=f; }

//member functions of template class ERROR__WRONG_ROW_OR_COL<l_imatrix_subv>
int ERROR__WRONG_ROW_OR_COL<l_imatrix_subv>::errnum() const throw() { return 11205; }
string ERROR__WRONG_ROW_OR_COL<l_imatrix_subv>::errtext() const throw() { return fkt+": ERROR__WRONG_ROW_OR_COL"; }
ERROR__WRONG_ROW_OR_COL<l_imatrix_subv>::ERROR__WRONG_ROW_OR_COL() throw() { fkt="<unknown function>"; }
ERROR__WRONG_ROW_OR_COL<l_imatrix_subv>::ERROR__WRONG_ROW_OR_COL(const string &f) throw() { fkt=f; }


typedef ERROR__WRONG_BOUNDARIES<rvector> ERROR_RVECTOR_WRONG_BOUNDARIES;
typedef ERROR__NO_MORE_MEMORY<rvector> ERROR_RVECTOR_NO_MORE_MEMORY;
typedef ERROR__SUB_ARRAY_TOO_BIG<rvector> ERROR_RVECTOR_SUB_ARRAY_TOO_BIG;
typedef ERROR__ELEMENT_NOT_IN_VEC<rvector> ERROR_RVECTOR_ELEMENT_NOT_IN_VEC;
typedef ERROR__TYPE_CAST_OF_THICK_OBJ<rvector> ERROR_RVECTOR_TYPE_CAST_OF_THICK_OBJ;
typedef ERROR__TYPE_CAST_OF_THICK_OBJ<rmatrix> ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ;
typedef ERROR__USE_OF_UNINITIALIZED_OBJ<rvector> ERROR_RVECTOR_USE_OF_UNINITIALIZED_OBJ;
typedef ERROR__USE_OF_UNINITIALIZED_OBJ<rmatrix> ERROR_RMATRIX_USE_OF_UNINITIALIZED_OBJ;
typedef ERROR__OP_WITH_WRONG_DIM<rvector> ERROR_RVECTOR_OP_WITH_WRONG_DIM;
typedef ERROR__OP_WITH_WRONG_DIM<rmatrix> ERROR_RMATRIX_OP_WITH_WRONG_DIM;
typedef ERROR__WRONG_BOUNDARIES<rmatrix> ERROR_RMATRIX_WRONG_BOUNDARIES;
typedef ERROR__NO_MORE_MEMORY<rmatrix> ERROR_RMATRIX_NO_MORE_MEMORY;
typedef ERROR__SUB_ARRAY_TOO_BIG<rmatrix> ERROR_RMATRIX_SUB_ARRAY_TOO_BIG;
typedef ERROR__ROW_OR_COL_NOT_IN_MAT<rmatrix> ERROR_RMATRIX_ROW_OR_COL_NOT_IN_MAT;
typedef ERROR__WRONG_ROW_OR_COL<rmatrix> ERROR_RMATRIX_WRONG_ROW_OR_COL;
typedef ERROR__WRONG_BOUNDARIES<cvector> ERROR_CVECTOR_WRONG_BOUNDARIES;
typedef ERROR__NO_MORE_MEMORY<cvector> ERROR_CVECTOR_NO_MORE_MEMORY;
typedef ERROR__SUB_ARRAY_TOO_BIG<cvector> ERROR_CVECTOR_SUB_ARRAY_TOO_BIG;
typedef ERROR__ELEMENT_NOT_IN_VEC<cvector> ERROR_CVECTOR_ELEMENT_NOT_IN_VEC;
typedef ERROR__TYPE_CAST_OF_THICK_OBJ<cvector> ERROR_CVECTOR_TYPE_CAST_OF_THICK_OBJ;
typedef ERROR__TYPE_CAST_OF_THICK_OBJ<cmatrix> ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ;
typedef ERROR__USE_OF_UNINITIALIZED_OBJ<cvector> ERROR_CVECTOR_USE_OF_UNINITIALIZED_OBJ;
typedef ERROR__USE_OF_UNINITIALIZED_OBJ<cmatrix> ERROR_CMATRIX_USE_OF_UNINITIALIZED_OBJ;
typedef ERROR__OP_WITH_WRONG_DIM<cvector> ERROR_CVECTOR_OP_WITH_WRONG_DIM;
typedef ERROR__OP_WITH_WRONG_DIM<cmatrix> ERROR_CMATRIX_OP_WITH_WRONG_DIM;
typedef ERROR__WRONG_BOUNDARIES<cmatrix> ERROR_CMATRIX_WRONG_BOUNDARIES;
typedef ERROR__NO_MORE_MEMORY<cmatrix> ERROR_CMATRIX_NO_MORE_MEMORY;
typedef ERROR__SUB_ARRAY_TOO_BIG<cmatrix> ERROR_CMATRIX_SUB_ARRAY_TOO_BIG;
typedef ERROR__ROW_OR_COL_NOT_IN_MAT<cmatrix> ERROR_CMATRIX_ROW_OR_COL_NOT_IN_MAT;
typedef ERROR__WRONG_ROW_OR_COL<cmatrix> ERROR_CMATRIX_WRONG_ROW_OR_COL;
typedef ERROR__WRONG_BOUNDARIES<ivector> ERROR_IVECTOR_WRONG_BOUNDARIES;
typedef ERROR__NO_MORE_MEMORY<ivector> ERROR_IVECTOR_NO_MORE_MEMORY;
typedef ERROR__SUB_ARRAY_TOO_BIG<ivector> ERROR_IVECTOR_SUB_ARRAY_TOO_BIG;
typedef ERROR__ELEMENT_NOT_IN_VEC<ivector> ERROR_IVECTOR_ELEMENT_NOT_IN_VEC;
typedef ERROR__TYPE_CAST_OF_THICK_OBJ<ivector> ERROR_IVECTOR_TYPE_CAST_OF_THICK_OBJ;
typedef ERROR__TYPE_CAST_OF_THICK_OBJ<imatrix> ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ;
typedef ERROR__USE_OF_UNINITIALIZED_OBJ<ivector> ERROR_IVECTOR_USE_OF_UNINITIALIZED_OBJ;
typedef ERROR__USE_OF_UNINITIALIZED_OBJ<imatrix> ERROR_IMATRIX_USE_OF_UNINITIALIZED_OBJ;
typedef ERROR__OP_WITH_WRONG_DIM<ivector> ERROR_IVECTOR_OP_WITH_WRONG_DIM;
typedef ERROR__OP_WITH_WRONG_DIM<imatrix> ERROR_IMATRIX_OP_WITH_WRONG_DIM;
typedef ERROR__WRONG_BOUNDARIES<imatrix> ERROR_IMATRIX_WRONG_BOUNDARIES;
typedef ERROR__NO_MORE_MEMORY<imatrix> ERROR_IMATRIX_NO_MORE_MEMORY;
typedef ERROR__SUB_ARRAY_TOO_BIG<imatrix> ERROR_IMATRIX_SUB_ARRAY_TOO_BIG;
typedef ERROR__ROW_OR_COL_NOT_IN_MAT<imatrix> ERROR_IMATRIX_ROW_OR_COL_NOT_IN_MAT;
typedef ERROR__WRONG_ROW_OR_COL<imatrix> ERROR_IMATRIX_WRONG_ROW_OR_COL;
typedef ERROR__WRONG_BOUNDARIES<civector> ERROR_CIVECTOR_WRONG_BOUNDARIES;
typedef ERROR__NO_MORE_MEMORY<civector> ERROR_CIVECTOR_NO_MORE_MEMORY;
typedef ERROR__SUB_ARRAY_TOO_BIG<civector> ERROR_CIVECTOR_SUB_ARRAY_TOO_BIG;
typedef ERROR__ELEMENT_NOT_IN_VEC<civector> ERROR_CIVECTOR_ELEMENT_NOT_IN_VEC;
typedef ERROR__TYPE_CAST_OF_THICK_OBJ<civector> ERROR_CIVECTOR_TYPE_CAST_OF_THICK_OBJ;
typedef ERROR__TYPE_CAST_OF_THICK_OBJ<cimatrix> ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ;
typedef ERROR__USE_OF_UNINITIALIZED_OBJ<civector> ERROR_CIVECTOR_USE_OF_UNINITIALIZED_OBJ;
typedef ERROR__USE_OF_UNINITIALIZED_OBJ<cimatrix> ERROR_CIMATRIX_USE_OF_UNINITIALIZED_OBJ;
typedef ERROR__OP_WITH_WRONG_DIM<civector> ERROR_CIVECTOR_OP_WITH_WRONG_DIM;
typedef ERROR__OP_WITH_WRONG_DIM<cimatrix> ERROR_CIMATRIX_OP_WITH_WRONG_DIM;
typedef ERROR__WRONG_BOUNDARIES<cimatrix> ERROR_CIMATRIX_WRONG_BOUNDARIES;
typedef ERROR__NO_MORE_MEMORY<cimatrix> ERROR_CIMATRIX_NO_MORE_MEMORY;
typedef ERROR__SUB_ARRAY_TOO_BIG<cimatrix> ERROR_CIMATRIX_SUB_ARRAY_TOO_BIG;
typedef ERROR__ROW_OR_COL_NOT_IN_MAT<cimatrix> ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT;
typedef ERROR__WRONG_ROW_OR_COL<cimatrix> ERROR_CIMATRIX_WRONG_ROW_OR_COL;
typedef ERROR__WRONG_BOUNDARIES<intvector> ERROR_INTVECTOR_WRONG_BOUNDARIES;
typedef ERROR__NO_MORE_MEMORY<intvector> ERROR_INTVECTOR_NO_MORE_MEMORY;
typedef ERROR__SUB_ARRAY_TOO_BIG<intvector> ERROR_INTVECTOR_SUB_ARRAY_TOO_BIG;
typedef ERROR__ELEMENT_NOT_IN_VEC<intvector> ERROR_INTVECTOR_ELEMENT_NOT_IN_VEC;
typedef ERROR__TYPE_CAST_OF_THICK_OBJ<intvector> ERROR_INTVECTOR_TYPE_CAST_OF_THICK_OBJ;
typedef ERROR__TYPE_CAST_OF_THICK_OBJ<intmatrix> ERROR_INTMATRIX_TYPE_CAST_OF_THICK_OBJ;
typedef ERROR__USE_OF_UNINITIALIZED_OBJ<intvector> ERROR_INTVECTOR_USE_OF_UNINITIALIZED_OBJ;
typedef ERROR__USE_OF_UNINITIALIZED_OBJ<intmatrix> ERROR_INTMATRIX_USE_OF_UNINITIALIZED_OBJ;
typedef ERROR__OP_WITH_WRONG_DIM<intvector> ERROR_INTVECTOR_OP_WITH_WRONG_DIM;
typedef ERROR__OP_WITH_WRONG_DIM<intmatrix> ERROR_INTMATRIX_OP_WITH_WRONG_DIM;
typedef ERROR__WRONG_BOUNDARIES<intmatrix> ERROR_INTMATRIX_WRONG_BOUNDARIES;
typedef ERROR__NO_MORE_MEMORY<intmatrix> ERROR_INTMATRIX_NO_MORE_MEMORY;
typedef ERROR__SUB_ARRAY_TOO_BIG<intmatrix> ERROR_INTMATRIX_SUB_ARRAY_TOO_BIG;
typedef ERROR__ROW_OR_COL_NOT_IN_MAT<intmatrix> ERROR_INTMATRIX_ROW_OR_COL_NOT_IN_MAT;
typedef ERROR__WRONG_ROW_OR_COL<intmatrix> ERROR_INTMATRIX_WRONG_ROW_OR_COL;

typedef ERROR__WRONG_BOUNDARIES<l_rvector> ERROR_LRVECTOR_WRONG_BOUNDARIES;
typedef ERROR__NO_MORE_MEMORY<l_rvector> ERROR_LRVECTOR_NO_MORE_MEMORY;
typedef ERROR__SUB_ARRAY_TOO_BIG<l_rvector> ERROR_LRVECTOR_SUB_ARRAY_TOO_BIG;
typedef ERROR__ELEMENT_NOT_IN_VEC<l_rvector> ERROR_LRVECTOR_ELEMENT_NOT_IN_VEC;
typedef ERROR__TYPE_CAST_OF_THICK_OBJ<l_rvector> ERROR_LRVECTOR_TYPE_CAST_OF_THICK_OBJ;
typedef ERROR__TYPE_CAST_OF_THICK_OBJ<l_rmatrix> ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ;
typedef ERROR__USE_OF_UNINITIALIZED_OBJ<l_rvector> ERROR_LRVECTOR_USE_OF_UNINITIALIZED_OBJ;
typedef ERROR__USE_OF_UNINITIALIZED_OBJ<l_rmatrix> ERROR_LRMATRIX_USE_OF_UNINITIALIZED_OBJ;
typedef ERROR__OP_WITH_WRONG_DIM<l_rvector> ERROR_LRVECTOR_OP_WITH_WRONG_DIM;
typedef ERROR__OP_WITH_WRONG_DIM<l_rmatrix> ERROR_LRMATRIX_OP_WITH_WRONG_DIM;
typedef ERROR__WRONG_BOUNDARIES<l_rmatrix> ERROR_LRMATRIX_WRONG_BOUNDARIES;
typedef ERROR__NO_MORE_MEMORY<l_rmatrix> ERROR_LRMATRIX_NO_MORE_MEMORY;
typedef ERROR__SUB_ARRAY_TOO_BIG<l_rmatrix> ERROR_LRMATRIX_SUB_ARRAY_TOO_BIG;
typedef ERROR__ROW_OR_COL_NOT_IN_MAT<l_rmatrix> ERROR_LRMATRIX_ROW_OR_COL_NOT_IN_MAT;
typedef ERROR__WRONG_ROW_OR_COL<l_rmatrix> ERROR_LRMATRIX_WRONG_ROW_OR_COL;

typedef ERROR__WRONG_BOUNDARIES<l_ivector> ERROR_LIVECTOR_WRONG_BOUNDARIES;
typedef ERROR__NO_MORE_MEMORY<l_ivector> ERROR_LIVECTOR_NO_MORE_MEMORY;
typedef ERROR__SUB_ARRAY_TOO_BIG<l_ivector> ERROR_LIVECTOR_SUB_ARRAY_TOO_BIG;
typedef ERROR__ELEMENT_NOT_IN_VEC<l_ivector> ERROR_LIVECTOR_ELEMENT_NOT_IN_VEC;
typedef ERROR__TYPE_CAST_OF_THICK_OBJ<l_ivector> ERROR_LIVECTOR_TYPE_CAST_OF_THICK_OBJ;
typedef ERROR__TYPE_CAST_OF_THICK_OBJ<l_imatrix> ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ;
typedef ERROR__USE_OF_UNINITIALIZED_OBJ<l_ivector> ERROR_LIVECTOR_USE_OF_UNINITIALIZED_OBJ;
typedef ERROR__USE_OF_UNINITIALIZED_OBJ<l_imatrix> ERROR_LIMATRIX_USE_OF_UNINITIALIZED_OBJ;
typedef ERROR__OP_WITH_WRONG_DIM<l_ivector> ERROR_LIVECTOR_OP_WITH_WRONG_DIM;
typedef ERROR__OP_WITH_WRONG_DIM<l_imatrix> ERROR_LIMATRIX_OP_WITH_WRONG_DIM;
typedef ERROR__WRONG_BOUNDARIES<l_imatrix> ERROR_LIMATRIX_WRONG_BOUNDARIES;
typedef ERROR__NO_MORE_MEMORY<l_imatrix> ERROR_LIMATRIX_NO_MORE_MEMORY;
typedef ERROR__SUB_ARRAY_TOO_BIG<l_imatrix> ERROR_LIMATRIX_SUB_ARRAY_TOO_BIG;
typedef ERROR__ROW_OR_COL_NOT_IN_MAT<l_imatrix> ERROR_LIMATRIX_ROW_OR_COL_NOT_IN_MAT;
typedef ERROR__WRONG_ROW_OR_COL<l_imatrix> ERROR_LIMATRIX_WRONG_ROW_OR_COL;

// --------------------------------------------------------------------------

//member functions of class ERROR_IDOTPRECISION_EMPTY_INTERVAL
int ERROR_IDOTPRECISION_EMPTY_INTERVAL::errnum() const throw() { return 2011; }
string ERROR_IDOTPRECISION_EMPTY_INTERVAL::errtext() const throw() { return fkt+": ERROR_IDOTPRECISION_EMPTY_INTERVAL"; }
ERROR_IDOTPRECISION_EMPTY_INTERVAL::ERROR_IDOTPRECISION_EMPTY_INTERVAL() throw() { fkt="<unknown function>"; }
ERROR_IDOTPRECISION_EMPTY_INTERVAL::ERROR_IDOTPRECISION_EMPTY_INTERVAL(const string &f) throw() { fkt=f; }

//member functions of class ERROR_CIDOTPRECISION_EMPTY_INTERVAL
int ERROR_CIDOTPRECISION_EMPTY_INTERVAL::errnum() const throw() { return 2011; }
string ERROR_CIDOTPRECISION_EMPTY_INTERVAL::errtext() const throw() { return fkt+": ERROR_CIDOTPRECISION_EMPTY_INTERVAL"; }
ERROR_CIDOTPRECISION_EMPTY_INTERVAL::ERROR_CIDOTPRECISION_EMPTY_INTERVAL() throw() { fkt="<unknown function>"; }
ERROR_CIDOTPRECISION_EMPTY_INTERVAL::ERROR_CIDOTPRECISION_EMPTY_INTERVAL(const string &f) throw() { fkt=f; }

//member functions of class ERROR_INTERVAL_EMPTY_INTERVAL
int ERROR_INTERVAL_EMPTY_INTERVAL::errnum() const throw() { return 4011; }
string ERROR_INTERVAL_EMPTY_INTERVAL::errtext() const throw() { return fkt+": ERROR_INTERVAL_EMPTY_INTERVAL"; }
ERROR_INTERVAL_EMPTY_INTERVAL::ERROR_INTERVAL_EMPTY_INTERVAL() throw() { fkt="<unknown function>"; }
ERROR_INTERVAL_EMPTY_INTERVAL::ERROR_INTERVAL_EMPTY_INTERVAL(const string &f) throw() { fkt=f; }

//member functions of class ERROR_CINTERVAL_EMPTY_INTERVAL
int ERROR_CINTERVAL_EMPTY_INTERVAL::errnum() const throw() { return 6011; }
string ERROR_CINTERVAL_EMPTY_INTERVAL::errtext() const throw() { return fkt+": ERROR_CINTERVAL_EMPTY_INTERVAL"; }
ERROR_CINTERVAL_EMPTY_INTERVAL::ERROR_CINTERVAL_EMPTY_INTERVAL() throw() { fkt="<unknown function>"; }
ERROR_CINTERVAL_EMPTY_INTERVAL::ERROR_CINTERVAL_EMPTY_INTERVAL(const string &f) throw() { fkt=f; }

//member functions of class ERROR_INTERVAL_STD_FKT_OUT_OF_DEF
int ERROR_INTERVAL_STD_FKT_OUT_OF_DEF::errnum() const throw() { return 4302; }
string ERROR_INTERVAL_STD_FKT_OUT_OF_DEF::errtext() const throw() { return fkt+": ERROR_INTERVAL_STD_FKT_OUT_OF_DEF"; }
ERROR_INTERVAL_STD_FKT_OUT_OF_DEF::ERROR_INTERVAL_STD_FKT_OUT_OF_DEF() throw() { fkt="<unknown function>"; }
ERROR_INTERVAL_STD_FKT_OUT_OF_DEF::ERROR_INTERVAL_STD_FKT_OUT_OF_DEF(const string &f) throw() { fkt=f; }

//member functions of class ERROR_LREAL_STD_FKT_OUT_OF_DEF
int ERROR_LREAL_STD_FKT_OUT_OF_DEF::errnum() const throw() { return 15302; }
string ERROR_LREAL_STD_FKT_OUT_OF_DEF::errtext() const throw() { return fkt+": ERROR_LREAL_STD_FKT_OUT_OF_DEF"; }
ERROR_LREAL_STD_FKT_OUT_OF_DEF::ERROR_LREAL_STD_FKT_OUT_OF_DEF() throw() { fkt="<unknown function>"; }
ERROR_LREAL_STD_FKT_OUT_OF_DEF::ERROR_LREAL_STD_FKT_OUT_OF_DEF(const string &f) throw() { fkt=f; }

//member functions of class ERROR_LINTERVAL_DIV_BY_ZERO
int ERROR_LINTERVAL_DIV_BY_ZERO::errnum() const throw() { return 16010; }
string ERROR_LINTERVAL_DIV_BY_ZERO::errtext() const throw() { return fkt+": ERROR_LINTERVAL_DIV_BY_ZERO"; }
ERROR_LINTERVAL_DIV_BY_ZERO::ERROR_LINTERVAL_DIV_BY_ZERO() throw() { fkt="<unknown function>"; }
ERROR_LINTERVAL_DIV_BY_ZERO::ERROR_LINTERVAL_DIV_BY_ZERO(const string &f) throw() { fkt=f; }

//member functions of class ERROR_LINTERVAL_EMPTY_INTERVAL
int ERROR_LINTERVAL_EMPTY_INTERVAL::errnum() const throw() { return 16011; }
string ERROR_LINTERVAL_EMPTY_INTERVAL::errtext() const throw() { return fkt+": ERROR_LINTERVAL_EMPTY_INTERVAL"; }
ERROR_LINTERVAL_EMPTY_INTERVAL::ERROR_LINTERVAL_EMPTY_INTERVAL() throw() { fkt="<unknown function>"; }
ERROR_LINTERVAL_EMPTY_INTERVAL::ERROR_LINTERVAL_EMPTY_INTERVAL(const string &f) throw() { fkt=f; }

//member functions of class ERROR_LINTERVAL_IN_EXACT_CH_OR_IS
int ERROR_LINTERVAL_IN_EXACT_CH_OR_IS::errnum() const throw() { return 16013; }
string ERROR_LINTERVAL_IN_EXACT_CH_OR_IS::errtext() const throw() { return fkt+": ERROR_LINTERVAL_IN_EXACT_CH_OR_IS"; }
ERROR_LINTERVAL_IN_EXACT_CH_OR_IS::ERROR_LINTERVAL_IN_EXACT_CH_OR_IS() throw() { fkt="<unknown function>"; }
ERROR_LINTERVAL_IN_EXACT_CH_OR_IS::ERROR_LINTERVAL_IN_EXACT_CH_OR_IS(const string &f) throw() { fkt=f; }

//member functions of class ERROR_LINTERVAL_WRONG_STAGPREC
int ERROR_LINTERVAL_WRONG_STAGPREC::errnum() const throw() { return 16300; }
string ERROR_LINTERVAL_WRONG_STAGPREC::errtext() const throw() { return fkt+": ERROR_LINTERVAL_WRONG_STAGPREC"; }
ERROR_LINTERVAL_WRONG_STAGPREC::ERROR_LINTERVAL_WRONG_STAGPREC() throw() { fkt="<unknown function>"; }
ERROR_LINTERVAL_WRONG_STAGPREC::ERROR_LINTERVAL_WRONG_STAGPREC(const string &f) throw() { fkt=f; }

//member functions of class ERROR_LINTERVAL_ELEMENT_NOT_IN_LONG
int ERROR_LINTERVAL_ELEMENT_NOT_IN_LONG::errnum() const throw() { return 16301; }
string ERROR_LINTERVAL_ELEMENT_NOT_IN_LONG::errtext() const throw() { return fkt+": ERROR_LINTERVAL_ELEMENT_NOT_IN_LONG"; }
ERROR_LINTERVAL_ELEMENT_NOT_IN_LONG::ERROR_LINTERVAL_ELEMENT_NOT_IN_LONG() throw() { fkt="<unknown function>"; }
ERROR_LINTERVAL_ELEMENT_NOT_IN_LONG::ERROR_LINTERVAL_ELEMENT_NOT_IN_LONG(const string &f) throw() { fkt=f; }

//member functions of class ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF
int ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF::errnum() const throw() { return 16302; }
string ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF::errtext() const throw() { return fkt+": ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF"; }
ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF::ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF() throw() { fkt="<unknown function>"; }
ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF::ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF(const string &f) throw() { fkt=f; }

//member functions of class ERROR_LINTERVAL_FAK_OVERFLOW
int ERROR_LINTERVAL_FAK_OVERFLOW::errnum() const throw() { return 16303; }
string ERROR_LINTERVAL_FAK_OVERFLOW::errtext() const throw() { return fkt+": ERROR_LINTERVAL_FAK_OVERFLOW"; }
ERROR_LINTERVAL_FAK_OVERFLOW::ERROR_LINTERVAL_FAK_OVERFLOW() throw() { fkt="<unknown function>"; }
ERROR_LINTERVAL_FAK_OVERFLOW::ERROR_LINTERVAL_FAK_OVERFLOW(const string &f) throw() { fkt=f; }

//member functions of class REAL_NOT_ALLOWED; Blomquist, 26.01.08;
int REAL_NOT_ALLOWED::errnum() const throw() { return 304; }
string REAL_NOT_ALLOWED::errtext() const throw() { return fkt+": REAL_NOT_ALLOWED"; }
REAL_NOT_ALLOWED::REAL_NOT_ALLOWED() throw() { fkt="<unknown function>"; }
REAL_NOT_ALLOWED::REAL_NOT_ALLOWED(const string &f) throw() { fkt=f; }

//member functions of class REAL_INT_OUT_OF_RANGE; Blomquist, 26.01.08;
int REAL_INT_OUT_OF_RANGE::errnum() const throw() { return 305; }
string REAL_INT_OUT_OF_RANGE::errtext() const throw() { return fkt+": REAL_INT_OUT_OF_RANGE"; }
REAL_INT_OUT_OF_RANGE::REAL_INT_OUT_OF_RANGE() throw() { fkt="<unknown function>"; }
REAL_INT_OUT_OF_RANGE::REAL_INT_OUT_OF_RANGE(const string &f) throw() { fkt=f; }

//member functions of class NO_BRACKETS_IN_STRING; Blomquist, 26.01.08;
int NO_BRACKETS_IN_STRING::errnum() const throw() { return 306; }
string NO_BRACKETS_IN_STRING::errtext() const throw() { return fkt+": NO_BRACKETS_IN_STRING"; }
NO_BRACKETS_IN_STRING::NO_BRACKETS_IN_STRING() throw() { fkt="<unknown function>"; }
NO_BRACKETS_IN_STRING::NO_BRACKETS_IN_STRING(const string &f) throw() { fkt=f; }


//explicit instantiation of template types
/*
template class ERROR__WRONG_BOUNDARIES<rvector>; 
template class ERROR__NO_MORE_MEMORY<rvector>;
template class ERROR__SUB_ARRAY_TOO_BIG<rvector>;
template class ERROR__ELEMENT_NOT_IN_VEC<rvector>;
template class ERROR__TYPE_CAST_OF_THICK_OBJ<rvector>;
template class ERROR__TYPE_CAST_OF_THICK_OBJ<rmatrix>;
template class ERROR__USE_OF_UNINITIALIZED_OBJ<rvector>;
template class ERROR__USE_OF_UNINITIALIZED_OBJ<rmatrix>;
template class ERROR__OP_WITH_WRONG_DIM<rvector>;
template class ERROR__OP_WITH_WRONG_DIM<rmatrix>;
template class ERROR__WRONG_BOUNDARIES<rmatrix>;
template class ERROR__NO_MORE_MEMORY<rmatrix>;
template class ERROR__SUB_ARRAY_TOO_BIG<rmatrix>;
template class ERROR__ROW_OR_COL_NOT_IN_MAT<rmatrix>;
template class ERROR__WRONG_ROW_OR_COL<rmatrix>;
template class ERROR__WRONG_BOUNDARIES<cvector>;
template class ERROR__NO_MORE_MEMORY<cvector>;
template class ERROR__SUB_ARRAY_TOO_BIG<cvector>;
template class ERROR__ELEMENT_NOT_IN_VEC<cvector>;
template class ERROR__TYPE_CAST_OF_THICK_OBJ<cvector>;
template class ERROR__TYPE_CAST_OF_THICK_OBJ<cmatrix>;
template class ERROR__USE_OF_UNINITIALIZED_OBJ<cvector>;
template class ERROR__USE_OF_UNINITIALIZED_OBJ<cmatrix>;
template class ERROR__OP_WITH_WRONG_DIM<cvector>;
template class ERROR__OP_WITH_WRONG_DIM<cmatrix>;
template class ERROR__WRONG_BOUNDARIES<cmatrix>;
template class ERROR__NO_MORE_MEMORY<cmatrix>;
template class ERROR__SUB_ARRAY_TOO_BIG<cmatrix>;
template class ERROR__ROW_OR_COL_NOT_IN_MAT<cmatrix>;
template class ERROR__WRONG_ROW_OR_COL<cmatrix>;
template class ERROR__WRONG_BOUNDARIES<ivector>;
template class ERROR__NO_MORE_MEMORY<ivector>;
template class ERROR__SUB_ARRAY_TOO_BIG<ivector>;
template class ERROR__ELEMENT_NOT_IN_VEC<ivector>;
template class ERROR__TYPE_CAST_OF_THICK_OBJ<ivector>;
template class ERROR__TYPE_CAST_OF_THICK_OBJ<imatrix>;
template class ERROR__USE_OF_UNINITIALIZED_OBJ<ivector>;
template class ERROR__USE_OF_UNINITIALIZED_OBJ<imatrix>;
template class ERROR__OP_WITH_WRONG_DIM<ivector>;
template class ERROR__OP_WITH_WRONG_DIM<imatrix>;
template class ERROR__WRONG_BOUNDARIES<imatrix>;
template class ERROR__NO_MORE_MEMORY<imatrix>;
template class ERROR__SUB_ARRAY_TOO_BIG<imatrix>;
template class ERROR__ROW_OR_COL_NOT_IN_MAT<imatrix>;
template class ERROR__WRONG_ROW_OR_COL<imatrix>;
template class ERROR__WRONG_BOUNDARIES<civector>;
template class ERROR__NO_MORE_MEMORY<civector>;
template class ERROR__SUB_ARRAY_TOO_BIG<civector>;
template class ERROR__ELEMENT_NOT_IN_VEC<civector>;
template class ERROR__TYPE_CAST_OF_THICK_OBJ<civector>;
template class ERROR__TYPE_CAST_OF_THICK_OBJ<cimatrix>;
template class ERROR__USE_OF_UNINITIALIZED_OBJ<civector>;
template class ERROR__USE_OF_UNINITIALIZED_OBJ<cimatrix>;
template class ERROR__OP_WITH_WRONG_DIM<civector>;
template class ERROR__OP_WITH_WRONG_DIM<cimatrix>;
template class ERROR__WRONG_BOUNDARIES<cimatrix>;
template class ERROR__NO_MORE_MEMORY<cimatrix>;
template class ERROR__SUB_ARRAY_TOO_BIG<cimatrix>;
template class ERROR__ROW_OR_COL_NOT_IN_MAT<cimatrix>;
template class ERROR__WRONG_ROW_OR_COL<cimatrix>;
template class ERROR__WRONG_BOUNDARIES<intvector>;
template class ERROR__NO_MORE_MEMORY<intvector>;
template class ERROR__SUB_ARRAY_TOO_BIG<intvector>;
template class ERROR__ELEMENT_NOT_IN_VEC<intvector>;
template class ERROR__TYPE_CAST_OF_THICK_OBJ<intvector>;
template class ERROR__TYPE_CAST_OF_THICK_OBJ<intmatrix>;
template class ERROR__USE_OF_UNINITIALIZED_OBJ<intvector>;
template class ERROR__USE_OF_UNINITIALIZED_OBJ<intmatrix>;
template class ERROR__OP_WITH_WRONG_DIM<intvector>;
template class ERROR__OP_WITH_WRONG_DIM<intmatrix>;
template class ERROR__WRONG_BOUNDARIES<intmatrix>;
template class ERROR__NO_MORE_MEMORY<intmatrix>;
template class ERROR__SUB_ARRAY_TOO_BIG<intmatrix>;
template class ERROR__ROW_OR_COL_NOT_IN_MAT<intmatrix>;
template class ERROR__WRONG_ROW_OR_COL<intmatrix>;

template class ERROR__WRONG_BOUNDARIES<l_rvector>;
template class ERROR__NO_MORE_MEMORY<l_rvector>;
template class ERROR__SUB_ARRAY_TOO_BIG<l_rvector>;
template class ERROR__ELEMENT_NOT_IN_VEC<l_rvector>;
template class ERROR__TYPE_CAST_OF_THICK_OBJ<l_rvector>;
template class ERROR__TYPE_CAST_OF_THICK_OBJ<l_rmatrix>;
template class ERROR__USE_OF_UNINITIALIZED_OBJ<l_rvector>;
template class ERROR__USE_OF_UNINITIALIZED_OBJ<l_rmatrix>;
template class ERROR__OP_WITH_WRONG_DIM<l_rvector>;
template class ERROR__OP_WITH_WRONG_DIM<l_rmatrix>;
template class ERROR__WRONG_BOUNDARIES<l_rmatrix>;
template class ERROR__NO_MORE_MEMORY<l_rmatrix>;
template class ERROR__SUB_ARRAY_TOO_BIG<l_rmatrix>;
template class ERROR__ROW_OR_COL_NOT_IN_MAT<l_rmatrix>;
template class ERROR__WRONG_ROW_OR_COL<l_rmatrix>;

template class ERROR__WRONG_BOUNDARIES<l_ivector>;
template class ERROR__NO_MORE_MEMORY<l_ivector>;
template class ERROR__SUB_ARRAY_TOO_BIG<l_ivector>;
template class ERROR__ELEMENT_NOT_IN_VEC<l_ivector>;
template class ERROR__TYPE_CAST_OF_THICK_OBJ<l_ivector>;
template class ERROR__TYPE_CAST_OF_THICK_OBJ<l_imatrix>;
template class ERROR__USE_OF_UNINITIALIZED_OBJ<l_ivector>;
template class ERROR__USE_OF_UNINITIALIZED_OBJ<l_imatrix>;
template class ERROR__OP_WITH_WRONG_DIM<l_ivector>;
template class ERROR__OP_WITH_WRONG_DIM<l_imatrix>;
template class ERROR__WRONG_BOUNDARIES<l_imatrix>;
template class ERROR__NO_MORE_MEMORY<l_imatrix>;
template class ERROR__SUB_ARRAY_TOO_BIG<l_imatrix>;
template class ERROR__ROW_OR_COL_NOT_IN_MAT<l_imatrix>;
template class ERROR__WRONG_ROW_OR_COL<l_imatrix>;
*/
// --------------------------------------------------------------------------

} // namespace cxsc 

