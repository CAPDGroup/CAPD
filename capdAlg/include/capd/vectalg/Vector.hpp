/////////////////////////////////////////////////////////////////////////////
/// @file Vector.hpp
///
/// @author Marian Mrozek, Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_VECTALG_VECTOR_HPP_
#define _CAPD_VECTALG_VECTOR_HPP_

#include <cmath>
#include <stack>
#include <stdexcept>
#include <cstdio>
#include <sstream>
#include <algorithm>
#include "capd/basicalg/minmax.h"
#include "capd/basicalg/power.h"
#include "capd/vectalg/Container.hpp"
#include "capd/vectalg/Vector.h"
#include "capd/vectalg/algebraicOperations.hpp"

namespace capd{
namespace vectalg{

//---------------------------constructors---------------------------//

template<typename Scalar,__size_type dim>
Vector<Scalar,dim>::Vector(size_type A_dimension, const ScalarType data[]) : ContainerType(A_dimension,true)
{
  std::copy(data,data+A_dimension,begin());
}

template<typename Scalar,__size_type dim>
Vector<Scalar,dim>::Vector(const char data[]) : ContainerType(dim,true)
{
   std::istringstream str(data);
   str >> *this;
}

template<typename Scalar,__size_type dim>
Vector<Scalar,dim>::Vector(const std::string & data) : ContainerType(dim,true)
{
   std::istringstream str(data);
   str >> *this;
}

template<typename Scalar,__size_type dim>
inline void setDimensionInArray(Vector<Scalar,dim>* , __size_type, __size_type){
}

template<typename Scalar>
inline void setDimensionInArray(Vector<Scalar,0>* t, __size_type N, __size_type _dim){
  for(__size_type i=0;i<N;++i) t[i].resize(_dim);
}

template<typename Scalar,__size_type dim>
Vector<Scalar,dim>* Vector<Scalar,dim>::makeArray(size_type N, size_type _dim)
{
  Vector* result = new Vector[N];
  setDimensionInArray(result,N,_dim);
  return result;
}

template<typename Scalar,__size_type dim>
void Vector<Scalar,dim>::sorting_permutation(typename rebind<int>::other& perm)
{
  typedef typename rebind<int>::other IntVectorType;
  difference_type i=0,j,k;

  if(dimension()!= perm.dimension())
     throw std::range_error("sorting_permutation: Incompatible vector dimensions");
  typename IntVectorType::iterator b=perm.begin(), e=perm.end();
  while(b!=e)
  {
    *b=i;
    ++i;
    ++b;
  }

  difference_type d = dimension();

  for(i=0;i<d;i++)
    for(j=d-1;j>i;j--)
    {
      if((*this)[perm[j]] > (*this)[perm[j-1]])
      {
        k=perm[j-1];
        perm[j-1]=perm[j];
        perm[j]=k;
      }
    }
}

template<typename Scalar,__size_type dim>
template<__size_type dataDim>
Vector<Scalar,dim>::Vector(const Scalar (&data)[dataDim]): ContainerType(dataDim,true)
{
  std::copy(data, data + dataDim, begin());
}

template<typename Scalar,__size_type dim>
template<typename Iterator>
Vector<Scalar,dim>::Vector(Iterator begin, Iterator end): ContainerType((end - begin), true)
{
  std::copy(begin, end, this->begin());
}

//----------------- input-output ----------------------------------//

template<typename Scalar,__size_type dim>
std::ostream& operator<<(std::ostream& out, const Vector<Scalar,dim>& v)
{
  typedef typename Vector<Scalar,dim>::size_type size_type;
  const size_type d = v.dimension();
  out << "{";
  if(d>0){
     //if(v[0]>=Scalar(0)) out << " "; /***** DW it does not work for complex vectors ***/
     out << v[0];
   }
   for(size_type i=1;i<d;i++)
   {
      out << ",";
      // if(v[i]>=Scalar(0)) out << " "; /***** DW it does not work for complex vectors ***/
      out << v[i];
   }
   out << "}";
   return out;
}

template<typename Vector>
std::string vectorToString( const Vector & v, int firstIndex /*= 0*/, int lastIndex /*= -1*/,  int precision /* = -1*/){
  std::ostringstream out;
  if(precision>0)
       out.precision(precision);
  print(out, v, firstIndex, lastIndex);
  return out.str();
}

template<typename Vector>
std::ostream & printVector(std::ostream & out, const Vector & v, int firstIndex /*= 0*/, int lastIndex /*= -1*/){

  const int d = v.dimension();

  if((lastIndex < 0) || (lastIndex >= d))
    lastIndex = d-1;

  if(firstIndex < d) {
    if(firstIndex < 0)
      firstIndex = 0;
    out << "{" << v[firstIndex];
    for(int i=firstIndex+1;i<=lastIndex;i++) {
      out << "," << v[i];
    }
    out << "}";
  } else {
    out << "{}";
  }
  return out;
}

template<typename Scalar,__size_type dim>
std::istream& operator>> (std::istream& inp, Vector<Scalar,dim>& v)
{
   std::deque<Scalar> st;
   Scalar s;
   int ch;

   while('{'!=(ch=inp.get()) && ch!=EOF)
     ;
   if(ch!= EOF)
   {
/*
      // -- begin of added lines for empty vectors
      while(' '==(ch=inp.peek())) ch=inp.get();
      if('}'==(ch=inp.peek())){
        ch=inp.get();
        return inp;
      }
      // -- end of added lines for empty vectors
*/
      inp >> s;
      st.push_back(s);
      do{
         do{
            ch=inp.get();
         }while(isspace(ch));
         if(ch==','){
            inp >> s;
            st.push_back(s);
         }
      }while(ch!='}' && ch!=EOF);
   }
   if(inp.eof())
       throw std::ios_base::failure("EOF encountered when reading a vector");
   v.resize(st.size());
   std::copy(st.begin(), st.end(), v.begin());
   return inp;
}

template<typename Scalar, __size_type dim>
std::string cppReprezentation(const Vector<Scalar,dim> & A, const std::string& varName,
			      const std::string& typeName)
{
  std::stringstream out;
  out << "capd::vectalg::Vector<" << typeName << ", " << dim << "> " << varName << "(";

  if (A.dimension() > 0) {
    out << "(" << typeName << "[" << A.dimension() << "])" << A;
  } else {
    out << "(capd::vectalg::__size_type)" << A.dimension();
  }

  out << ");";

  std::string str = out.str();
  std::replace(str.begin(), str.end(), '\n', ' ');

  return str;
}

}} // namespace capd::vectalg

#endif // _CAPD_VECTALG_VECTOR_HPP_
