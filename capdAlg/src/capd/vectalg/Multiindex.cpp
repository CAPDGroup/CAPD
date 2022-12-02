
/////////////////////////////////////////////////////////////////////////////
/// @file vectalg/Multiindex.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include "capd/basicalg/factrial.h"
#include "capd/vectalg/Multiindex.h"
#include "capd/vectalg/Container.hpp"

namespace capd{
namespace vectalg{

// static members
std::vector<Multipointer::IndicesSet> Multipointer::indices;
Multipointer::size_type Multipointer::maxKnownLevel=0;

// static functions

// ---------------------------------------------------------------------

const Multipointer::IndicesSet& Multipointer::generateList(size_type p, size_type k)
{
   if(k<0 || k>p || p<1)
      throw std::runtime_error("subIndices: wrong arguments! Call with 1<= second <= first");
   
   if(p>maxKnownLevel)
   {
      indices.resize(((p+1)*p)/2);
      for(size_type i=maxKnownLevel+1;i<=p;++i)
         computeNextLevel();
   }
   return getList(p,k);
}

// ---------------------------------------------------------------------

void Multipointer::computeNextLevel()
{
  if(maxKnownLevel==0)
  {
    Multipointer first(1);
    MultipointersVector mv;
    mv.push_back(first);
    indices[0].push_back(mv);
    maxKnownLevel=1;
  }else
  {
    size_type p = maxKnownLevel+1;
    Multipointer first(p,true);

    int i=0;
    iterator b=first.begin(), e=first.end();
    while(b!=e)
    {
       *b = i;
       ++b;
       ++i;
    }
    //for k=0
    MultipointersVector mv;
    mv.push_back(first);
    getList(p,1).push_back(mv);

    for(size_type k=2;k<=maxKnownLevel;++k)
    {

      Multipointer last(1,true);
      last[0] = maxKnownLevel;
      IndicesSet& current = getList(p,k);
      
      const IndicesSet& lower = getList(maxKnownLevel,k-1);
      IndicesSet::const_iterator b = lower.begin(), e=lower.end();
      while(b!=e)
      {
        MultipointersVector mv = *b;
        mv.push_back(last);
        current.push_back(mv);
        ++b;
      }

      const IndicesSet& higher = getList(maxKnownLevel,k);
      b = higher.begin();
      e = higher.end();
      while(b!=e)
      {
        for(unsigned int s=0;s<(*b).size();++s)
        {
          MultipointersVector copy = *b;
          copy[s] = sumMultipointers((*b)[s],last);
          current.push_back(copy);
        }
        ++b;
      }
    }

    // for k=maxKnownLevel+1
    MultipointersVector last;
    for(size_type s=0;s<p;++s)
    {
      Multipointer mv(1,true);
      mv[0] = s;
      last.push_back(mv);
    }
    getList(p,p).push_back(last);

     ++maxKnownLevel;
  }
}

// ---------------------------------------------------------------------

void Multiindex::generateList(size_type n, size_type k, IndicesSet& result)
{
   result.resize(k);
   //we storage the first level, i.e. k=1
   size_type i;
   for(i=0;i<n;++i){
     Multiindex mi(1,true);
     mi[0] = i;
     result[0].push_back(mi);
   }
   
   for(i=1;i<k;++i)
   {
      MultiindexVector::iterator b = result[i-1].begin(), e=result[i-1].end();
      while(b!=e)
      {
         for(size_type j=0;j<n;++j)
         {
            Multiindex m(i+1);
            for(size_type c=0;c<i;++c)
               m[c] = (*b)[c];
            m[i] = j;
            result[i].push_back(m);
         }
         ++b;
      }
   }
}

// ---------------------------------------------------------------------

bool Multipointer::hasNext(size_type dim){
  if(this->dimension())
  {
    iterator b=begin(), e=end(), e1=end();
    do
    {
      --e;
      if( ++(*e) % dim )
      {
        int p=*e;
        ++e;
        while(e!=e1)
        {
          *e=p;
          ++e;
        }
        return true;
      }
    }while(b!=e);
  }
  return false;
}

// ---------------------------------------------------------------------

long Multipointer::factorial() const{
   const_iterator b=begin(), e=end();
   long result=1;
   while(b!=e)
   {
      const_iterator temp=b;
      int p = *b;
      do{
         ++b;
      }while(b!=e && *b==p);
      size_t n = b-temp;
      if(n>1)
         result *= ::factorial(n);
   }
   return result;
}

// ---------------------------------------------------------------------
/// Returns multipointer containing entries which indices are in mp
///
/// e.g. for a = (1,3,3,6,7)  mp=(1,2,4)
///   a.subMultipointer(mp)  returns (3,3,7)
///
Multipointer Multipointer::subMultipointer(const Multipointer& mp) const
{
   Multipointer result(mp.dimension(),true);
   iterator i=result.begin();
   const_iterator j=begin();
   const_iterator b=mp.begin(), e=mp.end();
   while(b!=e)
   {
      (*i) = *(j+(*b));
      ++i;
      ++b;
   }
   return result;
}

// ---------------------------------------------------------------------
/// returns sum of the multiindex entries
Multiindex::size_type Multiindex::module() const{
   const_iterator b=begin(), e=end();
   int result=0;
   while(b!=e)
   {
      result += (*b); // assume Multiindex has nonnegative coordinates only
      ++b;
   }
   return result;
}

// ---------------------------------------------------------------------
/// for multiindex (a,b,..,n) returns a!b!...n!
long Multiindex::factorial() const{
   const_iterator b=begin(), e=end();
   long result=1;
   while(b!=e)
   {
      if((*b)>1)
         result *= ::factorial(*b);
      ++b;
   }
   return result;
}

// ---------------------------------------------------------------------

Multiindex::Multiindex(size_type dim, const Multipointer& mp) : Vector<int,0>(dim)
{
   Multipointer::const_iterator b=mp.begin(), e=mp.end();
   while(b!=e)
   {
      ++ ((*this))[*b];
      ++b;
   }
}


// ---------------------------------------------------------------------

Multipointer::Multipointer(const Multiindex& mi) : Vector<int,0>(mi.module(),true)
{
   iterator i=begin();
   for(size_type j=0;j<mi.dimension();++j)
   {
      for(int r=0;r<mi[j]; ++r)
      {
         (*i) = j;
         ++i;
      }
   }
}

// ---------------------------------------------------------------------

// the following function computes the next multipointer after mp
// it returns false if 'a' is zero multiindex

bool Multiindex::hasNext()
{
  if(this->dimension()<2) return false;
  if(this->data[0]!=0) {
    this->data[0]--;
    this->data[1]++;
    return true;
  }
  for(size_type i=1;i<this->dimension()-1;++i)
  {
    if(this->data[i]!=0){
      this->data[0] = this->data[i]-1;
      this->data[i]=0;
      this->data[i+1]++;
      return true;
    }
  }
  return false;
}

// ---------------------------------------------------------------------

bool Multiindex::hasNext(int* a, int* b) const {
  for(size_type i=0;i<this->dimension();++i){
    if(b[i]>0){
      b[i]--;
      a[i]++;
      return true;
    }
    b[i] = (*this)[i];
    a[i] = 0;
  }
  return false;
}

// ---------------------------------------------------------------------

bool Multiindex::hasNext(int* a, int* b, size_type j) const{
  for(size_type i=0;i<this->dimension();++i){
    if(b[i]> (i==j)){
      b[i]--;
      a[i]++;
      return true;
    }
    a[i] = (i==j);
    b[i] = (*this)[i]-a[i];
  }
  return false;
}

// ---------------------------------------------------------------------

// this function just concatenate sorted Multipointers to the another sorted Multipointer
Multipointer sumMultipointers(const Multipointer& x, const Multipointer& y)
{
   Multipointer result(x.module()+y.module());
   Multipointer::const_iterator xb=x.begin(), xe=x.end(), yb=y.begin(), ye=y.end();
   Multipointer::iterator b=result.begin();
   while(xb!=xe && yb!=ye)
   {
      if((*xb)<(*yb))
      {
         (*b) = (*xb);
         ++xb;
      }else{
         (*b) = (*yb);
         ++yb;
      }
      ++b;
   }
   if(xb==xe)
   {
      xb = yb;
      xe = ye;
   }
   while(xb!=xe)
   {
      (*b) = (*xb);
      ++b;
      ++xb;
   }
   return result;
}

Multipointer addIndex(const Multipointer & mp, int index) {
  Multipointer result(mp.dimension() + 1);
  Multipointer::iterator res = result.begin();
  Multipointer::const_iterator src = mp.begin(), end = mp.end();
  while(src != end){
    if((index>=0) && (index < *src)) {
      *res = index;
      index = -1;
      res++;
    }
    *res = *src;
    ++res; ++src;
  }
  if(res != result.end()){
    *res = index;
  }
  return result;
}

Multipointer push_back(const Multipointer & mp, int index) {
  Multipointer result(mp.dimension() + 1);
  Multipointer::iterator res = result.begin();
  Multipointer::const_iterator src = mp.begin(), end = mp.end();
  while(src != end){
    *res = *src;
    ++res; ++src;
  }
  *res = index;
  return result;
}

// ----------------------------------------------------------

inline int computeNewton(int d,int l)
{
  return binomial(d+l-1,l);
}

// ----------------------------------------------------------

// The following procedure computes an index of of an element in array that corresponds to the multipointer.
Multipointer::size_type Multipointer::index(size_type dim, size_type maxDegree) const
{
  size_type level = this->module(); // in fact dimension
  if (level<=0) return 0;
  if(level>maxDegree){
    throw std::range_error("Multipointer::index(int dim, int maxDegree): requested degree is to large");
  }
  size_type result=0,i;
  size_type prev = 0;

  for(i=0;i<level;++i)
  {
    if((*this)[i]-prev){
      result += (computeNewton(dim-prev,level-i) - computeNewton(dim-(*this)[i],level-i));
      prev = (*this)[i];
    }
  }
  return result;
}

// ----------------------------------------------------------

Multipointer::size_type Multipointer::index(size_type dim, size_type maxDegree, const Multipointer& sub) const
{
  size_type level = sub.module();
  if (level<=0) return 0;
  if(level>maxDegree)
    throw std::range_error("Multipointer::index(size_type dim, size_type maxDegree,Multipointer): requested degree is to large");

  size_type result=0,i;
  size_type prev = 0;

  for(i=0;i<level;++i)
  {
    size_type s = (*this)[sub[i]];
    if(s-prev){
      result += (computeNewton(dim-prev,level-i) - computeNewton(dim-s,level-i));
      prev = s;
    }
  }
  return result;
}

// ----------------------------------------------------------

Multiindex::size_type Multiindex::index(size_type maxDegree) const
{
  size_type level = this->module(); // sum norm
  if (level<=0) return 0;
  if(level>maxDegree)
    throw std::range_error("Multiindex::index(size_type maxDegree): requested degree is to large");

  size_type result=0,i, prev=0;
  for(i=0;i<this->dimension();++i)
  {
    if((*this)[i]!=0){
      result+= (computeNewton(this->dimension()-prev,level) - computeNewton(this->dimension()-i,level));
      prev = i;
      level -= (*this)[i];
    }
  }
  return result;
}

// has is not really important it does a bit less than indexCount, but the only time
// hasIndex will perform less operations than indexCount will be only when the indexes
// will actually be found and than we will still need to compute indexCount
bool hasIndex(const  Multipointer & mp, int index) {
  for (Multipointer::const_iterator it = mp.begin() ; it != mp.end() && *it <= index; ++it) {
    if (*it == index) return true;
  }
  return false;
}


int indexCount(const Multipointer &mp, int index) {

  int count = 0;
  for (Multipointer::const_iterator it = mp.begin() ; it != mp.end() && *it <= index; ++it) {
    count += (*it == index);
  }
  return count;
}

Multipointer removeIndex(const Multipointer & mp, int index) {
  Multipointer result(mp.dimension() -1 );
  Multipointer::iterator res = result.begin();
  Multipointer::const_iterator src = mp.begin(), end = mp.end();
  bool removed = false;
  while(src != end){
    if (!removed && index == *src) {
    removed = true;
    ++src;
  }
  else{
    *res = *src;
    ++res; ++src;
  }
  }
  return result;
}

}} // namespace capd::vectalg

