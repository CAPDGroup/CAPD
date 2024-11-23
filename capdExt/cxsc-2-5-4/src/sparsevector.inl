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

/* CVS $Id: sparsevector.inl,v 1.13 2014/01/30 17:23:49 cxsc Exp $ */

#ifndef _CXSC_SPARSEVECTOR_INL_INCLUDED
#define _CXSC_SPARSEVECTOR_INL_INCLUDED

namespace cxsc {

template<class T>
inline bool is_interval() {
  return false;
}

template<> inline bool is_interval<interval>() {
  return true;
}

template<> inline bool is_interval<cinterval>() {
  return true;
}

template<class Tx, class Ty, class Tres, class TDot>
inline Tres spsp_vv_mult(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator*(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  if(opdotprec == 1 && !is_interval<Tres>()) {

    Tres dot(0.0);

    for(int i=0 ; i<v1.get_nnz() ; i++) {
      for(int j=0 ; j<v2.get_nnz() && v2.p[j]<=v1.p[i] ; j++) {
         if(v1.p[i] == v2.p[j])
           dot += v1.x[i] * v2.x[j];
      }
    }

    return dot;

  } else {

    TDot dot(opdotprec);

    for(int i=0 ; i<v1.get_nnz() ; i++) {
      for(int j=0 ; j<v2.get_nnz() && v2.p[j]<=v1.p[i] ; j++) {
         if(v1.p[i] == v2.p[j])
           dot.add_dot(v1.x[i], v2.x[j]);
      }
    }

    return dot.result();
  }

}

template<class Tx, class Ty, class Tres, class TDot>
inline Tres slsp_vv_mult(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator*(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  
  if(opdotprec == 1 && !is_interval<Tres>()) {
    Tres dot(0.0);

    for(int i=v1.start ; i<=v1.end ; i++) {
      for(int j=0 ; j<v2.get_nnz() && v2.p[j]<=v1.p[i]-v1.offset ; j++) {
         if(v1.p[i]-v1.offset == v2.p[j])
           dot += v1.x[i] * v2.x[j];
      }
    }

    return dot;

  } else {

    TDot dot(opdotprec);

    for(int i=v1.start ; i<=v1.end ; i++) {
      for(int j=0 ; j<v2.get_nnz() && v2.p[j]<=v1.p[i]-v1.offset ; j++) {
         if(v1.p[i]-v1.offset == v2.p[j])
           dot.add_dot(v1.x[i], v2.x[j]);
      }
    }

    return dot.result();
  }

}


template<class Tx, class Ty, class Tres, class TDot>
inline Tres spsl_vv_mult(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator*(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif

  if(opdotprec == 1 && !is_interval<Tres>()) {
    Tres dot(0.0);

    for(int i=0 ; i<v1.get_nnz() ; i++) {
      for(int j=v2.start ; j<v2.end+1 && v2.p[j]-v2.offset<=v1.p[i] ; j++) {
         if(v1.p[i] == v2.p[j]-v2.offset)
           dot += v1.x[i] * v2.x[j];
      }
    }

    return dot;

  } else {

    TDot dot(opdotprec);

    for(int i=0 ; i<v1.get_nnz() ; i++) {
      for(int j=v2.start ; j<v2.end+1 && v2.p[j]-v2.offset<=v1.p[i] ; j++) {
         if(v1.p[i] == v2.p[j]-v2.offset)
           dot.add_dot(v1.x[i], v2.x[j]);
      }
    }

    return dot.result();
  }
}

template<class Tx, class Ty, class Tres, class TDot>
inline Tres spf_vv_mult(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=VecLen(v2)) cxscthrow(OP_WITH_WRONG_DIM("operator*(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  
  if(opdotprec == 1 && !is_interval<Tres>()) {
    Tres dot(0.0);
    int lb = Lb(v2);

    for(int i=0 ; i<v1.get_nnz() ; i++)
      dot += v1.x[i] * v2[v1.p[i]+lb];

    return dot;

  } else {

    TDot dot(opdotprec);
    int lb = Lb(v2);

    for(int i=0 ; i<v1.get_nnz() ; i++)
      dot.add_dot(v1.x[i], v2[v1.p[i]+lb]);

    return dot.result();
  }
}

template<class Tx, class Ty, class Tres, class TDot>
inline Tres fsp_vv_mult(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(VecLen(v1)!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator*(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  
  if(opdotprec == 1 && !is_interval<Tres>()) {

    Tres dot(0.0);
    int lb = Lb(v1);

    for(int i=0 ; i<v2.get_nnz() ; i++)
      dot += v1[v2.p[i]+lb] * v2.x[i];

    return dot;

  } else {

    TDot dot(opdotprec);
    int lb = Lb(v1);

    for(int i=0 ; i<v2.get_nnz() ; i++)
      dot.add_dot(v1[v2.p[i]+lb], v2.x[i]);

    return dot.result();
  }

}

template<class Tx, class Ty, class Tres, class TDot>
inline Tres slf_vv_mult(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=VecLen(v2)) cxscthrow(OP_WITH_WRONG_DIM("operator*(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  if(opdotprec == 1 && !is_interval<Tres>()) {
    Tres dot(0.0);
    int lb = Lb(v2);

    for(int i=v1.start ; i<=v1.end ; i++)
      dot += v1.x[i] * v2[v1.p[i]-v1.offset+lb];

    return dot;

  } else {

    TDot dot(opdotprec);
    int lb = Lb(v2);

    for(int i=v1.start ; i<=v1.end ; i++)
      dot.add_dot(v1.x[i], v2[v1.p[i]-v1.offset+lb]);

    return dot.result();
  }
}

template<class Tx, class Ty, class Tres, class TDot>
inline Tres fsl_vv_mult(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(VecLen(v1)!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator*(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  if(opdotprec == 1 && !is_interval<Tres>()) {
    Tres dot(0.0);
    int lb = Lb(v1);

    for(int i=v2.start ; i<=v2.end ; i++)
      dot += v1[v2.p[i]-v2.offset+lb] * v2.x[i];

    return dot;

  } else {

    TDot dot(opdotprec);
    int lb = Lb(v1);

    for(int i=v2.start ; i<=v2.end ; i++)
      dot.add_dot(v1[v2.p[i]-v2.offset+lb], v2.x[i]);

    return dot.result();
  }
}

template<class Tx, class Ty, class Tres, class TDot>
inline Tres slsl_vv_mult(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator*(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  if(opdotprec == 1 && !is_interval<Tres>()) {

    Tres dot(0.0);

    for(int i=v1.start ; i<=v1.end ; i++) {
      for(int j=v2.start ; j<v2.end+1 && v2.p[j]-v2.offset<=v1.p[i]-v1.offset ; j++) {
         if(v1.p[i]-v1.offset == v2.p[j]-v2.offset)
           dot += v1.x[i] * v2.x[j];
      }
    }

    return dot;

  } else {

    TDot dot(opdotprec);

    for(int i=v1.start ; i<=v1.end ; i++) {
      for(int j=v2.start ; j<v2.end+1 && v2.p[j]-v2.offset<=v1.p[i]-v1.offset ; j++) {
         if(v1.p[i]-v1.offset == v2.p[j]-v2.offset)
           dot.add_dot(v1.x[i], v2.x[j]);
      }
    }

    return dot.result();
  }

}

//-----------------------------------------------------------------------------------------

template<class Tv, class Ts, class Tres>
inline Tres sp_vs_div(const Tv& v, const Ts& s) {
  Tres res(v);

  for(int i=0 ; i<res.get_nnz() ; i++)
    res.x[i] /= s;

  return res;
}

template<class Tv, class Ts, class Tres>
inline Tres sp_vs_mult(const Tv& v, const Ts& s) {
  Tres res(v);

  for(int i=0 ; i<res.get_nnz() ; i++)
    res.x[i] *= s;

  return res;
}

template<class Ts, class Tv, class Tres>
inline Tres sp_sv_mult(const Ts& s, const Tv& v) {
  Tres res(v);

  for(int i=0 ; i<res.get_nnz() ; i++)
    res.x[i] *= s;

  return res;
}

//----------------------------------------------------------------------------------------

template<class Tx, class Ty, class Tres, class TType>
inline Tres spsp_vv_add(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator+(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v1.n, v1.get_nnz()+v2.get_nnz());

  int i=0, j=0;

  while(i<v1.get_nnz() && j<v2.get_nnz()) {
    if(v1.p[i] < v2.p[j]) {
      res.p.push_back(v1.p[i]);
      res.x.push_back(TType(v1.x[i]));
      i++;
    } else if (v1.p[i] == v2.p[j]) {
      res.p.push_back(v1.p[i]);
      res.x.push_back(TType(v1.x[i]+v2.x[j]));
      i++; j++;
    } else if (v1.p[i] > v2.p[j]) {
      res.p.push_back(v2.p[j]);
      res.x.push_back(TType(v2.x[j]));
      j++;
    }
  }

  for( ; i<v1.get_nnz() ; i++) {
      res.p.push_back(v1.p[i]);
      res.x.push_back(TType(v1.x[i]));
  }      

  for( ; j<v2.get_nnz() ; j++) {
      res.p.push_back(v2.p[j]);
      res.x.push_back(TType(v2.x[j]));
  }      

  return res;
}


template<class Tx, class Ty, class Tres, class TType>
inline Tres slsp_vv_add(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator+(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v1.n, v1.get_nnz()+v2.get_nnz());
  
  int i=v1.start, j=0;

  while(i<=v1.end && j<v2.get_nnz()) {

    if(v1.p[i]-v1.offset < v2.p[j]) {

      res.p.push_back(v1.p[i]-v1.offset);
      res.x.push_back(TType(v1.x[i]));
      i++;

    } else if(v1.p[i]-v1.offset == v2.p[j]){

      res.p.push_back(v2.p[j]);
      res.x.push_back(TType(v1.x[i]+v2.x[j]));
      i++; j++;

    } else if(v1.p[i]-v1.offset > v2.p[j]) {

      res.p.push_back(v2.p[j]);
      res.x.push_back(TType(v2.x[j]));
      j++;
    }
  }

  for( ; i<=v1.end ; i++) {
    res.p.push_back(v1.p[i]-v1.offset);
    res.x.push_back(TType(v1.x[i]));
  }      

  for( ; j<v2.get_nnz() ; j++) {
    res.p.push_back(v2.p[j]);
    res.x.push_back(TType(v2.x[j]));
  }      

  return res;
}

template<class Tx, class Ty, class Tres, class TType>
inline Tres spsl_vv_add(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator+(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v1.n, v1.get_nnz()+v2.get_nnz());
  
  int i=0, j=v2.start;

  while(i<v1.get_nnz() && j<=v2.end) {

    if(v1.p[i] < v2.p[j]-v2.offset) {

      res.p.push_back(v1.p[i]);
      res.x.push_back(TType(v1.x[i]));
      i++;

    } else if(v1.p[i] == v2.p[j]-v2.offset) {

      res.p.push_back(v1.p[i]);
      res.x.push_back(TType(v1.x[i]+v2.x[j]));
      i++; j++;

    } else if(v1.p[i] > v2.p[j]-v2.offset) {

      res.p.push_back(v2.p[j]-v2.offset);
      res.x.push_back(TType(v2.x[j]));
      j++;

    }

  }

  for( ; i<v1.get_nnz() ; i++) {
    res.p.push_back(v1.p[i]);
    res.x.push_back(TType(v1.x[i]));
  }      

  for( ; j<=v2.end ; j++) {
    res.p.push_back(v2.p[j]-v2.offset);
    res.x.push_back(TType(v2.x[j]));
  }      

  return res;
}

template<class Tx, class Ty, class Tres>
inline Tres spf_vv_add(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=VecLen(v2)) cxscthrow(OP_WITH_WRONG_DIM("operator+(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v2);
  int lb = Lb(res);

  for(int i=0 ; i<v1.get_nnz() ; i++) 
    res[v1.p[i]+lb] += v1.x[i];
  
  return res;
}

template<class Tx, class Ty, class Tres>
inline Tres fsp_vv_add(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(VecLen(v1)!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator+(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v1);
  int lb = Lb(res);

  for(int i=0 ; i<v2.get_nnz() ; i++) 
    res[v2.p[i]+lb] += v2.x[i];
  
  return res;
}

template<class Tx, class Ty, class Tres>
inline Tres slf_vv_add(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=VecLen(v2)) cxscthrow(OP_WITH_WRONG_DIM("operator+(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v2);
  int lb = Lb(res);

  for(int i=v1.start ; i<=v1.end ; i++) 
    res[v1.p[i]-v1.offset+lb] += v1.x[i];
  
  return res;
}

template<class Tx, class Ty, class Tres>
inline Tres fsl_vv_add(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(VecLen(v1)!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator+(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v1);
  int lb = Lb(res);

  for(int i=v2.start ; i<=v2.end ; i++) 
    res[v2.p[i]-v2.offset+lb] += v2.x[i];
  
  return res;

}

template<class Tx, class Ty, class Tres, class TType>
inline Tres slsl_vv_add(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator+(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v1.n, v1.get_nnz()+v2.get_nnz());
  
  int i=v1.start, j=v2.start;

  while(i<=v1.end && j<=v2.end) {

    if(v1.p[i]-v1.offset < v2.p[j]-v2.offset) {

      res.p.push_back(v1.p[i]-v1.offset);
      res.x.push_back(TType(v1.x[i]));
      i++;

    } else if (v1.p[i]-v1.offset == v2.p[j]-v2.offset) {

      res.p.push_back(v2.p[j]-v2.offset);
      res.x.push_back(TType(v1.x[i]+v2.x[j]));
      i++; j++;

    } else if (v1.p[i]-v1.offset > v2.p[j]-v2.offset) {
    
      res.p.push_back(v2.p[j]-v2.offset);
      res.x.push_back(TType(v2.x[j]));
      j++;

    }

  }

  for( ; i<=v1.end ; i++) {
    res.p.push_back(v1.p[i]-v1.offset);
    res.x.push_back(TType(v1.x[i]));
  }      

  for( ; j<=v2.end ; j++) {
    res.p.push_back(v2.p[j]-v2.offset);
    res.x.push_back(TType(v2.x[j]));
  }      

  return res;
}


//----------------------------------------------------------------------------------------

template<class Tx, class Ty, class Tres, class TType>
inline Tres spsp_vv_sub(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator-(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v1.n, v1.get_nnz()+v2.get_nnz());

  int i=0, j=0;

  while(i<v1.get_nnz() && j<v2.get_nnz()) {
    if(v1.p[i] < v2.p[j]) {
      res.p.push_back(v1.p[i]);
      res.x.push_back(TType(v1.x[i]));
      i++;
    } else if (v1.p[i] == v2.p[j]) {
      res.p.push_back(v1.p[i]);
      res.x.push_back(TType(v1.x[i]-v2.x[j]));
      i++; j++;
    } else if (v1.p[i] > v2.p[j]) {
      res.p.push_back(v2.p[j]);
      res.x.push_back(TType(-v2.x[j]));
      j++;
    }
  }

  for( ; i<v1.get_nnz() ; i++) {
      res.p.push_back(v1.p[i]);
      res.x.push_back(TType(v1.x[i]));
  }      

  for( ; j<v2.get_nnz() ; j++) {
      res.p.push_back(v2.p[j]);
      res.x.push_back(TType(-v2.x[j]));
  }      

  return res;
}


template<class Tx, class Ty, class Tres, class TType>
inline Tres slsp_vv_sub(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator-(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v1.n, v1.get_nnz()+v2.get_nnz());
  
  int i=v1.start, j=0;

  while(i<=v1.end && j<v2.get_nnz()) {

    if(v1.p[i]-v1.offset < v2.p[j]) {

      res.p.push_back(v1.p[i]-v1.offset);
      res.x.push_back(TType(v1.x[i]));
      i++;

    } else if(v1.p[i]-v1.offset == v2.p[j]){

      res.p.push_back(v2.p[j]);
      res.x.push_back(TType(v1.x[i]-v2.x[j]));
      i++; j++;

    } else if(v1.p[i]-v1.offset > v2.p[j]) {

      res.p.push_back(v2.p[j]);
      res.x.push_back(TType(-v2.x[j]));
      j++;
    }
  }

  for( ; i<=v1.end ; i++) {
    res.p.push_back(v1.p[i]-v1.offset);
    res.x.push_back(TType(v1.x[i]));
  }      

  for( ; j<v2.get_nnz() ; j++) {
    res.p.push_back(v2.p[j]);
    res.x.push_back(TType(-v2.x[j]));
  }      

  return res;
}

template<class Tx, class Ty, class Tres, class TType>
inline Tres spsl_vv_sub(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator-(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v1.n, v1.get_nnz()+v2.get_nnz());
  
  int i=0, j=v2.start;

  while(i<v1.get_nnz() && j<=v2.end) {

    if(v1.p[i] < v2.p[j]-v2.offset) {

      res.p.push_back(v1.p[i]);
      res.x.push_back(TType(v1.x[i]));
      i++;

    } else if(v1.p[i] == v2.p[j]-v2.offset) {

      res.p.push_back(v1.p[i]);
      res.x.push_back(TType(v1.x[i]-v2.x[j]));
      i++; j++;

    } else if(v1.p[i] > v2.p[j]-v2.offset) {

      res.p.push_back(v2.p[j]-v2.offset);
      res.x.push_back(TType(-v2.x[j]));
      j++;

    }

  }

  for( ; i<v1.get_nnz() ; i++) {
    res.p.push_back(v1.p[i]);
    res.x.push_back(TType(v1.x[i]));
  }      

  for( ; j<=v2.end ; j++) {
    res.p.push_back(v2.p[j]-v2.offset);
    res.x.push_back(TType(-v2.x[j]));
  }      

  return res;
}

template<class Tx, class Ty, class Tres>
inline Tres spf_vv_sub(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=VecLen(v2)) cxscthrow(OP_WITH_WRONG_DIM("operator-(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(-v2);
  int lb = Lb(res);

  for(int i=0 ; i<v1.get_nnz() ; i++) 
    res[v1.p[i]+lb] += v1.x[i];
  
  return res;
}

template<class Tx, class Ty, class Tres>
inline Tres fsp_vv_sub(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(VecLen(v1)!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator-(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v1);
  int lb = Lb(res);

  for(int i=0 ; i<v2.get_nnz() ; i++) 
    res[v2.p[i]+lb] -= v2.x[i];
  
  return res;
}

template<class Tx, class Ty, class Tres>
inline Tres slf_vv_sub(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=VecLen(v2)) cxscthrow(OP_WITH_WRONG_DIM("operator-(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(-v2);
  int lb = Lb(res);

  for(int i=v1.start ; i<=v1.end ; i++) 
    res[v1.p[i]-v1.offset+lb] += v1.x[i];
  
  return res;
}

template<class Tx, class Ty, class Tres>
inline Tres fsl_vv_sub(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(VecLen(v1)!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator-(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v1);
  int lb = Lb(res);

  for(int i=v2.start ; i<=v2.end ; i++) 
    res[v2.p[i]-v2.offset+lb] -= v2.x[i];
  
  return res;

}

template<class Tx, class Ty, class Tres, class TType>
inline Tres slsl_vv_sub(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator-(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v1.n, v1.get_nnz()+v2.get_nnz());
  
  int i=v1.start, j=v2.start;

  while(i<=v1.end && j<=v2.end) {

    if(v1.p[i]-v1.offset < v2.p[j]-v2.offset) {

      res.p.push_back(v1.p[i]-v1.offset);
      res.x.push_back(TType(v1.x[i]));
      i++;

    } else if (v1.p[i]-v1.offset == v2.p[j]-v2.offset) {

      res.p.push_back(v2.p[j]-v2.offset);
      res.x.push_back(TType(v1.x[i]-v2.x[j]));
      i++; j++;

    } else if (v1.p[i]-v1.offset > v2.p[j]-v2.offset) {
    
      res.p.push_back(v2.p[j]-v2.offset);
      res.x.push_back(TType(-v2.x[j]));
      j++;

    }

  }

  for( ; i<=v1.end ; i++) {
    res.p.push_back(v1.p[i]-v1.offset);
    res.x.push_back(TType(v1.x[i]));
  }      

  for( ; j<=v2.end ; j++) {
    res.p.push_back(v2.p[j]-v2.offset);
    res.x.push_back(TType(-v2.x[j]));
  }      

  return res;
}



//----------------------------------------------------------------------------------------

template<class Tx, class Ty, class Tres, class TType>
inline Tres spsp_vv_hull(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator|(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v1.n, v1.get_nnz()+v2.get_nnz());

  int i=0, j=0;

  while(i<v1.get_nnz() && j<v2.get_nnz()) {
    if(v1.p[i] < v2.p[j]) {
      res.p.push_back(v1.p[i]);
      res.x.push_back(TType(v1.x[i]) | TType(0.0));
      i++;
    } else if (v1.p[i] == v2.p[j]) {
      res.p.push_back(v1.p[i]);
      res.x.push_back(TType(v1.x[i] | v2.x[j]));
      i++; j++;
    } else if (v1.p[i] > v2.p[j]) {
      res.p.push_back(v2.p[j]);
      res.x.push_back(TType(0.0) | TType(v2.x[j]));
      j++;
    }
  }

  for( ; i<v1.get_nnz() ; i++) {
      res.p.push_back(v1.p[i]);
      res.x.push_back(TType(v1.x[i]) | TType(0.0));
  }      

  for( ; j<v2.get_nnz() ; j++) {
      res.p.push_back(v2.p[j]);
      res.x.push_back(TType(0.0) | TType(v2.x[j]));
  }      

  return res;
}


template<class Tx, class Ty, class Tres, class TType>
inline Tres slsp_vv_hull(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator|(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v1.n, v1.get_nnz()+v2.get_nnz());
  
  int i=v1.start, j=0;

  while(i<=v1.end && j<v2.get_nnz()) {

    if(v1.p[i]-v1.offset < v2.p[j]) {

      res.p.push_back(v1.p[i]-v1.offset);
      res.x.push_back(TType(v1.x[i]) | TType(0.0));
      i++;

    } else if(v1.p[i]-v1.offset == v2.p[j]){

      res.p.push_back(v2.p[j]);
      res.x.push_back(TType(v1.x[i] | v2.x[j]));
      i++; j++;

    } else if(v1.p[i]-v1.offset > v2.p[j]) {

      res.p.push_back(v2.p[j]);
      res.x.push_back(TType(0.0) | TType(v2.x[j]));
      j++;
    }
  }

  for( ; i<=v1.end ; i++) {
    res.p.push_back(v1.p[i]-v1.offset);
    res.x.push_back(TType(v1.x[i]) | TType(0.0));
  }      

  for( ; j<v2.get_nnz() ; j++) {
    res.p.push_back(v2.p[j]);
    res.x.push_back(TType(0.0) | TType(v2.x[j]));
  }      

  return res;
}

template<class Tx, class Ty, class Tres, class TType>
inline Tres spsl_vv_hull(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator|(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v1.n, v1.get_nnz()+v2.get_nnz());
  
  int i=0, j=v2.start;

  while(i<v1.get_nnz() && j<=v2.end) {

    if(v1.p[i] < v2.p[j]-v2.offset) {

      res.p.push_back(v1.p[i]);
      res.x.push_back(TType(v1.x[i]) | TType(0.0));
      i++;

    } else if(v1.p[i] == v2.p[j]-v2.offset) {

      res.p.push_back(v1.p[i]);
      res.x.push_back(TType(v1.x[i] | v2.x[j]));
      i++; j++;

    } else if(v1.p[i] > v2.p[j]-v2.offset) {

      res.p.push_back(v2.p[j]-v2.offset);
      res.x.push_back(TType(0.0) | TType(v2.x[j]));
      j++;

    }

  }

  for( ; i<v1.get_nnz() ; i++) {
    res.p.push_back(v1.p[i]);
    res.x.push_back(TType(v1.x[i]) | TType(0.0));
  }      

  for( ; j<=v2.end ; j++) {
    res.p.push_back(v2.p[j]-v2.offset);
    res.x.push_back(TType(0.0) | TType(v2.x[j]));
  }      

  return res;
}

template<class Tx, class Ty, class Tres>
inline Tres spf_vv_hull(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=VecLen(v2)) cxscthrow(OP_WITH_WRONG_DIM("operator|(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v2);
  int lb = Lb(res);

  for(int i=0 ; i<v1.get_nnz() ; i++) 
    res[v1.p[i]+lb] |= v1.x[i];
  
  return res;
}

template<class Tx, class Ty, class Tres>
inline Tres fsp_vv_hull(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(VecLen(v1)!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator|(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v1);
  int lb = Lb(res);

  for(int i=0 ; i<v2.get_nnz() ; i++) 
    res[v2.p[i]+lb] |= v2.x[i];
  
  return res;
}

template<class Tx, class Ty, class Tres>
inline Tres slf_vv_hull(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=VecLen(v2)) cxscthrow(OP_WITH_WRONG_DIM("operator|(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v2);
  int lb = Lb(res);

  for(int i=v1.start ; i<=v1.end ; i++) 
    res[v1.p[i]-v1.offset+lb] |= v1.x[i];
  
  return res;
}

template<class Tx, class Ty, class Tres>
inline Tres fsl_vv_hull(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(VecLen(v1)!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator|(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v1);
  int lb = Lb(res);

  for(int i=v2.start ; i<=v2.end ; i++) 
    res[v2.p[i]-v2.offset+lb] |= v2.x[i];
  
  return res;

}

template<class Tx, class Ty, class Tres, class TType>
inline Tres slsl_vv_hull(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator|(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v1.n, v1.get_nnz()+v2.get_nnz());
  
  int i=v1.start, j=v2.start;

  while(i<=v1.end && j<=v2.end) {

    if(v1.p[i]-v1.offset < v2.p[j]-v2.offset) {

      res.p.push_back(v1.p[i]-v1.offset);
      res.x.push_back(TType(v1.x[i]) | TType(0.0));
      i++;

    } else if (v1.p[i]-v1.offset == v2.p[j]-v2.offset) {

      res.p.push_back(v2.p[j]-v2.offset);
      res.x.push_back(TType(v1.x[i] | v2.x[j]));
      i++; j++;

    } else if (v1.p[i]-v1.offset > v2.p[j]-v2.offset) {
    
      res.p.push_back(v2.p[j]-v2.offset);
      res.x.push_back(TType(0.0) | TType(v2.x[j]));
      j++;

    }

  }

  for( ; i<=v1.end ; i++) {
    res.p.push_back(v1.p[i]-v1.offset);
    res.x.push_back(TType(v1.x[i]) | TType(0.0));
  }      

  for( ; j<=v2.end ; j++) {
    res.p.push_back(v2.p[j]-v2.offset);
    res.x.push_back(TType(0.0) | TType(v2.x[j]));
  }      

  return res;
}



//----------------------------------------------------------------------------------------

template<class Tx, class Ty, class Tres, class TType>
inline Tres spsp_vv_intersect(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator&(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v1.n, v1.get_nnz()+v2.get_nnz());

  int i=0, j=0;

  while(i<v1.get_nnz() && j<v2.get_nnz()) {
    if(v1.p[i] < v2.p[j]) {
      res.p.push_back(v1.p[i]);
      res.x.push_back(TType(v1.x[i]) & TType(0.0));
      i++;
    } else if (v1.p[i] == v2.p[j]) {
      res.p.push_back(v1.p[i]);
      res.x.push_back(TType(v1.x[i] & v2.x[j]));
      i++; j++;
    } else if (v1.p[i] > v2.p[j]) {
      res.p.push_back(v2.p[j]);
      res.x.push_back(TType(0.0) & TType(v2.x[j]));
      j++;
    }
  }

  for( ; i<v1.get_nnz() ; i++) {
      res.p.push_back(v1.p[i]);
      res.x.push_back(TType(v1.x[i]) & TType(0.0));
  }      

  for( ; j<v2.get_nnz() ; j++) {
      res.p.push_back(v2.p[j]);
      res.x.push_back(TType(0.0) & TType(v2.x[j]));
  }      

  return res;
}


template<class Tx, class Ty, class Tres, class TType>
inline Tres slsp_vv_intersect(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator&(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v1.n, v1.get_nnz()+v2.get_nnz());
  
  int i=v1.start, j=0;

  while(i<=v1.end && j<v2.get_nnz()) {

    if(v1.p[i]-v1.offset < v2.p[j]) {

      res.p.push_back(v1.p[i]-v1.offset);
      res.x.push_back(TType(v1.x[i]) & TType(0.0));
      i++;

    } else if(v1.p[i]-v1.offset == v2.p[j]){

      res.p.push_back(v2.p[j]);
      res.x.push_back(TType(v1.x[i] & v2.x[j]));
      i++; j++;

    } else if(v1.p[i]-v1.offset > v2.p[j]) {

      res.p.push_back(v2.p[j]);
      res.x.push_back(TType(0.0) & TType(v2.x[j]));
      j++;
    }
  }

  for( ; i<=v1.end ; i++) {
    res.p.push_back(v1.p[i]-v1.offset);
    res.x.push_back(TType(v1.x[i]) & TType(0.0));
  }      

  for( ; j<v2.get_nnz() ; j++) {
    res.p.push_back(v2.p[j]);
    res.x.push_back(TType(0.0) & TType(v2.x[j]));
  }      

  return res;
}

template<class Tx, class Ty, class Tres, class TType>
inline Tres spsl_vv_intersect(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator&(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v1.n, v1.get_nnz()+v2.get_nnz());
  
  int i=0, j=v2.start;

  while(i<v1.get_nnz() && j<=v2.end) {

    if(v1.p[i] < v2.p[j]-v2.offset) {

      res.p.push_back(v1.p[i]);
      res.x.push_back(TType(v1.x[i]) & TType(0.0));
      i++;

    } else if(v1.p[i] == v2.p[j]-v2.offset) {

      res.p.push_back(v1.p[i]);
      res.x.push_back(TType(v1.x[i] & v2.x[j]));
      i++; j++;

    } else if(v1.p[i] > v2.p[j]-v2.offset) {

      res.p.push_back(v2.p[j]-v2.offset);
      res.x.push_back(TType(0.0) & TType(v2.x[j]));
      j++;

    }

  }

  for( ; i<v1.get_nnz() ; i++) {
    res.p.push_back(v1.p[i]);
    res.x.push_back(TType(v1.x[i]) & TType(0.0));
  }      

  for( ; j<=v2.end ; j++) {
    res.p.push_back(v2.p[j]-v2.offset);
    res.x.push_back(TType(0.0) & TType(v2.x[j]));
  }      

  return res;
}

template<class Tx, class Ty, class Tres>
inline Tres spf_vv_intersect(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=VecLen(v2)) cxscthrow(OP_WITH_WRONG_DIM("operator&(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v2);
  int lb = Lb(res);

  for(int i=0 ; i<v1.get_nnz() ; i++) 
    res[v1.p[i]+lb] &= v1.x[i];
  
  return res;
}

template<class Tx, class Ty, class Tres>
inline Tres fsp_vv_intersect(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(VecLen(v1)!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator&(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v1);
  int lb = Lb(res);

  for(int i=0 ; i<v2.get_nnz() ; i++) 
    res[v2.p[i]+lb] &= v2.x[i];
  
  return res;
}

template<class Tx, class Ty, class Tres>
inline Tres slf_vv_intersect(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=VecLen(v2)) cxscthrow(OP_WITH_WRONG_DIM("operator&(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v2);
  int lb = Lb(res);

  for(int i=v1.start ; i<=v1.end ; i++) 
    res[v1.p[i]-v1.offset+lb] &= v1.x[i];
  
  return res;
}

template<class Tx, class Ty, class Tres>
inline Tres fsl_vv_intersect(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(VecLen(v1)!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator&(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v1);
  int lb = Lb(res);

  for(int i=v2.start ; i<=v2.end ; i++) 
    res[v2.p[i]-v2.offset+lb] &= v2.x[i];
  
  return res;

}

template<class Tx, class Ty, class Tres, class TType>
inline Tres slsl_vv_intersect(const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator&(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  Tres res(v1.n, v1.get_nnz()+v2.get_nnz());
  
  int i=v1.start, j=v2.start;

  while(i<=v1.end && j<=v2.end) {

    if(v1.p[i]-v1.offset < v2.p[j]-v2.offset) {

      res.p.push_back(v1.p[i]-v1.offset);
      res.x.push_back(TType(v1.x[i]) & TType(0.0));
      i++;

    } else if (v1.p[i]-v1.offset == v2.p[j]-v2.offset) {

      res.p.push_back(v2.p[j]-v2.offset);
      res.x.push_back(TType(v1.x[i] & v2.x[j]));
      i++; j++;

    } else if (v1.p[i]-v1.offset > v2.p[j]-v2.offset) {
    
      res.p.push_back(v2.p[j]-v2.offset);
      res.x.push_back(TType(0.0) & TType(v2.x[j]));
      j++;

    }

  }

  for( ; i<=v1.end ; i++) {
    res.p.push_back(v1.p[i]-v1.offset);
    res.x.push_back(TType(v1.x[i]) & TType(0.0));
  }      

  for( ; j<=v2.end ; j++) {
    res.p.push_back(v2.p[j]-v2.offset);
    res.x.push_back(TType(0.0) & TType(v2.x[j]));
  }      

  return res;
}

//---------------------------------------------------------------------------

template<class Tx, class Ty>
inline bool spsp_vv_comp(const Tx& v1, const Ty& v2) {
  if(v1.n != v2.n) 
    return false;

  int j=0;
  for(int i=0 ; i<v1.get_nnz() ; i++) {
    if(v2.p[j] < v1.p[i]) {
      if(v2.x[j] != 0.0) return false;
      j++;
    } else if(v2.p[j] == v1.p[i]) {
      if(v2.x[j] != v1.x[i]) return false;
      j++;
    } else if(v2.p[j] > v1.p[i]) {
      if(v1.x[i] != 0.0) return false;
    }
  }

  for( ; j<v2.get_nnz() ; j++) {
    if(v2.x[j] != 0.0) return false;
  }

  return true;
}

template<class Tx, class Ty>
inline bool slsp_vv_comp(const Tx& v1, const Ty& v2) {
  if(v1.n != v2.n) 
    return false;

  int j=0;
  for(int i=v1.start ; i<=v1.end ; i++) {
    if(v2.p[j] < v1.p[i]-v1.offset) {
      if(v2.x[j] != 0.0) return false;
      j++;
    } else if(v2.p[j] == v1.p[i]-v1.offset) {
      if(v2.x[j] != v1.x[i]) return false;
      j++;
    } else if(v2.p[j] > v1.p[i]-v1.offset) {
      if(v1.x[i] != 0.0) return false;
    }
  }

  for( ; j<v2.get_nnz() ; j++) {
    if(v2.x[j] != 0.0) return false;
  }

  return true;
}

template<class Tx, class Ty>
inline bool spsl_vv_comp(const Tx& v1, const Ty& v2) {
  if(v1.n != v2.n) 
    return false;

  int j=v2.start;
  for(int i=0 ; i<v1.get_nnz() ; i++) {
    if(v2.p[j]-v2.offset < v1.p[i]) {
      if(v2.x[j] != 0.0) return false;
      j++;
    } else if(v2.p[j]-v2.offset == v1.p[i]) {
      if(v2.x[j] != v1.x[i]) return false;
      j++;
    } else if(v2.p[j]-v2.offset > v1.p[i]) {
      if(v1.x[i] != 0.0) return false;
    }
  }

  for( ; j<=v2.end ; j++) {
    if(v2.x[j] != 0.0) return false;
  }

  return true;
}

template<class Tx, class Ty>
inline bool spf_vv_comp(const Tx& v1, const Ty& v2) {
  if(v1.n != VecLen(v2)) return false;

  int j = 0;
  for(int i=0 ; i<v1.n ; i++) {
    if(v1.p[j] == i) {
      if(v1.x[j] == v2[Lb(v2)+i])
        j++;
      else
        return false;
    } else if(v2[Lb(v2)+i] != 0.0) {
      return false;
    }
  }

  return true;
}

template<class Tx, class Ty>
inline bool fsp_vv_comp(const Tx& v2, const Ty& v1) {
  if(v1.n != VecLen(v2)) return false;

  int j = 0;
  for(int i=0 ; i<v1.n ; i++) {
    if(v1.p[j] == i) {
      if(v1.x[j] == v2[Lb(v2)+i])
        j++;
      else
        return false;
    } else if(v2[Lb(v2)+i] != 0.0) {
      return false;
    }
  }

  return true;
}

template<class Tx, class Ty>
inline bool slf_vv_comp(const Tx& v1, const Ty& v2) {
  if(v1.n != VecLen(v2)) return false;

  int j = v1.start;
  for(int i=0 ; i<v1.n ; i++) {
    if(v1.p[j]-v1.offset == i) {
      if(v1.x[j] == v2[Lb(v2)+i])
        j++;
      else
        return false;
    } else if(v2[Lb(v2)+i] != 0.0) {
      return false;
    }
  }

  return true;
}

template<class Tx, class Ty>
inline bool fsl_vv_comp(const Tx& v2, const Ty& v1) {
  if(v1.n != VecLen(v2)) return false;

  int j = v1.start;
  for(int i=0 ; i<v1.n ; i++) {
    if(v1.p[j]-v1.offset == i) {
      if(v1.x[j] == v2[Lb(v2)+i])
        j++;
      else
        return false;
    } else if(v2[Lb(v2)+i] != 0.0) {
      return false;
    }
  }

  return true;
}

template<class Tx, class Ty>
inline bool slsl_vv_comp(const Tx& v1, const Ty& v2) {
  if(v1.n != v2.n) 
    return false;

  int j=v2.start;
  for(int i=v1.start ; i<=v1.end ; i++) {
    if(v2.p[j]-v2.offset < v1.p[i]-v2.offset) {
      if(v2.x[j] != 0.0) return false;
      j++;
    } else if(v2.p[j]-v2.offset == v1.p[i]-v1.offset) {
      if(v2.x[j] != v1.x[i]) return false;
      j++;
    } else if(v2.p[j]-v2.offset > v1.p[i]-v1.offset) {
      if(v1.x[i] != 0.0) return false;
    }
  }

  for( ; j<=v2.end ; j++) {
    if(v2.x[j] != 0.0) return false;
  }

  return true;

}

//---------------------------------------------------------------------------

template<class Tx, class Ty, class TType>
inline bool spsp_vv_less(const Tx& v1, const Ty& v2) {
  if(v1.n != v2.n) 
    return false;

  int j=0;
  for(int i=0 ; i<v1.get_nnz() ; i++) {
    if(v2.p[j] < v1.p[i]) {
      if(!(TType(0.0) < v2.x[j])) return false;
      j++;
    } else if(v2.p[j] == v1.p[i]) {
      if(!(v1.x[i] < v2.x[j])) return false;
      j++;
    } else if(v2.p[j] > v1.p[i]) {
      if(!(v1.x[i] < TType(0.0))) return false;
    }
  }

  for( ; j<v2.get_nnz() ; j++) {
    if(!(TType(0.0) < v2.x[j])) return false;
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool slsp_vv_less(const Tx& v1, const Ty& v2) {
  if(v1.n != v2.n) 
    return false;

  int j=0;
  for(int i=v1.start ; i<=v1.end ; i++) {
    if(v2.p[j] < v1.p[i]-v1.offset) {
      if(!(TType(0.0) < v2.x[j])) return false;
      j++;
    } else if(v2.p[j] == v1.p[i]-v1.offset) {
      if(!(v1.x[i] < v2.x[j])) return false;
      j++;
    } else if(v2.p[j] > v1.p[i]-v1.offset) {
      if(!(v1.x[i] < TType(0.0))) return false;
    }
  }

  for( ; j<v2.get_nnz() ; j++) {
    if(!(TType(0.0) < v2.x[j])) return false;
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool spsl_vv_less(const Tx& v1, const Ty& v2) {
  if(v1.n != v2.n) 
    return false;

  int j=v2.start;
  for(int i=0 ; i<v1.get_nnz() ; i++) {
    if(v2.p[j]-v2.offset < v1.p[i]) {
      if(!(TType(0.0) < v2.x[j])) return false;
      j++;
    } else if(v2.p[j]-v2.offset == v1.p[i]) {
      if(!(v1.x[i] < v2.x[j])) return false;
      j++;
    } else if(v2.p[j]-v2.offset > v1.p[i]) {
      if(!(v1.x[i] < TType(0.0))) return false;
    }
  }

  for( ; j<=v2.end ; j++) {
    if(!(TType(0.0) < v2.x[j])) return false;
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool spf_vv_less(const Tx& v1, const Ty& v2) {
  if(v1.n != VecLen(v2)) return false;

  int j = 0;
  for(int i=0 ; i<v1.n ; i++) {
    if(v1.p[j] == i) {
      if(v1.x[j] < v2[Lb(v2)+i])
        j++;
      else
        return false;
    } else if(!(TType(0.0) < v2[Lb(v2)+i])) {
      return false;
    }
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool fsp_vv_less(const Tx& v2, const Ty& v1) {
  if(v1.n != VecLen(v2)) return false;

  int j = 0;
  for(int i=0 ; i<v1.n ; i++) {
    if(v1.p[j] == i) {
      if(v1.x[j] < v2[Lb(v2)+i])
        j++;
      else
        return false;
    } else if(!(TType(0.0) < v2[Lb(v2)+i])) {
      return false;
    }
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool slf_vv_less(const Tx& v1, const Ty& v2) {
  if(v1.n != VecLen(v2)) return false;

  int j = v1.start;
  for(int i=0 ; i<v1.n ; i++) {
    if(v1.p[j]-v1.offset == i) {
      if(v1.x[j] < v2[Lb(v2)+i])
        j++;
      else
        return false;
    } else if(!(TType(0.0) < v2[Lb(v2)+i])) {
      return false;
    }
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool fsl_vv_less(const Tx& v2, const Ty& v1) {
  if(v1.n != VecLen(v2)) return false;

  int j = v1.start;
  for(int i=0 ; i<v1.n ; i++) {
    if(v1.p[j]-v1.offset == i) {
      if(v2[Lb(v2)+i] < v1.x[j])
        j++;
      else
        return false;
    } else if(!(v2[Lb(v2)+i]) < TType(0.0)) {
      return false;
    }
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool slsl_vv_less(const Tx& v1, const Ty& v2) {
  if(v1.n != v2.n) 
    return false;

  int j=v2.start;
  for(int i=v1.start ; i<=v1.end ; i++) {
    if(v2.p[j]-v2.offset < v1.p[i]-v2.offset) {
      if(!(TType(0.0) < v2.x[j])) return false;
      j++;
    } else if(v2.p[j]-v2.offset == v1.p[i]-v1.offset) {
      if(!(v1.x[i] < v2.x[j])) return false;
      j++;
    } else if(v2.p[j]-v2.offset > v1.p[i]-v1.offset) {
      if(!(v1.x[i] < TType(0.0))) return false;
    }
  }

  for( ; j<=v2.end ; j++) {
    if(!(TType(0.0) < v2.x[j])) return false;
  }

  return true;

}


//---------------------------------------------------------------------------

template<class Tx, class Ty, class TType>
inline bool spsp_vv_leq(const Tx& v1, const Ty& v2) {
  if(v1.n != v2.n) 
    return false;

  int j=0;
  for(int i=0 ; i<v1.get_nnz() ; i++) {
    if(v2.p[j] < v1.p[i]) {
      if(!(TType(0.0) <= v2.x[j])) return false;
      j++;
    } else if(v2.p[j] == v1.p[i]) {
      if(!(v1.x[i] <= v2.x[j])) return false;
      j++;
    } else if(v2.p[j] > v1.p[i]) {
      if(!(v1.x[i] <= TType(0.0))) return false;
    }
  }

  for( ; j<v2.get_nnz() ; j++) {
    if(!(TType(0.0) <= v2.x[j])) return false;
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool slsp_vv_leq(const Tx& v1, const Ty& v2) {
  if(v1.n != v2.n) 
    return false;

  int j=0;
  for(int i=v1.start ; i<=v1.end ; i++) {
    if(v2.p[j] <= v1.p[i]-v1.offset) {
      if(!(TType(0.0) <= v2.x[j])) return false;
      j++;
    } else if(v2.p[j] == v1.p[i]-v1.offset) {
      if(!(v1.x[i] <= v2.x[j])) return false;
      j++;
    } else if(v2.p[j] > v1.p[i]-v1.offset) {
      if(!(v1.x[i] <= TType(0.0))) return false;
    }
  }

  for( ; j<v2.get_nnz() ; j++) {
    if(!(TType(0.0) <= v2.x[j])) return false;
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool spsl_vv_leq(const Tx& v1, const Ty& v2) {
  if(v1.n != v2.n) 
    return false;

  int j=v2.start;
  for(int i=0 ; i<v1.get_nnz() ; i++) {
    if(v2.p[j]-v2.offset < v1.p[i]) {
      if(!(TType(0.0) <= v2.x[j])) return false;
      j++;
    } else if(v2.p[j]-v2.offset == v1.p[i]) {
      if(!(v1.x[i] <= v2.x[j])) return false;
      j++;
    } else if(v2.p[j]-v2.offset > v1.p[i]) {
      if(!(v1.x[i] <= TType(0.0))) return false;
    }
  }

  for( ; j<=v2.end ; j++) {
    if(!(TType(0.0) <= v2.x[j])) return false;
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool spf_vv_leq(const Tx& v1, const Ty& v2) {
  if(v1.n != VecLen(v2)) return false;

  int j = 0;
  for(int i=0 ; i<v1.n ; i++) {
    if(v1.p[j] == i) {
      if(v1.x[j] <= v2[Lb(v2)+i])
        j++;
      else
        return false;
    } else if(!(TType(0.0) <= v2[Lb(v2)+i])) {
      return false;
    }
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool fsp_vv_leq(const Tx& v2, const Ty& v1) {
  if(v1.n != VecLen(v2)) return false;

  int j = 0;
  for(int i=0 ; i<v1.n ; i++) {
    if(v1.p[j] == i) {
      if(v1.x[j] <= v2[Lb(v2)+i])
        j++;
      else
        return false;
    } else if(!(TType(0.0) <= v2[Lb(v2)+i])) {
      return false;
    }
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool slf_vv_leq(const Tx& v1, const Ty& v2) {
  if(v1.n != VecLen(v2)) return false;

  int j = v1.start;
  for(int i=0 ; i<v1.n ; i++) {
    if(v1.p[j]-v1.offset == i) {
      if(v1.x[j] <= v2[Lb(v2)+i])
        j++;
      else
        return false;
    } else if(!(TType(0.0) <= v2[Lb(v2)+i])) {
      return false;
    }
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool fsl_vv_leq(const Tx& v2, const Ty& v1) {
  if(v1.n != VecLen(v2)) return false;

  int j = v1.start;
  for(int i=0 ; i<v1.n ; i++) {
    if(v1.p[j]-v1.offset == i) {
      if(v2[Lb(v2)+i] <= v1.x[j])
        j++;
      else
        return false;
    } else if(!(v2[Lb(v2)+i] <= TType(0.0))) {
      return false;
    }
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool slsl_vv_leq(const Tx& v1, const Ty& v2) {
  if(v1.n != v2.n) 
    return false;

  int j=v2.start;
  for(int i=v1.start ; i<=v1.end ; i++) {
    if(v2.p[j]-v2.offset < v1.p[i]-v2.offset) {
      if(!(TType(0.0) <= v2.x[j])) return false;
      j++;
    } else if(v2.p[j]-v2.offset == v1.p[i]-v1.offset) {
      if(!(v1.x[i] <= v2.x[j])) return false;
      j++;
    } else if(v2.p[j]-v2.offset > v1.p[i]-v1.offset) {
      if(!(v1.x[i] <= TType(0.0))) return false;
    }
  }

  for( ; j<=v2.end ; j++) {
    if(!(TType(0.0) <= v2.x[j])) return false;
  }

  return true;

}

//--------------------------------------------------------------------------

template<class Tx, class Ty, class TType>
inline bool spsp_vv_greater(const Tx& v1, const Ty& v2) {
  if(v1.n != v2.n) 
    return false;

  int j=0;
  for(int i=0 ; i<v1.get_nnz() ; i++) {
    if(v2.p[j] < v1.p[i]) {
      if(!(TType(0.0) > v2.x[j])) return false;
      j++;
    } else if(v2.p[j] == v1.p[i]) {
      if(!(v1.x[i] > v2.x[j])) return false;
      j++;
    } else if(v2.p[j] > v1.p[i]) {
      if(!(v1.x[i] > TType(0.0))) return false;
    }
  }

  for( ; j<v2.get_nnz() ; j++) {
    if(!(TType(0.0) > v2.x[j])) return false;
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool slsp_vv_greater(const Tx& v1, const Ty& v2) {
  if(v1.n != v2.n) 
    return false;

  int j=0;
  for(int i=v1.start ; i<=v1.end ; i++) {
    if(v2.p[j] < v1.p[i]-v1.offset) {
      if(!(TType(0.0) > v2.x[j])) return false;
      j++;
    } else if(v2.p[j] == v1.p[i]-v1.offset) {
      if(!(v1.x[i] > v2.x[j])) return false;
      j++;
    } else if(v2.p[j] > v1.p[i]-v1.offset) {
      if(!(v1.x[i] > TType(0.0))) return false;
    }
  }

  for( ; j<v2.get_nnz() ; j++) {
    if(!(TType(0.0) > v2.x[j])) return false;
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool spsl_vv_greater(const Tx& v1, const Ty& v2) {
  if(v1.n != v2.n) 
    return false;

  int j=v2.start;
  for(int i=0 ; i<v1.get_nnz() ; i++) {
    if(v2.p[j]-v2.offset < v1.p[i]) {
      if(!(TType(0.0) > v2.x[j])) return false;
      j++;
    } else if(v2.p[j]-v2.offset == v1.p[i]) {
      if(!(v1.x[i] > v2.x[j])) return false;
      j++;
    } else if(v2.p[j]-v2.offset > v1.p[i]) {
      if(!(v1.x[i] > TType(0.0))) return false;
    }
  }

  for( ; j<=v2.end ; j++) {
    if(!(TType(0.0) < v2.x[j])) return false;
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool spf_vv_greater(const Tx& v1, const Ty& v2) {
  if(v1.n != VecLen(v2)) return false;

  int j = 0;
  for(int i=0 ; i<v1.n ; i++) {
    if(v1.p[j] == i) {
      if(v1.x[j] > v2[Lb(v2)+i])
        j++;
      else
        return false;
    } else if(!(TType(0.0) > v2[Lb(v2)+i])) {
      return false;
    }
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool fsp_vv_greater(const Tx& v2, const Ty& v1) {
  if(v1.n != VecLen(v2)) return false;

  int j = 0;
  for(int i=0 ; i<v1.n ; i++) {
    if(v1.p[j] == i) {
      if(v1.x[j] > v2[Lb(v2)+i])
        j++;
      else
        return false;
    } else if(!(TType(0.0) > v2[Lb(v2)+i])) {
      return false;
    }
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool slf_vv_greater(const Tx& v1, const Ty& v2) {
  if(v1.n != VecLen(v2)) return false;

  int j = v1.start;
  for(int i=0 ; i<v1.n ; i++) {
    if(v1.p[j]-v1.offset == i) {
      if(v1.x[j] > v2[Lb(v2)+i])
        j++;
      else
        return false;
    } else if(!(TType(0.0) > v2[Lb(v2)+i])) {
      return false;
    }
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool fsl_vv_greater(const Tx& v2, const Ty& v1) {
  if(v1.n != VecLen(v2)) return false;

  int j = v1.start;
  for(int i=0 ; i<v1.n ; i++) {
    if(v1.p[j]-v1.offset == i) {
      if(v2[Lb(v2)+i] > v1.x[j])
        j++;
      else
        return false;
    } else if(!(v2[Lb(v2)+i] > TType(0.0))) {
      return false;
    }
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool slsl_vv_greater(const Tx& v1, const Ty& v2) {
  if(v1.n != v2.n) 
    return false;

  int j=v2.start;
  for(int i=v1.start ; i<=v1.end ; i++) {
    if(v2.p[j]-v2.offset < v1.p[i]-v2.offset) {
      if(!(TType(0.0) > v2.x[j])) return false;
      j++;
    } else if(v2.p[j]-v2.offset == v1.p[i]-v1.offset) {
      if(!(v1.x[i] > v2.x[j])) return false;
      j++;
    } else if(v2.p[j]-v2.offset > v1.p[i]-v1.offset) {
      if(!(v1.x[i] > TType(0.0))) return false;
    }
  }

  for( ; j<=v2.end ; j++) {
    if(!(TType(0.0) > v2.x[j])) return false;
  }

  return true;

}

//---------------------------------------------------------------------------

template<class Tx, class Ty, class TType>
inline bool spsp_vv_geq(const Tx& v1, const Ty& v2) {
  if(v1.n != v2.n) 
    return false;

  int j=0;
  for(int i=0 ; i<v1.get_nnz() ; i++) {
    if(v2.p[j] < v1.p[i]) {
      if(!(TType(0.0) >= v2.x[j])) return false;
      j++;
    } else if(v2.p[j] == v1.p[i]) {
      if(!(v1.x[i] >= v2.x[j])) return false;
      j++;
    } else if(v2.p[j] > v1.p[i]) {
      if(!(v1.x[i] >= TType(0.0))) return false;
    }
  }

  for( ; j<v2.get_nnz() ; j++) {
    if(!(TType(0.0) >= v2.x[j])) return false;
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool slsp_vv_geq(const Tx& v1, const Ty& v2) {
  if(v1.n != v2.n) 
    return false;

  int j=0;
  for(int i=v1.start ; i<=v1.end ; i++) {
    if(v2.p[j] < v1.p[i]-v1.offset) {
      if(!(TType(0.0) >= v2.x[j])) return false;
      j++;
    } else if(v2.p[j] == v1.p[i]-v1.offset) {
      if(!(v1.x[i] >= v2.x[j])) return false;
      j++;
    } else if(v2.p[j] > v1.p[i]-v1.offset) {
      if(!(v1.x[i] >= TType(0.0))) return false;
    }
  }

  for( ; j<v2.get_nnz() ; j++) {
    if(!(TType(0.0) >= v2.x[j])) return false;
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool spsl_vv_geq(const Tx& v1, const Ty& v2) {
  if(v1.n != v2.n) 
    return false;

  int j=v2.start;
  for(int i=0 ; i<v1.get_nnz() ; i++) {
    if(v2.p[j]-v2.offset < v1.p[i]) {
      if(!(TType(0.0) >= v2.x[j])) return false;
      j++;
    } else if(v2.p[j]-v2.offset == v1.p[i]) {
      if(!(v1.x[i] >= v2.x[j])) return false;
      j++;
    } else if(v2.p[j]-v2.offset > v1.p[i]) {
      if(!(v1.x[i] >= TType(0.0))) return false;
    }
  }

  for( ; j<=v2.end ; j++) {
    if(!(TType(0.0) >= v2.x[j])) return false;
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool spf_vv_geq(const Tx& v1, const Ty& v2) {
  if(v1.n != VecLen(v2)) return false;

  int j = 0;
  for(int i=0 ; i<v1.n ; i++) {
    if(v1.p[j] == i) {
      if(v1.x[j] >= v2[Lb(v2)+i])
        j++;
      else
        return false;
    } else if(!(TType(0.0) >= v2[Lb(v2)+i])) {
      return false;
    }
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool fsp_vv_geq(const Tx& v2, const Ty& v1) {
  if(v1.n != VecLen(v2)) return false;

  int j = 0;
  for(int i=0 ; i<v1.n ; i++) {
    if(v1.p[j] == i) {
      if(v1.x[j] >= v2[Lb(v2)+i])
        j++;
      else
        return false;
    } else if(!(TType(0.0) >= v2[Lb(v2)+i])) {
      return false;
    }
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool slf_vv_geq(const Tx& v1, const Ty& v2) {
  if(v1.n != VecLen(v2)) return false;

  int j = v1.start;
  for(int i=0 ; i<v1.n ; i++) {
    if(v1.p[j]-v1.offset == i) {
      if(v1.x[j] >= v2[Lb(v2)+i])
        j++;
      else
        return false;
    } else if(!(TType(0.0) >= v2[Lb(v2)+i])) {
      return false;
    }
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool fsl_vv_geq(const Tx& v2, const Ty& v1) {
  if(v1.n != VecLen(v2)) return false;

  int j = v1.start;
  for(int i=0 ; i<v1.n ; i++) {
    if(v1.p[j]-v1.offset == i) {
      if(v2[Lb(v2)+i] >= v1.x[j])
        j++;
      else
        return false;
    } else if(!(v2[Lb(v2)+i] >= TType(0.0))) {
      return false;
    }
  }

  return true;
}

template<class Tx, class Ty, class TType>
inline bool slsl_vv_geq(const Tx& v1, const Ty& v2) {
  if(v1.n != v2.n) 
    return false;

  int j=v2.start;
  for(int i=v1.start ; i<=v1.end ; i++) {
    if(v2.p[j]-v2.offset < v1.p[i]-v2.offset) {
      if(!(TType(0.0) >= v2.x[j])) return false;
      j++;
    } else if(v2.p[j]-v2.offset == v1.p[i]-v1.offset) {
      if(!(v1.x[i] >= v2.x[j])) return false;
      j++;
    } else if(v2.p[j]-v2.offset > v1.p[i]-v1.offset) {
      if(!(v1.x[i] >= TType(0.0))) return false;
    }
  }

  for( ; j<=v2.end ; j++) {
    if(!(TType(0.0) >= v2.x[j])) return false;
  }

  return true;

}

//-------------------------------------------------------------------------------

template<class Tx, class Ts, class TType>
inline Tx& sp_vs_assign(Tx& v1, const Ts& s) 
{
  v1.p.clear();
  v1.x.clear();

  if(s != 0.0) {
     for(int i=v1.lb ; i<=v1.ub ; i++) {
        v1.p.push_back(i-v1.lb);
        v1.x.push_back(TType(s));
     }
  }

  return v1;
}

template<class Tx, class Ts, class TType, class TIt>
inline Tx& sl_vs_assign(Tx& v1, const Ts& s) 
{
  v1.x.reserve(v1.p.capacity()+v1.n);
  v1.p.reserve(v1.x.capacity()+v1.n);

  std::vector<int>::iterator p_it = v1.p.erase(v1.p.begin()+v1.start, v1.p.begin()+v1.end+1);
  TIt x_it = v1.x.erase(v1.x.begin()+v1.start, v1.x.begin()+v1.end+1);

  v1.nnz = v1.n;
  v1.end = v1.start+v1.nnz-1;

  for(int i=v1.n-1 ; i>=0 ; i--) {
    v1.p.insert(p_it, i+v1.offset);
    v1.x.insert(x_it, TType(s));
  }

  return v1;
}

template<class Tx, class Ty, class TType>
inline Tx& spf_vv_assign(Tx& v1, const Ty& v2) 
{
  v1.n = VecLen(v2);
  v1.p.clear();
  v1.x.clear();

  for(int i=0 ; i<v1.n ; i++) {
     if(v2[i+Lb(v2)] != 0.0) {
        v1.p.push_back(i);
        v1.x.push_back(TType(v2[i+Lb(v2)]));
     }
  }

  return v1;
}

template<class Tx, class Ty, class TType>
inline Tx& spsl_vv_assign(Tx& v1, const Ty& v2) 
{
  v1.n = v2.n;
  v1.p.clear();
  v1.x.clear();

  for(int i=v2.start ; i<=v2.end ; i++) {
        v1.p.push_back(v2.p[i]);
        v1.x.push_back(TType(v2.x[i]));
  }

  return v1;
}

template<class Tx, class Ty, class TType>
inline Tx& fsp_vv_assign(Tx& v1, const Ty& v2) {
  v1 = Tx(v2.lb,v2.ub);
  v1 = TType(0.0);
  for(int i=0 ; i<v2.get_nnz() ; i++) {
    v1[v2.lb + v2.p[i]] = TType(v2.x[i]);
  }
  return v1;
}

template<class Tx, class Ty, class TType>
inline Tx& fsl_vv_assign(Tx& v1, const Ty& v2) {
  v1 = Tx(v2.lb,v2.ub);
  v1 = TType(0.0);
  for(int i=v2.start ; i<=v2.end ; i++) {
    v1[v2.lb + v2.p[i]] = TType(v2.x[i]);
  }
  return v1;
}

//full slice assignment, dimensions must fit
template<class Tx, class Ty, class TType>
inline Tx& fssp_vv_assign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  v1 = TType(0.0);
  for(int i=0 ; i<v2.get_nnz() ; i++) {
    v1[v2.lb + v2.p[i]] = TType(v2.x[i]);
  }
  return v1;
}

//full slice assignment, dimensions must fit
template<class Tx, class Ty, class TType>
inline Tx& fssl_vv_assign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  v1 = TType(0.0);
  for(int i=v2.start ; i<=v2.end ; i++) {
    v1[v2.lb + v2.p[i]] = TType(v2.x[i]);
  }
  return v1;
}

template<class Tx, class Ty, class TType, class TIt>
inline Tx& slsl_vv_assign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  std::vector<int>::iterator p_it = v1.p.erase(v1.p.begin()+v1.start, v1.p.begin()+v1.end+1);
  TIt x_it = v1.x.erase(v1.x.begin()+v1.start, v1.x.begin()+v1.end+1);

  for(int i=v2.get_nnz()-1 ; i>=0 ; i--) {
    v1.p.insert(p_it, v2.p[i]-v2.offset+v1.offset);
    v1.x.insert(x_it, TType(v2.x[i]));
  }

  v1.nnz = v2.nnz;
  v1.end = v1.start+v1.nnz-1;

  return v1;
}

template<class Tx, class Ty, class TType, class TIt>
inline Tx& slsp_vv_assign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  std::vector<int>::iterator p_it = v1.p.erase(v1.p.begin()+v1.start, v1.p.begin()+v1.end+1);
  TIt x_it = v1.x.erase(v1.x.begin()+v1.start, v1.x.begin()+v1.end+1);

  for(int i=v2.get_nnz()-1 ; i>=0 ; i--) {
    v1.p.insert(p_it, v2.p[i]+v1.offset);
    v1.x.insert(x_it, TType(v2.x[i]));
  }

  v1.nnz = v2.get_nnz();
  v1.end = v1.start+v1.nnz-1;

  return v1;
}

template<class Tx, class Ty, class TType, class TIt>
inline Tx& slf_vv_assign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=VecLen(v2)) cxscthrow(OP_WITH_WRONG_DIM("operator=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  std::vector<int>::iterator p_it = v1.p.erase(v1.p.begin()+v1.start, v1.p.begin()+v1.end+1);
  TIt x_it = v1.x.erase(v1.x.begin()+v1.start, v1.x.begin()+v1.end+1);

  v1.nnz = VecLen(v2);
  v1.end = v1.start+v1.nnz-1;

  for(int i=VecLen(v2)-1 ; i>=0 ; i--) {
    v1.p.insert(p_it, i+v1.offset);
    v1.x.insert(x_it, TType(v2[i+Lb(v2)]));
  }

  return v1;
}

//----------------------------------------------------------------------------------------------

template<class Tx, class TType>
inline std::ostream& sp_v_output(std::ostream& os, const Tx& v) {
  if(ioflags.isset(IOFlags::fullinout)) {
    int j=0;

    for(int i=0 ; i<v.n ; i++) {
      if(j < v.get_nnz()  && i == v.p[j]) {
        os << v.x[j] << std::endl;
        j++;
      } else {
        os << (TType)0.0 << std::endl;
      }
    } 

  } else {
    os << v.n << std::endl;
    os << v.get_nnz() << std::endl;
    for(unsigned int i=0 ; i<v.p.size() ; i++) 
      os << v.p[i] << " ";
    os << std::endl;
    for(unsigned int i=0 ; i<v.x.size() ; i++) 
      os << v.x[i] << " ";
    os << std::endl;
  }

  return os;
}

template<class Tx, class TType>
inline std::ostream& sl_v_output(std::ostream& os, const Tx& v) {
  if(ioflags.isset(IOFlags::fullinout)) {
    int j=v.start;

    for(int i=0 ; i<v.n ; i++) {
      if(j <= v.end  && i == v.p[j]) {
        os << v.x[j] << std::endl;
        j++;
      } else {
        os << (TType)0.0 << std::endl;
      }
    } 

  } else {
    os << v.n << std::endl;
    os << v.get_nnz() << std::endl;
    for(int i=v.start ; i<=v.end ; i++) 
      os << v.p[i] << " ";
    os << std::endl;
    for(int i=v.start ; i<=v.end ; i++) 
      os << v.x[i] << " ";
    os << std::endl;
  }

  return os;
}

template<class Tx, class TType>
inline std::istream& sp_v_input(std::istream& is, Tx& v) {
  if(ioflags.isset(IOFlags::fullinout)) {

    v.p.clear();
    v.x.clear();
    TType tmp;

    for(int i=0 ; i<v.n ; i++) {
      is >> tmp;
      if(tmp != 0.0) {
        v.p.push_back(i);
        v.x.push_back(tmp);
      }
    }
      
  } else {

    v.p.clear();
    v.x.clear();
    TType tmp;
    int ind;

    is >> v.n;
    int nnz;
    is >> nnz;

    for(int i=0 ; i<nnz ; i++) {
      is >> ind;
      v.p.push_back(ind);
    }

    for(int i=0 ; i<nnz ; i++) {
      is >> tmp;
      v.x.push_back(tmp);
    }

  }

  return is;
}

template<class Tx, class TType>
inline std::istream& sl_v_input(std::istream& is, Tx& v) {
  if(ioflags.isset(IOFlags::fullinout)) {

    v.p.erase(v.p.begin()+v.start, v.p.begin()+v.end+1);
    v.x.erase(v.x.begin()+v.start, v.x.begin()+v.end+1);
    v.nnz = 0; 
    TType tmp;

    for(int i=v.start ; i<=v.end ; i++) {
      is >> tmp;
      if(tmp != 0.0) {
        v.p.push_back(i);
        v.x.push_back(tmp);
        v.nnz++;
      }
    }

  } else {

    v.p.erase(v.p.begin()+v.start, v.p.begin()+v.end+1);
    v.x.erase(v.x.begin()+v.start, v.x.begin()+v.end+1);
    TType tmp;
    int ind;

    is >> v.n;
    is >> v.nnz;

    for(int i=0 ; i<v.nnz ; i++) {
      is >> ind;
      v.p.insert(v.p.begin()+v.start+i,ind);
    }

    for(int i=0 ; i<v.nnz ; i++) {
      is >> tmp;
      v.x.insert(v.x.begin()+v.start+i,tmp);
    }


  }

  return is;
}


//---------------------------------------------------------------------------

template<class Tx>
inline Tx sp_v_negative(const Tx& v) {
  Tx res(v);
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x[i] = -res.x[i];
  return res;
}

template<class Tx, class Tres>
inline Tres sl_v_negative(const Tx& v) {
  Tres res(v);
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x[i] = -res.x[i];
  return res;
}

//---------------------------------------------------------------------------

template<class Tx, class Ty>
inline Tx& spf_vv_addassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=VecLen(v2)) cxscthrow(OP_WITH_WRONG_DIM("operator+=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  v1 = Tx(v1+v2);
  return v1;
}

template<class Tx, class Ty>
inline Tx& spsp_vv_addassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator+=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  v1 = Tx(v1+v2);
  return v1;
}

template<class Tx, class Ty>
inline Tx& spsl_vv_addassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator+=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  v1 = Tx(v1+v2);
  return v1;
}

template<class Tx, class Ty>
inline Tx& fsp_vv_addassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(VecLen(v1)!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator+=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  for(int i=0 ; i<v2.get_nnz() ; i++) 
    v1[v2.p[i]+Lb(v1)] += v2.x[i];
  
  return v1;
}

template<class Tx, class Ty>
inline Tx& fsl_vv_addassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(VecLen(v1)!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator+=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  for(int i=v2.start ; i<=v2.end ; i++) 
    v1[v2.p[i]-v2.offset+Lb(v1)] += v2.x[i];
  
  return v1;
}

template<class Tx, class Ty>
inline Tx& spf_vv_subassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=VecLen(v2)) cxscthrow(OP_WITH_WRONG_DIM("operator-=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  v1 = Tx(v1-v2);
  return v1;
}

template<class Tx, class Ty>
inline Tx& spsp_vv_subassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator-=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  v1 = Tx(v1-v2);
  return v1;
}

template<class Tx, class Ty>
inline Tx& spsl_vv_subassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator-=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  v1 = Tx(v1-v2);
  return v1;
}

template<class Tx, class Ty>
inline Tx& fsp_vv_subassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(VecLen(v1)!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator-=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  for(int i=0 ; i<v2.get_nnz() ; i++) 
    v1[v2.p[i]+Lb(v1)] -= v2.x[i];
  
  return v1;
}

template<class Tx, class Ty>
inline Tx& fsl_vv_subassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(VecLen(v1)!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator-=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  for(int i=v2.start ; i<=v2.end ; i++) 
    v1[v2.p[i]-v2.offset+Lb(v1)] -= v2.x[i];
  
  return v1;
}

template<class Tx, class Ty>
inline Tx& spf_vv_hullassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=VecLen(v2)) cxscthrow(OP_WITH_WRONG_DIM("operator|=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  v1 = Tx(v1|v2);
  return v1;
}

template<class Tx, class Ty>
inline Tx& spsp_vv_hullassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator|=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  v1 = Tx(v1|v2);
  return v1;
}

template<class Tx, class Ty>
inline Tx& spsl_vv_hullassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator|=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  v1 = Tx(v1|v2);
  return v1;
}

template<class Tx, class Ty>
inline Tx& fsp_vv_hullassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(VecLen(v1)!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator|=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  for(int i=0 ; i<v2.get_nnz() ; i++) 
    v1[v2.p[i]+Lb(v1)] |= v2.x[i];
  
  return v1;
}

template<class Tx, class Ty>
inline Tx& fsl_vv_hullassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(VecLen(v1)!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator|=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  for(int i=v2.start ; i<=v2.end ; i++) 
    v1[v2.p[i]-v2.offset+Lb(v1)] |= v2.x[i];
  
  return v1;
}

template<class Tx, class Ty>
inline Tx& spf_vv_intersectassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=VecLen(v2)) cxscthrow(OP_WITH_WRONG_DIM("operator&=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  v1 = Tx(v1&v2);
  return v1;
}

template<class Tx, class Ty>
inline Tx& spsp_vv_intersectassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator&=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  v1 = Tx(v1&v2);
  return v1;
}

template<class Tx, class Ty>
inline Tx& spsl_vv_intersectassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator&=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  v1 = Tx(v1&v2);
  return v1;
}

template<class Tx, class Ty>
inline Tx& fsp_vv_intersectassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(VecLen(v1)!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator&=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  for(int i=0 ; i<v2.get_nnz() ; i++) 
    v1[v2.p[i]+Lb(v1)] &= v2.x[i];
  
  return v1;
}

template<class Tx, class Ty>
inline Tx& fsl_vv_intersectassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(VecLen(v1)!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator&=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  for(int i=v2.start ; i<=v2.end ; i++) 
    v1[v2.p[i]-v2.offset+Lb(v1)] &= v2.x[i];
  
  return v1;
}

template<class Tx, class Ty>
inline Tx& sp_vs_multassign(Tx& v, const Ty& s) {
  for(int i=0 ; i<v.get_nnz() ; i++)
    v.x[i] *= s;
  return v;
}

template<class Tx, class Ty>
inline Tx& sp_vs_divassign(Tx& v, const Ty& s) {
  for(int i=0 ; i<v.get_nnz() ; i++)
    v.x[i] /= s;
  return v;
}


//---------------------------------------------------------------------------

template<class Tx, class Ty, class TType>
inline Tx& slf_vv_addassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=VecLen(v2)) cxscthrow(OP_WITH_WRONG_DIM("operator+=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  int j = v1.start;
 
  for(int i=0 ; i<v1.n ; i++) {
    if(v1.p[j]-v1.offset == i) {
      v1.x[j] += v2[i+Lb(v2)];
      j++;
    } else if(v2[i+Lb(v2)] != 0.0) {
      v1.p.insert(v1.p.begin()+j, i);
      v1.x.insert(v1.x.begin()+j, TType(v2[i+Lb(v2)]));
      v1.nnz++;
      v1.end++;
    }
  }

  return v1;
}

template<class Tx, class Ty>
inline Tx& slsp_vv_addassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator+=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  v1 = v1+v2;
  return v1;
}

template<class Tx, class Ty>
inline Tx& slsl_vv_addassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator+=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  v1 = v1+v2;
  return v1;
}

template<class Tx, class Ty, class TType>
inline Tx& slf_vv_subassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=VecLen(v2)) cxscthrow(OP_WITH_WRONG_DIM("operator-=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  int j = v1.start;
 
  for(int i=0 ; i<v1.n ; i++) {
    if(v1.p[j]-v1.offset == i) {
      v1.x[j] -= v2[i+Lb(v2)];
      j++;
    } else if(v2[i+Lb(v2)] != 0.0) {
      v1.p.insert(v1.p.begin()+j, i);
      v1.x.insert(v1.x.begin()+j, TType(-v2[i+Lb(v2)]));
      v1.nnz++;
      v1.end++;
    }
  }

  return v1;
}

template<class Tx, class Ty>
inline Tx& slsp_vv_subassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator-=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  v1 = v1-v2;
  return v1;
}

template<class Tx, class Ty>
inline Tx& slsl_vv_subassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator-=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  v1 = v1-v2;
  return v1;
}

template<class Tx, class Ty>
inline Tx& sl_vs_multassign(Tx& v, const Ty& s) {
  for(int i=v.start ; i<=v.end ; i++)
    v.x[i] *= s;

  return v;
}

template<class Tx, class Ty>
inline Tx& sl_vs_divassign(Tx& v, const Ty& s) {
  for(int i=v.start ; i<=v.end ; i++)
    v.x[i] = s;

  return v;
}


template<class Tx, class Ty, class TType>
inline Tx& slf_vv_hullassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=VecLen(v2)) cxscthrow(OP_WITH_WRONG_DIM("operator|=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  int j = v1.start;
 
  for(int i=0 ; i<v1.n ; i++) {
    if(v1.p[j]-v1.offset == i) {
      v1.x[j] |= v2[i+Lb(v2)];
      j++;
    } else {
      v1.p.insert(v1.p.begin()+j, i);
      v1.x.insert(v1.x.begin()+j, TType(0.0) | TType(v2[i+Lb(v2)]));
      v1.nnz++;
      v1.end++;
    }
  }

  return v1;
}

template<class Tx, class Ty>
inline Tx& slsp_vv_hullassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator|=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  v1 = v1 | v2;
  return v1;
}

template<class Tx, class Ty>
inline Tx& slsl_vv_hullassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator|=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  v1 = v1 | v2;
  return v1;
}

template<class Tx, class Ty, class TType>
inline Tx& slf_vv_intersectassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=VecLen(v2)) cxscthrow(OP_WITH_WRONG_DIM("operator&=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  int j = v1.start;
 
  for(int i=0 ; i<v1.n ; i++) {
    if(v1.p[j]-v1.offset == i) {
      v1.x[j] &= v2[i+Lb(v2)];
      j++;
    } else {
      v1.p.insert(v1.p.begin()+j, i);
      v1.x.insert(v1.x.begin()+j, TType(0.0) & TType(v2[i+Lb(v2)]));
      v1.nnz++;
      v1.end++;
    }
  }

  return v1;
}

template<class Tx, class Ty>
inline Tx& slsp_vv_intersectassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator&=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  v1 = v1 & v2;
  return v1;
}

template<class Tx, class Ty>
inline Tx& slsl_vv_intersectassign(Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("operator&=(const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  v1 = v1 & v2;
  return v1;
}

//--------------------------------------------------------------

template<class Tx>
inline bool sp_v_not(const Tx& v) {
  bool ret = true;
  for(int i=0 ; i<v.get_nnz() ; i++)
    ret = ret && (!v.x[i]);
  return ret;
}

template<class Tx>
inline bool sl_v_not(const Tx& v) {
  bool ret = true;
  for(int i=v.start ; i<=v.end ; i++)
    ret = ret && (!v.x[i]);
  return ret;
}

//--------------------------------------------------------------

template<class TDot, class Tx, class Ty, class TSparseDot>
inline void spsp_vv_accu(TDot& dot, const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
  throw(OP_WITH_WRONG_DIM)
#else
  throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("accumulate(" + nameof(dot) + "&, const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  TSparseDot sdot(dot.get_k());

  for(int i=0 ; i<v1.get_nnz() ; i++) {
    for(int j=0 ; j<v2.get_nnz() && v2.p[j]<=v1.p[i] ; j++) {
       if(v1.p[i] == v2.p[j])
         sdot.add_dot_err(v1.x[i], v2.x[j]);
    }
  }

  sdot.result(dot);
}

template<class TDot, class Tx, class Ty, class TSparseDot>
inline void spf_vv_accu(TDot& dot, const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
  throw(OP_WITH_WRONG_DIM)
#else
  throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=VecLen(v2)) cxscthrow(OP_WITH_WRONG_DIM("accumulate(" + nameof(dot) + "&, const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  TSparseDot sdot(opdotprec);
  int lb = Lb(v2);

  for(int i=0 ; i<v1.get_nnz() ; i++)
    sdot.add_dot_err(v1.x[i], v2[v1.p[i]+lb]);

  sdot.result(dot);
}

template<class TDot, class Tx, class Ty, class TSparseDot>
inline void fsp_vv_accu(TDot& dot, const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
  throw(OP_WITH_WRONG_DIM)
#else
  throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(VecLen(v1)!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("accumulate(" + nameof(dot) + "&, const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  TSparseDot sdot(opdotprec);
  int lb = Lb(v1);

  for(int i=0 ; i<v2.get_nnz() ; i++)
    sdot.add_dot_err(v1[v2.p[i]+lb], v2.x[i]);

  sdot.result(dot);
}

template<class TDot, class Tv1, class Tv2, class TSparseDot>
inline void slsl_vv_accu(TDot& dot, const Tv1& v1, const Tv2& v2) 
#if(CXSC_INDEX_CHECK)
  throw(OP_WITH_WRONG_DIM)
#else
  throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("accumulate(" + nameof(dot) + "&, const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  TSparseDot sdot(opdotprec);

  for(int i=v1.start ; i<=v1.end ; i++) {
    for(int j=v2.start ; j<v2.end+1 && v2.p[j]-v2.offset<=v1.p[i]-v1.offset ; j++) {
       if(v1.p[i]-v1.offset == v2.p[j]-v2.offset)
         sdot.add_dot_err(v1.x[i], v2.x[j]);
    }
  }

  sdot.result(dot);

}

template<class TDot, class Tv1, class Tv2, class TSparseDot>
inline void spsl_vv_accu(TDot& dot, const Tv1& v1, const Tv2& v2) 
#if(CXSC_INDEX_CHECK)
  throw(OP_WITH_WRONG_DIM)
#else
  throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("accumulate(" + nameof(dot) + "&, const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  TSparseDot sdot(opdotprec);

  for(int i=0 ; i<v1.get_nnz() ; i++) {
    for(int j=v2.start ; j<v2.end+1 && v2.p[j]-v2.offset<=v1.p[i] ; j++) {
       if(v1.p[i] == v2.p[j]-v2.offset)
         sdot.add_dot_err(v1.x[i], v2.x[j]);
    }
  }

  sdot.result(dot);
}

template<class TDot, class Tv1, class Tv2, class TSparseDot>
inline void slsp_vv_accu(TDot& dot, const Tv1& v1, const Tv2& v2) 
#if(CXSC_INDEX_CHECK)
  throw(OP_WITH_WRONG_DIM)
#else
  throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("accumulate(" + nameof(dot) + "&, const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  TSparseDot sdot(opdotprec);

  for(int i=v1.start ; i<=v1.end ; i++) {
    for(int j=0 ; j<v2.get_nnz() && v2.p[j]<=v1.p[i]-v1.offset ; j++) {
       if(v1.p[i]-v1.offset == v2.p[j])
         sdot.add_dot_err(v1.x[i], v2.x[j]);
    }
  }

  sdot.result(dot);
}

template<class TDot, class Tv1, class Tv2, class TSparseDot>
inline void slf_vv_accu(TDot& dot, const Tv1& v1, const Tv2& v2) 
#if(CXSC_INDEX_CHECK)
  throw(OP_WITH_WRONG_DIM)
#else
  throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=VecLen(v2)) cxscthrow(OP_WITH_WRONG_DIM("accumulate(" + nameof(dot) + "&, const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  TSparseDot sdot(opdotprec);
  int lb = Lb(v2);

  for(int i=v1.start ; i<=v1.end ; i++)
    sdot.add_dot_err(v1.x[i], v2[v1.p[i]-v1.offset+lb]);

  sdot.result(dot);
}

template<class TDot, class Tv1, class Tv2, class TSparseDot>
inline void fsl_vv_accu(TDot& dot, const Tv1& v1, const Tv2& v2) 
#if(CXSC_INDEX_CHECK)
  throw(OP_WITH_WRONG_DIM)
#else
  throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(VecLen(v1)!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("accumulate(" + nameof(dot) + "&, const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  TSparseDot sdot(opdotprec);
  int lb = Lb(v1);

  for(int i=v2.start ; i<=v2.end ; i++)
    sdot.add_dot_err(v1[v2.p[i]-v2.offset+lb], v2.x[i]);

  sdot.result(dot);
}


template<class TDot, class Tx, class Ty, class TSparseDot>
inline void spsp_vv_accuapprox(TDot& dot, const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
  throw(OP_WITH_WRONG_DIM)
#else
  throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("accumulate_approx(" + nameof(dot) + "&, const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  TSparseDot sdot(dot.get_k());

  for(int i=0 ; i<v1.get_nnz() ; i++) {
    for(int j=0 ; j<v2.get_nnz() && v2.p[j]<=v1.p[i] ; j++) {
       if(v1.p[i] == v2.p[j])
         sdot.add_dot(v1.x[i], v2.x[j]);
    }
  }

  sdot.result(dot);
}

template<class TDot, class Tx, class Ty, class TSparseDot>
inline void spf_vv_accuapprox(TDot& dot, const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
  throw(OP_WITH_WRONG_DIM)
#else
  throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=VecLen(v2)) cxscthrow(OP_WITH_WRONG_DIM("accumulate_approx(" + nameof(dot) + "&, const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  TSparseDot sdot(opdotprec);
  int lb = Lb(v2);

  for(int i=0 ; i<v1.get_nnz() ; i++)
    sdot.add_dot(v1.x[i], v2[v1.p[i]+lb]);

  sdot.result(dot);
}

template<class TDot, class Tx, class Ty, class TSparseDot>
inline void fsp_vv_accuapprox(TDot& dot, const Tx& v1, const Ty& v2) 
#if(CXSC_INDEX_CHECK)
  throw(OP_WITH_WRONG_DIM)
#else
  throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(VecLen(v1)!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("accumulate_approx(" + nameof(dot) + "&, const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  TSparseDot sdot(opdotprec);
  int lb = Lb(v1);

  for(int i=0 ; i<v2.get_nnz() ; i++)
    sdot.add_dot(v1[v2.p[i]+lb], v2.x[i]);

  sdot.result(dot);
}

template<class TDot, class Tv1, class Tv2, class TSparseDot>
inline void slsl_vv_accuapprox(TDot& dot, const Tv1& v1, const Tv2& v2) 
#if(CXSC_INDEX_CHECK)
  throw(OP_WITH_WRONG_DIM)
#else
  throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("accumulate_approx(" + nameof(dot) + "&, const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  TSparseDot sdot(opdotprec);

  for(int i=v1.start ; i<=v1.end ; i++) {
    for(int j=v2.start ; j<v2.end+1 && v2.p[j]-v2.offset<=v1.p[i]-v1.offset ; j++) {
       if(v1.p[i]-v1.offset == v2.p[j]-v2.offset)
         sdot.add_dot(v1.x[i], v2.x[j]);
    }
  }

  sdot.result(dot);

}

template<class TDot, class Tv1, class Tv2, class TSparseDot>
inline void spsl_vv_accuapprox(TDot& dot, const Tv1& v1, const Tv2& v2) 
#if(CXSC_INDEX_CHECK)
  throw(OP_WITH_WRONG_DIM)
#else
  throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("accumulate_approx(" + nameof(dot) + "&, const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  TSparseDot sdot(opdotprec);

  for(int i=0 ; i<v1.get_nnz() ; i++) {
    for(int j=v2.start ; j<v2.end+1 && v2.p[j]-v2.offset<=v1.p[i] ; j++) {
       if(v1.p[i] == v2.p[j]-v2.offset)
         sdot.add_dot(v1.x[i], v2.x[j]);
    }
  }

  sdot.result(dot);
}

template<class TDot, class Tv1, class Tv2, class TSparseDot>
inline void slsp_vv_accuapprox(TDot& dot, const Tv1& v1, const Tv2& v2) 
#if(CXSC_INDEX_CHECK)
  throw(OP_WITH_WRONG_DIM)
#else
  throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("accumulate_approx(" + nameof(dot) + "&, const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  TSparseDot sdot(opdotprec);

  for(int i=v1.start ; i<=v1.end ; i++) {
    for(int j=0 ; j<v2.get_nnz() && v2.p[j]<=v1.p[i]-v1.offset ; j++) {
       if(v1.p[i]-v1.offset == v2.p[j])
         sdot.add_dot(v1.x[i], v2.x[j]);
    }
  }

  sdot.result(dot);
}

template<class TDot, class Tv1, class Tv2, class TSparseDot>
inline void slf_vv_accuapprox(TDot& dot, const Tv1& v1, const Tv2& v2) 
#if(CXSC_INDEX_CHECK)
  throw(OP_WITH_WRONG_DIM)
#else
  throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(v1.n!=VecLen(v2)) cxscthrow(OP_WITH_WRONG_DIM("accumulate_approx(" + nameof(dot) + "&, const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  TSparseDot sdot(opdotprec);
  int lb = Lb(v2);

  for(int i=v1.start ; i<=v1.end ; i++)
    sdot.add_dot(v1.x[i], v2[v1.p[i]-v1.offset+lb]);

  sdot.result(dot);
}

template<class TDot, class Tv1, class Tv2, class TSparseDot>
inline void fsl_vv_accuapprox(TDot& dot, const Tv1& v1, const Tv2& v2) 
#if(CXSC_INDEX_CHECK)
  throw(OP_WITH_WRONG_DIM)
#else
  throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(VecLen(v1)!=v2.n) cxscthrow(OP_WITH_WRONG_DIM("accumulate_approx(" + nameof(dot) + "&, const "+nameof(v1)+" &, const "+nameof(v2)+" &)"));
#endif
  TSparseDot sdot(opdotprec);
  int lb = Lb(v1);

  for(int i=v2.start ; i<=v2.end ; i++)
    sdot.add_dot(v1[v2.p[i]-v2.offset+lb], v2.x[i]);

  sdot.result(dot);
}


//---------------------------------------------------------------------------------------

template<class Tx>
inline void sp_v_resize(Tx &v) throw() {
  v.n = v.ub = 0;
  v.lb = 0;
  v.p.clear();
  v.x.clear();
}

template <class Tx>
inline void sp_v_resize(Tx &v, const int &len)
#if(CXSC_INDEX_CHECK)
    throw(WRONG_BOUNDARIES)
#else
    throw()
#endif
{
   if(v.n == len)
     SetLb(v,1);
   else {
#if(CXSC_INDEX_CHECK)
     if(len<0) cxscthrow(WRONG_BOUNDARIES(" Resize("+nameof(v)+" &, const int &)"));
#endif
     for(int i=0 ; i<v.get_nnz() ; i++) {
       if(v.p[i] >= len) {
         v.p.erase(v.p.begin()+i, v.p.end());
         v.x.erase(v.x.begin()+i, v.x.end());
         break;
       }
     }
     v.lb = 1;
     v.n = v.ub = len;
   }
}

template<class Tx>
inline void sp_v_resize(Tx &v, const int &lb, const int &ub)
#if(CXSC_INDEX_CHECK)
	throw(WRONG_BOUNDARIES)
#else
	throw()
#endif
{
     if(v.n == ub-lb+1)
       SetUb(v,ub);
     else {
#if(CXSC_INDEX_CHECK)
       if(ub<lb) cxscthrow(WRONG_BOUNDARIES("void Resize("+nameof(v)+" &, const int &, const int &)"));
#endif
       sp_v_resize(v, ub-lb+1);
       SetLb(v, lb);
     }
}

//---------------------------------------------------------------------------------


} //namespace cxsc

#endif 
