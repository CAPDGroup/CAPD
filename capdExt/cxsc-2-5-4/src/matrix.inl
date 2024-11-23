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

/* CVS $Id: matrix.inl,v 1.40 2014/01/30 17:23:47 cxsc Exp $ */

#ifndef _CXSC_MATRIX_INL_INCLUDED
#define _CXSC_MATRIX_INL_INCLUDED

#include <iostream>
#include <sstream>
#include <cctype>
#include <string>
#include <cinterval.hpp>
#include "vector.inl"

#ifdef CXSC_USE_BLAS
#include "cxsc_blas.hpp"
#endif

namespace cxsc {

class l_real;
class l_interval;
	
template <class S,class M>
TINLINE void _smconstr(S &s,const M &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__TYPE_CAST_OF_THICK_OBJ<M>,ERROR__USE_OF_UNINITIALIZED_OBJ<M>)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if((m.xsize>1)||(m.ysize>1)) cxscthrow(ERROR__TYPE_CAST_OF_THICK_OBJ<M>(nameof(s)+"::"+nameof(s)+"(const "+nameof(m)+" &)"));
	if((m.xsize<1)||(m.ysize<1)) cxscthrow(ERROR__USE_OF_UNINITIALIZED_OBJ<M>(nameof(s)+"::"+nameof(s)+"(const "+nameof(m)+" &)"));
#endif
	s=m.dat[0];
}

template <class V,class M,class S>
TINLINE void _vmconstr(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__TYPE_CAST_OF_THICK_OBJ<M>)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(m.xsize!=1&&m.ysize!=1) cxscthrow(ERROR__TYPE_CAST_OF_THICK_OBJ<M>(nameof(v)+"::"+nameof(v)+"(const "+nameof(m)+" &)"));
#endif
	if(m.ysize==1)
	{
		v.dat=new S[m.xsize];
		for(int i=0;i<m.xsize;i++)
			v.dat[i]=m.dat[i];
		v.l=m.lb2;
		v.u=m.ub2;
		v.size=m.xsize;
	}
	else
	{
		v.dat=new S[m.ysize];
		for(int i=0;i<m.ysize;i++)
			v.dat[i]=m.dat[i];
		v.l=m.lb1;
		v.u=m.ub1;
		v.size=m.ysize;
	}
}

template <class V,class MS,class S>
TINLINE void _vmsconstr(V &v,const MS &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__TYPE_CAST_OF_THICK_OBJ<MS>)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(m.sxsize!=1&&m.sysize!=1) cxscthrow(ERROR__TYPE_CAST_OF_THICK_OBJ<MS>(nameof(S())+"::"+nameof(S())+"(const "+nameof(m)+" &)"));
#endif
	if(m.sysize==1)
	{
		v.dat=new S[m.sxsize];
		for(int i=0;i<m.sxsize;i++)
			v.dat[i]=m.dat[m.offset1*m.mxsize+i+m.offset2];
		v.l=m.start2;
		v.u=m.end2;
		v.size=m.sxsize;
	}
	else
	{
		v.dat=new S[m.sysize];
		for(int i=0;i<m.sysize;i++)
			v.dat[i]=m.dat[(i+m.offset1)*m.mxsize+m.offset2];
		v.l=m.start1;
		v.u=m.end1;
		v.size=m.sysize;
	}
}

enum MM_FORMAT {mm_coordinate, mm_array};
enum MM_FIELD {mm_real, mm_integer, mm_complex, mm_pattern, mm_interval, mm_cinterval};
enum MM_SYMMETRY {mm_general, mm_symmetric, mm_skew_symmetric, mm_hermitian};

inline string ElementName(const intmatrix& A) {
  return "integer";
}

inline string ElementName(const rmatrix& A) {
  return "real";
}

inline string ElementName(const l_rmatrix& A) {
  return "real";
}

inline string ElementName(const cmatrix& A) {
  return "complex";
}

inline string ElementName(const imatrix& A) {
  return "interval";
}

inline string ElementName(const l_imatrix& A) {
  return "interval";
}

inline string ElementName(const cimatrix& A) {
  return "cinterval";
}


inline real MatrixMarketElement(const real& x) {
  return x;
}

inline real MatrixMarketElement(const l_real& x) {
  return _real(x);
}

inline string MatrixMarketElement(const interval& x) {
  std::stringstream temp;
  temp << SaveOpt << RndDown << Inf(x) << " " 
       << RndUp << Sup(x) << RestoreOpt;
  return temp.str();
}

inline string MatrixMarketElement(const l_interval& x) {
  std::stringstream temp;
  temp << SaveOpt << RndDown << _real(Inf(x)) << " " 
       << RndUp << _real(Sup(x)) << RestoreOpt;
  return temp.str();
}

inline string MatrixMarketElement(const complex& x) {
  std::stringstream temp;
  temp << Re(x) << " " << Im(x);
  return temp.str();
}

inline string MatrixMarketElement(const cinterval& x) {
  std::stringstream temp;
  temp << SaveOpt << RndDown << InfRe(x) << " " 
       << RndUp << SupRe(x) << " " 
       << RndDown << InfIm(x) << " " 
       << RndUp << SupIm(x) << RestoreOpt;
  return temp.str();
}


template <class M>
std::ostream &_mout(std::ostream &s,const M &r) throw()
{
	if(ioflags.isset(IOFlags::matrixmarketinout)) {
		s << "%%MatrixMarket matrix array " << ElementName(r) <<" general" << std::endl;
		s << "%Generated by C-XSC" << std::endl;
		s << r.ysize << " " << r.xsize << std::endl;
		for (int j=0;j<r.xsize;j++)
		{
			for (int i=0;i<r.ysize;i++)
			{
				s << MatrixMarketElement(r.dat[i*r.xsize+j]) << std::endl;
			}
		}
		return s;
	} else {
		int i,j;
		for (i=0;i<r.ysize;i++)
		{
			for (j=0;j<r.xsize;j++)
			{
				s << r.dat[i*r.xsize+j]<<" ";
			}
			s<<std::endl;
		}
		return s;
	}
}

inline void toLower(string& s) {
  for(unsigned int i=0 ; i<s.size() ; i++)
    s[i] = tolower(s[i]);
}

inline int conj_el(const int& r) {
  return r;
}

inline real conj_el(const real& r) {
  return r;
}

inline l_real conj_el(const l_real& r) {
  return r;
}

inline interval conj_el(const interval& r) {
  return r;
}

inline l_interval conj_el(const l_interval& r) {
  return r;
}

inline complex conj_el(const complex& r) {
  return conj(r);
}

inline cinterval conj_el(const cinterval& r) {
  return conj(r);
}


inline void mm_add_element(int& i, std::string& tmp, MM_FIELD& field) {
  std::stringstream ss(tmp);
  real r;
  ss >> r;
  i = (int)_double(r);
}

inline void mm_add_element(real& r, std::string& tmp, MM_FIELD& field) {
  //std::stringstream ss(tmp);
  tmp >> r;
}

inline void mm_add_element(l_real& r, std::string& tmp, MM_FIELD& field) {
  //std::stringstream ss(tmp);
  //ss >> r;
  tmp >> r;
}

inline void mm_add_element(complex& c, std::string& tmp, MM_FIELD& field) {
  std::stringstream ss(tmp);
  real val_re, val_im;

  if(field == mm_real || field == mm_integer || field == mm_interval) {
    ss >> val_re;
    c = val_re;
  } else if(field == mm_complex || field == mm_cinterval) {
    ss >> val_re >> val_im;
    c = complex(val_re,val_im);
  }
}


inline void mm_add_element(interval& i, std::string& tmp, MM_FIELD& field) {
  std::stringstream ss(tmp);
  real val_inf, val_sup;

  if(field == mm_interval || field == mm_cinterval) {
    ss >> val_inf >> val_sup;
    i = interval(val_inf,val_sup);
  } else {
    ss >> val_inf;
    i = val_inf;
  }
}

inline void mm_add_element(l_interval& i, std::string& tmp, MM_FIELD& field) {
  std::stringstream ss(tmp);
  l_real val_inf, val_sup;

  if(field == mm_interval || field == mm_cinterval) {
    ss >> val_inf >> val_sup;
    i = l_interval(val_inf,val_sup);
  } else {
    ss >> val_inf;
    i = val_inf;
  }
}


inline void mm_add_element(cinterval& ci, std::string& tmp, MM_FIELD& field) {
  std::stringstream ss(tmp);
  real val_inf_re, val_sup_re, val_inf_im, val_sup_im;

  if(field == mm_real || field == mm_integer) {
    ss >> val_inf_re;
    ci = val_inf_re;
  } else if(field == mm_complex) {
    ss >> val_inf_re >> val_inf_im;
    ci = complex(val_inf_re,val_inf_im);
  } else if(field == mm_interval) {
    ss >> val_inf_re >> val_sup_re;
    ci = interval(val_inf_re,val_sup_re);
  } else if(field == mm_cinterval) {
    ss >> val_inf_re >> val_sup_re >> val_inf_im >> val_sup_im;
    ci = cinterval(interval(val_inf_re,val_sup_re),interval(val_inf_im,val_sup_im));
  }
}


template <class M>
std::istream &_min(std::istream &s,M &r) throw()
{
	if(ioflags.isset(IOFlags::matrixmarketinout)) {
		//MatrixMarket
		std::string header;
		std::getline(s,header);
		std::stringstream ss(header);
		std::string tmp, s_format, s_field, s_symmetry;
//		MM_FORMAT format;
		MM_FIELD field;
		MM_SYMMETRY symmetry;

		//Read start tag
		ss >> tmp;
		toLower(tmp);
		if(tmp != "%%matrixmarket") return s;

		//Read object tag
		ss >> tmp;
		toLower(tmp);
		if(tmp != "matrix") return s;

		//Read data format tag
		ss >> s_format;
		toLower(s_format);
		if(s_format != "array" && s_format != "coordinate")
			return s;


		//Read number field tag
		ss >> s_field;
		toLower(s_field);
		if(s_field == "real")
			field = mm_real;
		else if(s_field == "integer")
			field = mm_integer;
		else if(s_field == "complex")
			field = mm_complex;
		else if(s_field == "pattern") 
			field = mm_pattern;
		else if(s_field == "interval") 
			field = mm_interval;
		else if(s_field == "cinterval") 
			field = mm_cinterval;
		else 
			return s;

		//Read symmetry tag
		ss >> s_symmetry;
		toLower(s_symmetry);
		if(s_symmetry == "general")
			symmetry = mm_general;
		else if(s_symmetry == "symmetric")
			symmetry = mm_symmetric;
		else if(s_symmetry == "skew-symmetric")
			symmetry = mm_skew_symmetric;
		else if(s_symmetry == "hermitian")
			symmetry = mm_hermitian;
		else
			return s;

		//Read past commentary
		std::getline(s,tmp);
		while(tmp[0] == '%')
			std::getline(s,tmp);

		//Read past empty lines
		while(tmp.size() == 0)
			std::getline(s,tmp);


		//Read size information
		ss.clear();
		ss.str(tmp);
		ss >> r.ysize >> r.xsize;

		Resize(r,r.ysize,r.xsize);
		r.lb1=1; r.ub1=r.ysize;
		r.lb2=1; r.ub2=r.xsize;
		

		//Read past empty lines
		while(tmp.size() == 0)
			std::getline(s,tmp);

		//Read Data
		int i=0, j=0;
                if(s_format == "array") {
			//Dense data
			while(s.good()) {
				std::getline(s,tmp);
				if(tmp.size() != 0) {
					mm_add_element(r.dat[i*r.xsize+j],tmp,field);
					i++;
					if(i==r.ysize) { 
						if(symmetry==mm_general) {
						  i=0; j++; 
						} else {
						  j++; i=j;
						}
					} 
				}
			}
		} else {
			//Sparse data
			r = 0.0;
			while(s.good()) {
				std::getline(s,tmp);
				if(tmp.size() != 0) {
					std::stringstream ss(tmp);
					ss >> i >> j >> tmp;
					i--; j--;
					mm_add_element(r.dat[i*r.xsize+j],tmp,field);
				}
			}
		}

		if(symmetry == mm_symmetric) {
			int i,j;
			for (i=0;i<r.ysize;i++) {
				for (j=0;j<i;j++)	{
					r.dat[j*r.xsize+i] = r.dat[i*r.xsize+j];
				}
			}
		} else if(symmetry == mm_skew_symmetric) {
			int i,j;
			for (i=0;i<r.ysize;i++) {
				for (j=0;j<i;j++)	{
					r.dat[j*r.xsize+i] = -r.dat[i*r.xsize+j];
				}
			}
		} else if(symmetry == mm_hermitian) {
			int i,j;
			for (i=0;i<r.ysize;i++) {
				for (j=0;j<i;j++)	{
					r.dat[j*r.xsize+i] = conj_el(r.dat[i*r.xsize+j]);
				}
			}
		}
		return s;

	} else {
		int i,j;
		for (i=0;i<r.ysize;i++)
		{
			for (j=0;j<r.xsize;j++)
			{
				s >> r.dat[i*r.xsize+j];
			}
		}
		return s;
	}
}

	template <class MS>
	std::ostream &_msout(std::ostream &s,const MS &r) throw()
	{
		int i,j;
		for (i=0;i<r.sysize;i++)
		{
			for (j=0;j<r.sxsize;j++)
			{
				s << r.dat[(i+r.offset1)*r.mxsize+j+r.offset2]<<" ";
			}
			s<<std::endl;
		}
		return s;
	}

	template <class MS>
	std::istream &_msin(std::istream &s,MS &r) throw()
	{
		int i,j;
		for (i=0;i<r.sysize;i++)
		{
			for (j=0;j<r.sxsize;j++)
			{
				s >> r.dat[(i+r.offset1)*r.mxsize+j+r.offset2];
			}
		}
		return s;
	}

	template <class M1,class M2,class S>
	TINLINE M1 &_mmassign(M1 &m1,const M2 &m,S ms) throw()
	{
		S *ndat=new S[m.xsize*m.ysize];
		for(int i=0;i<m.xsize*m.ysize;i++)
			ndat[i]=m.dat[i];
		delete [] m1.dat;
		m1.dat=ndat;

		m1.xsize=m.xsize;
		m1.ysize=m.ysize;
		m1.lb1=m.lb1; m1.ub1=m.ub1;
		m1.lb2=m.lb2; m1.ub2=m.ub2;
		return m1;
	}

	template <class M,class MS2,class S>
	TINLINE M &_mmsassign(M &m,const MS2 &ms) throw()
	{
		S *ndat=new S[ms.sxsize*ms.sysize];
		m.lb1=ms.start1;
		m.ub1=ms.end1;
		m.lb2=ms.start2;
		m.ub2=ms.end2;
		m.xsize=ms.sxsize;
		m.ysize=ms.sysize;
		for(int i=0;i<m.ysize;i++)
		{
			for(int j=0;j<m.xsize;j++)
				ndat[i*m.xsize+j]=ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2];
		}
		delete [] m.dat;
		m.dat=ndat;
		return m;
	}

	template <class MS,class M>
	TINLINE MS &_msmassign(MS &ms,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m.xsize!=ms.sxsize)||(m.ysize!=ms.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS>(nameof(ms)+" &"+nameof(ms)+"::operator =(const "+nameof(m)+" &)"));
#endif
		for(int i=0;i<ms.sysize;i++)
		{
			for(int j=0;j<ms.sxsize;j++)
				ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2]=m.dat[i*m.xsize+j];
		}
		return ms;
	}

	template <class MS1,class MS2>
	TINLINE MS1 &_msmsassign(MS1 &ms1,const MS2 &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((ms.sxsize!=ms1.sxsize)||(ms.sysize!=ms1.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS1>(nameof(ms1)+" &"+nameof(ms1)+"::operator =(const "+nameof(ms)+" &)"));
#endif
		for(int i=0;i<ms1.sysize;i++)
		{
			for(int j=0;j<ms1.sxsize;j++)
				ms1.dat[(i+ms1.offset1)*ms1.mxsize+j+ms1.offset2]=ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2];
		}
		return ms1;
	}

	template <class M,class S>
	TINLINE M &_msassign(M &m,const S &r) throw()
	{
		for(int i=0;i<m.ysize;i++)
		{
			for(int j=0;j<m.xsize;j++)
				m.dat[i*m.xsize+j]=r;
		}
		return m;
	}

	template <class MS,class S>
	TINLINE MS &_mssassign(MS &ms,const S &r) throw()
	{
		for(int i=0;i<ms.sysize;i++)
		{
			for(int j=0;j<ms.sxsize;j++)
				ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2]=r;
		}
		return ms;
	}

template <class V,class M,class S>
TINLINE V &_vmassign(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__TYPE_CAST_OF_THICK_OBJ<M>)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(m.xsize!=1&&m.ysize!=1) cxscthrow(ERROR__TYPE_CAST_OF_THICK_OBJ<M>(nameof(v)+" &"+nameof(v)+"::operator =(const "+nameof(m)+" &)"));
#endif
	if(m.ysize==1)
	{
		delete [] v.dat;
		v.dat=new S[m.xsize];
		for(int i=0;i<m.xsize;i++)
			v.dat[i]=m.dat[i];
		v.l=m.lb2;
		v.u=m.ub2;
		v.size=m.xsize;
	}
	else
	{
		delete [] v.dat;
		v.dat=new S[m.ysize];
		for(int i=0;i<m.ysize;i++)
			v.dat[i]=m.dat[i];
		v.l=m.lb1;
		v.u=m.ub1;
		v.size=m.ysize;
	}
	return v;
}

// result depends on the matrix (left arg), as defined in the reference to C-XSC
template <class M,class V,class S>
TINLINE M &_mvassign(M &m,const V &v) throw()
{
	delete [] m.dat;
	m.dat=new S[v.size];
	for(int i=0;i<v.size;i++)
		m.dat[i]=v.dat[i];
	if(m.ysize==1&&m.xsize!=1)
	{
		m.lb1=m.ub1=m.ysize=1;
		m.lb2=v.l;
		m.ub2=v.u;
		m.xsize=v.size;
	}
	else
	{
		m.lb1=v.l;
		m.ub1=v.u;
		m.ysize=v.size;
		m.lb2=m.ub2=m.xsize=1;
	}
	return m;
}

	template <class M>
	TINLINE int _mlb(const M &m, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_ROW_OR_COL<M>)
#else
	throw()
#endif
	{
		switch(i)
		{
			case 1:
				return m.lb1;
			case 2:
				return m.lb2;
			default:
#if(CXSC_INDEX_CHECK)
				throw ERROR__WRONG_ROW_OR_COL<M>("int Lb(const "+nameof(m)+" &, const int &)");
#endif
				return 0;
		}
	}
	
	template <class M>
	TINLINE int _mub(const M &m, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_ROW_OR_COL<M>)
#else
	throw()
#endif
	{
		switch(i)
		{
			case 1:
				return m.ub1;
			case 2:
				return m.ub2;
			default:
#if(CXSC_INDEX_CHECK)
				throw ERROR__WRONG_ROW_OR_COL<M>("int Ub(const "+nameof(m)+" &, const int &)");
#endif
				return 0;
		}
	}
	
	template <class MS>
	TINLINE int _mslb(const MS &ms, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_ROW_OR_COL<MS>)
#else
	throw()
#endif
	{
		switch(i)
		{
			case 1:
				return ms.start1;
			case 2:
				return ms.start2;
			default:
#if(CXSC_INDEX_CHECK)
				throw ERROR__WRONG_ROW_OR_COL<MS>("int Lb(const "+nameof(ms)+" &, const int &)");
#endif
				return 0;
		}
	}
	
	template <class MS>
	TINLINE int _msub(const MS &ms, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_ROW_OR_COL<MS>)
#else
	throw()
#endif
	{
		switch(i)
		{
			case 1:
				return ms.end1;
			case 2:
				return ms.end2;
			default:
#if(CXSC_INDEX_CHECK)
				throw ERROR__WRONG_ROW_OR_COL<MS>("int Ub(const "+nameof(ms)+" &, const int &)");
#endif
				return 0;
		}
	}
	
	template <class M>
	TINLINE M &_msetlb(M &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_ROW_OR_COL<M>)
#else
	throw()
#endif
	{
		switch(i)
		{
			case 1:
				m.lb1=j;
				m.ub1=j+m.ysize-1;
				break;
			case 2:
				m.lb2=j;
				m.ub2=j+m.xsize-1;
				break;
			default:
#if(CXSC_INDEX_CHECK)
				throw ERROR__WRONG_ROW_OR_COL<M>("voidnameof(m)+  &SetLb("+nameof(m)+" &, const int &,const int &)");
#endif
				break;
		}
		return m;
	}
	
	template <class M>
	TINLINE M &_msetub(M &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_ROW_OR_COL<M>)
#else
	throw()
#endif
	{
		switch(i)
		{
			case 1:
				m.ub1=j;
				m.lb1=j-m.ysize+1;
				break;
			case 2:
				m.ub2=j;
				m.lb2=j-m.xsize+1;
				break;
			default:
#if(CXSC_INDEX_CHECK)
				throw ERROR__WRONG_ROW_OR_COL<M>("voidnameof(m)+  &SetUb("+nameof(m)+" &, const int &,const int &)");
#endif
				break;
		}
		return m;
	}
	
	template <class M>
	TINLINE void _mresize(M &A) throw()
	{
		A.xsize=A.ysize=A.ub1=A.ub2=0;
		A.lb1=A.lb2=1;
		delete [] A.dat;
		A.dat=NULL;
	}

	template <class M,class S>
	TINLINE void _mresize(M &A,const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__WRONG_BOUNDARIES<M>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
	if((m<0)||(n<0)) cxscthrow(ERROR__WRONG_BOUNDARIES<M>("void Resize("+nameof(A)+" &,const int &, const int &)"));
#endif
		S *ndat=new S[m*n];
		int l1,u1,l2,u2;
		int i,j;
		l1=(A.lb1>1)?A.lb1:1;
		u1=(A.ub1<m)?A.ub1:m;
		l2=(A.lb2>1)?A.lb2:1;
		u2=(A.ub2<n)?A.ub2:n;
		for(i=0;i<l1-1;i++)
		{
			for(j=0;j<n;j++)
				ndat[i*n+j]=S(0);
		}
		for(;i<u1;i++)
		{
			for(j=0;j<l2-1;j++)
				ndat[i*n+j]=S(0);
			for(;j<u2;j++)
				ndat[i*n+j]=A.dat[(i-A.lb1+1)*A.xsize+j];
			for(;j<n;j++)
				ndat[i*n+j]=S(0);
		}
		for(;i<m;i++)
		{
			for(j=0;j<n;j++)
				ndat[i*n+j]=S(0);
		}
		delete [] A.dat;
		A.dat=ndat;
		A.xsize=A.ub2=n;
		A.ysize=A.ub1=m;
		A.lb1=A.lb2=1;
	}

	template <class M,class S>
	TINLINE void _mresize(M &A,const int &m1, const int &m2,const int &n1,const int &n2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__WRONG_BOUNDARIES<M>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
	if((m2<m1)||(n2<n1)) cxscthrow(ERROR__WRONG_BOUNDARIES<M>("void Resize("+nameof(A)+" &,const int &, const int &,const int &,const int &)"));
#endif
/*		int i,j,m,n;
		n=A.xsize=n2-n1+1;
		m=A.ysize=m2-m1+1;*/
		int m = m2-m1+1;
		int n = n2-n1+1;
		S *ndat=new S[m*n];

		//Initialize with zero
		for(int i=0 ; i<m ; i++)
			for(int j=0 ; j<n ; j++)
				ndat[i*n+j] = S(0);

		int new_l1, new_l2; //start and end indices in ndat
		int old_l1, old_u1, old_l2, old_u2; //start and end indices in A.dat

		//Find start and end indices for the rows
		if(A.lb1 > m1) {
		  old_l1 = 0;
		  new_l1 = A.lb1 - m1;
		} else if (A.lb1 == m1) {
		  old_l1 = 0;
		  new_l1 = 0;
		} else {
		  old_l1 = m1 - A.lb1;
		  new_l1 = 0;
		}

		if(A.ub1 < m2) {
		  old_u1 = A.ub1 - A.lb1 + 1;
		  //new_u1 = m2 - A.ub1 + 1;
		} else if (A.ub1 == m2) {
		  old_u1 = A.ub1 - A.lb1 + 1;
		  //new_u1 = m;
		} else {
		  old_u1 = A.ub1 - m2 + 1;
		  //new_u1 = m;
		}


		//Find start and end indices for the columns
		if(A.lb2 > n1) {
		  old_l2 = 0;
		  new_l2 = A.lb2 - n1;
		} else if (A.lb2 == n1) {
		  old_l2 = 0;
		  new_l2 = 0;
		} else {
		  old_l2 = n1 - A.lb2;
		  new_l2 = 0;
		}

		if(A.ub2 < n2) {
		  old_u2 = A.ub2 - A.lb2 + 1;
		  //new_u2 = n2 - A.ub2 + 1;
		} else if (A.ub2 == n2) {
		  old_u2 = A.ub2 - A.lb2 + 1;
		  //new_u2 = n;
		} else {
		  old_u2 = A.ub2 - n2 + 1;
		  //new_u2 = n;
		}


		//Copy values that are left after resize
		for(int i=old_l1 ; i<old_u1 ; i++)
			for(int j=old_l2 ; j<old_u2 ; j++)
				ndat[(i-old_l1+new_l1)*n+(j-old_l2+new_l2)] = A.dat[i*A.xsize+j];

		delete [] A.dat;

		A.dat = ndat;
		A.xsize = n;
		A.ysize = m;
		A.lb1 = m1;
		A.ub1 = m2;
		A.lb2 = n1;
		A.ub2 = n2;
	}
	
	template <class M,class E>
	TINLINE E _mabs(const M &m) throw()
	{
		E r(m.lb1,m.ub1,m.lb2,m.ub2);
		for(int i=0;i<m.xsize*m.ysize;i++)
			r.dat[i]=abs(m.dat[i]);
		return r;
	}

	template <class MS,class E>
	TINLINE E _msabs(const MS &ms) throw()
	{
		E r(ms.start1,ms.end1,ms.start2,ms.end2);
		for(int i=0;i<ms.sysize;i++)
		{
			for(int j=0;j<ms.sxsize;j++)
				r.dat[i*ms.sxsize+j]=abs(ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2]);
		}
		return r;
	}
	
	template <class M,class E>
	TINLINE E _mdiam(const M &m) throw()
	{
		E r(m.lb1,m.ub1,m.lb2,m.ub2);
		for(int i=0;i<m.xsize*m.ysize;i++)
			r.dat[i]=diam(m.dat[i]);
		return r;
	}

	template <class MS,class E>
	TINLINE E _msdiam(const MS &ms) throw()
	{
		E r(ms.start1,ms.end1,ms.start2,ms.end2);
		for(int i=0;i<ms.sysize;i++)
		{
			for(int j=0;j<ms.sxsize;j++)
				r.dat[i*ms.sxsize+j]=diam(ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2]);
		}
		return r;
	}
	
	template <class M,class E>
	TINLINE E _mmid(const M &m) throw()
	{
		E r(m.lb1,m.ub1,m.lb2,m.ub2);
		for(int i=0;i<m.xsize*m.ysize;i++)
			r.dat[i]=mid(m.dat[i]);
		return r;
	}

	template <class MS,class E>
	TINLINE E _msmid(const MS &ms) throw()
	{
		E r(ms.start1,ms.end1,ms.start2,ms.end2);
		for(int i=0;i<ms.sysize;i++)
		{
			for(int j=0;j<ms.sxsize;j++)
				r.dat[i*ms.sxsize+j]=mid(ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2]);
		}
		return r;
	}
	
	template <class M,class E>
	TINLINE E _minf(const M &m) throw()
	{
		E r(m.lb1,m.ub1,m.lb2,m.ub2);
		for(int i=0;i<m.xsize*m.ysize;i++)
			r.dat[i]=Inf(m.dat[i]);
		return r;
	}

	template <class MS,class E>
	TINLINE E _msinf(const MS &ms) throw()
	{
		E r(ms.start1,ms.end1,ms.start2,ms.end2);
		for(int i=0;i<ms.sysize;i++)
		{
			for(int j=0;j<ms.sxsize;j++)
				r.dat[i*ms.sxsize+j]=Inf(ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2]);
		}
		return r;
	}
	
	template <class M,class E>
	TINLINE E _msup(const M &m) throw()
	{
		E r(m.lb1,m.ub1,m.lb2,m.ub2);
		for(int i=0;i<m.xsize*m.ysize;i++)
			r.dat[i]=Sup(m.dat[i]);
		return r;
	}

	template <class MS,class E>
	TINLINE E _mssup(const MS &ms) throw()
	{
		E r(ms.start1,ms.end1,ms.start2,ms.end2);
		for(int i=0;i<ms.sysize;i++)
		{
			for(int j=0;j<ms.sxsize;j++)
				r.dat[i*ms.sxsize+j]=Sup(ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2]);
		}
		return r;
	}
	
	template <class M1,class M2>
	TINLINE M1 &_mmsetinf(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=m2.xsize)||(m1.ysize!=m2.ysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M1>(nameof(m1)+" &SetInf("+nameof(m1)+" &,const "+nameof(m2)+" &)"));
#endif
		for(int i=0;i<m1.xsize*m1.ysize;i++)
			SetInf(m1.dat[i],m2.dat[i]);
		return m1;
	}

	template <class MS1,class M2>
	TINLINE MS1 &_msmsetinf(MS1 &ms1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((ms1.sxsize!=m2.xsize)||(ms1.sysize!=m2.ysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS1>(nameof(ms1)+" &SetInf("+nameof(ms1)+" &,const "+nameof(m2)+" &)"));
#endif
		for(int i=0;i<ms1.sysize;i++)
		{
			for(int j=0;j<ms1.sxsize;j++)
				SetInf(ms1.dat[(i+ms1.offset1)*ms1.mxsize+j+ms1.offset2],m2.dat[i*m2.xsize+j]);
		}
		return ms1;
	}
	
	template <class M1,class MS2>
	TINLINE M1 &_mmssetinf(M1 &m1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=ms2.sxsize)||(m1.ysize!=ms2.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M1>(nameof(m1)+" &SetInf("+nameof(m1)+" &,const "+nameof(ms2)+" &)"));
#endif
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<m1.xsize;j++)
				SetInf(m1.dat[i*m1.xsize+j],ms2.dat[(i+ms2.offset1)*ms2.mxsize+j+ms2.offset2]);
		}
		return m1;
	}
	
	template <class MS1,class MS2>
	TINLINE MS1 &_msmssetinf(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((ms1.sxsize!=ms2.sxsize)||(ms1.sysize!=ms2.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS1>(nameof(ms1)+" &SetInf("+nameof(ms1)+" &,const "+nameof(ms2)+" &)"));
#endif
		for(int i=0;i<ms1.sysize;i++)
		{
			for(int j=0;j<ms1.sxsize;j++)
				SetInf(ms1.dat[(i+ms1.offset1)*ms1.mxsize+j+ms1.offset2],ms2.dat[(i+ms2.offset1)*ms2.mxsize+j+ms2.offset2]);
		}
		return ms1;
	}
	
	template <class M1,class M2>
	TINLINE M1 &_mmsetsup(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=m2.xsize)||(m1.ysize!=m2.ysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M1>(nameof(m1)+" &SetSup("+nameof(m1)+" &,const "+nameof(m2)+" &)"));
#endif
		for(int i=0;i<m1.xsize*m1.ysize;i++)
			SetSup(m1.dat[i],m2.dat[i]);
		return m1;
	}

	template <class MS1,class M2>
	TINLINE MS1 &_msmsetsup(MS1 &ms1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((ms1.sxsize!=m2.xsize)||(ms1.sysize!=m2.ysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS1>(nameof(ms1)+" &SetSup("+nameof(ms1)+" &,const "+nameof(m2)+" &)"));
#endif
		for(int i=0;i<ms1.sysize;i++)
		{
			for(int j=0;j<ms1.sxsize;j++)
				SetSup(ms1.dat[(i+ms1.offset1)*ms1.mxsize+j+ms1.offset2],m2.dat[i*m2.xsize+j]);
		}
		return ms1;
	}
	
	template <class M1,class MS2>
	TINLINE M1 &_mmssetsup(M1 &m1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=ms2.sxsize)||(m1.ysize!=ms2.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M1>(nameof(m1)+" &SetSup("+nameof(m1)+" &,const "+nameof(ms2)+" &)"));
#endif
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<m1.xsize;j++)
				SetSup(m1.dat[i*m1.xsize+j],ms2.dat[(i+ms2.offset1)*ms2.mxsize+j+ms2.offset2]);
		}
		return m1;
	}
	
	template <class MS1,class MS2>
	TINLINE MS1 &_msmssetsup(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((ms1.sxsize!=ms2.sxsize)||(ms1.sysize!=ms2.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS1>(nameof(ms1)+" &SetSup("+nameof(ms1)+" &,const "+nameof(ms2)+" &)"));
#endif
		for(int i=0;i<ms1.sysize;i++)
		{
			for(int j=0;j<ms1.sxsize;j++)
				SetSup(ms1.dat[(i+ms1.offset1)*ms1.mxsize+j+ms1.offset2],ms2.dat[(i+ms2.offset1)*ms2.mxsize+j+ms2.offset2]);
		}
		return ms1;
	}
	
	template <class M1,class M2>
	TINLINE M1 &_mmusetinf(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=m2.xsize)||(m1.ysize!=m2.ysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M1>(nameof(m1)+" &UncheckedSetInf("+nameof(m1)+" &,const "+nameof(m2)+" &)"));
#endif
		for(int i=0;i<m1.xsize*m1.ysize;i++)
			UncheckedSetInf(m1.dat[i],m2.dat[i]);
		return m1;
	}

	template <class MS1,class M2>
	TINLINE MS1 &_msmusetinf(MS1 &ms1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((ms1.sxsize!=m2.xsize)||(ms1.sysize!=m2.ysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS1>(nameof(ms1)+" &UncheckedSetInf("+nameof(ms1)+" &,const "+nameof(m2)+" &)"));
#endif
		for(int i=0;i<ms1.sysize;i++)
		{
			for(int j=0;j<ms1.sxsize;j++)
				UncheckedSetInf(ms1.dat[(i+ms1.offset1)*ms1.mxsize+j+ms1.offset2],m2.dat[i*m2.xsize+j]);
		}
		return ms1;
	}
	
	template <class M1,class MS2>
	TINLINE M1 &_mmsusetinf(M1 &m1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=ms2.sxsize)||(m1.ysize!=ms2.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M1>(nameof(m1)+" &UncheckedSetInf("+nameof(m1)+" &,const "+nameof(ms2)+" &)"));
#endif
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<m1.xsize;j++)
				UncheckedSetInf(m1.dat[i*m1.xsize+j],ms2.dat[(i+ms2.offset1)*ms2.mxsize+j+ms2.offset2]);
		}
		return m1;
	}
	
	template <class MS1,class MS2>
	TINLINE MS1 &_msmsusetinf(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((ms1.sxsize!=ms2.sxsize)||(ms1.sysize!=ms2.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS1>(nameof(ms1)+" &UncheckedSetInf("+nameof(ms1)+" &,const "+nameof(ms2)+" &)"));
#endif
		for(int i=0;i<ms1.sysize;i++)
		{
			for(int j=0;j<ms1.sxsize;j++)
				UncheckedSetInf(ms1.dat[(i+ms1.offset1)*ms1.mxsize+j+ms1.offset2],ms2.dat[(i+ms2.offset1)*ms2.mxsize+j+ms2.offset2]);
		}
		return ms1;
	}
	
	template <class M1,class M2>
	TINLINE M1 &_mmusetsup(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=m2.xsize)||(m1.ysize!=m2.ysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M1>(nameof(m1)+" &UncheckedSetSup("+nameof(m1)+" &,const "+nameof(m2)+" &)"));
#endif
		for(int i=0;i<m1.xsize*m1.ysize;i++)
			UncheckedSetSup(m1.dat[i],m2.dat[i]);
		return m1;
	}

	template <class MS1,class M2>
	TINLINE MS1 &_msmusetsup(MS1 &ms1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((ms1.sxsize!=m2.xsize)||(ms1.sysize!=m2.ysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS1>(nameof(ms1)+" &UncheckedSetSup("+nameof(ms1)+" &,const "+nameof(m2)+" &)"));
#endif
		for(int i=0;i<ms1.sysize;i++)
		{
			for(int j=0;j<ms1.sxsize;j++)
				UncheckedSetSup(ms1.dat[(i+ms1.offset1)*ms1.mxsize+j+ms1.offset2],m2.dat[i*m2.xsize+j]);
		}
		return ms1;
	}
	
	template <class M1,class MS2>
	TINLINE M1 &_mmsusetsup(M1 &m1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=ms2.sxsize)||(m1.ysize!=ms2.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M1>(nameof(m1)+" &UncheckedSetSup("+nameof(m1)+" &,const "+nameof(ms2)+" &)"));
#endif
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<m1.xsize;j++)
				UncheckedSetSup(m1.dat[i*m1.xsize+j],ms2.dat[(i+ms2.offset1)*ms2.mxsize+j+ms2.offset2]);
		}
		return m1;
	}
	
	template <class MS1,class MS2>
	TINLINE MS1 &_msmsusetsup(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((ms1.sxsize!=ms2.sxsize)||(ms1.sysize!=ms2.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS1>(nameof(ms1)+" &UncheckedSetSup("+nameof(ms1)+" &,const "+nameof(ms2)+" &)"));
#endif
		for(int i=0;i<ms1.sysize;i++)
		{
			for(int j=0;j<ms1.sxsize;j++)
				UncheckedSetSup(ms1.dat[(i+ms1.offset1)*ms1.mxsize+j+ms1.offset2],ms2.dat[(i+ms2.offset1)*ms2.mxsize+j+ms2.offset2]);
		}
		return ms1;
	}
	
	template <class M1,class M2>
	TINLINE M1 &_mmsetre(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=m2.xsize)||(m1.ysize!=m2.ysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M1>(nameof(m1)+" &SetRe("+nameof(m1)+" &,const "+nameof(m2)+" &)"));
#endif
		for(int i=0;i<m1.xsize*m1.ysize;i++)
			SetRe(m1.dat[i],m2.dat[i]);
		return m1;
	}

	template <class MS1,class M2>
	TINLINE MS1 &_msmsetre(MS1 &ms1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((ms1.sxsize!=m2.xsize)||(ms1.sysize!=m2.ysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS1>(nameof(ms1)+" &SetRe("+nameof(ms1)+" &,const "+nameof(m2)+" &)"));
#endif
		for(int i=0;i<ms1.sysize;i++)
		{
			for(int j=0;j<ms1.sxsize;j++)
				SetRe(ms1.dat[(i+ms1.offset1)*ms1.mxsize+j+ms1.offset2],m2.dat[i*m2.xsize+j]);
		}
		return ms1;
	}
	
	template <class M1,class MS2>
	TINLINE M1 &_mmssetre(M1 &m1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=ms2.sxsize)||(m1.ysize!=ms2.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M1>(nameof(m1)+" &SetRe("+nameof(m1)+" &,const "+nameof(ms2)+" &)"));
#endif
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<m1.xsize;j++)
				SetRe(m1.dat[i*m1.xsize+j],ms2.dat[(i+ms2.offset1)*ms2.mxsize+j+ms2.offset2]);
		}
		return m1;
	}
	
	template <class MS1,class MS2>
	TINLINE MS1 &_msmssetre(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((ms1.sxsize!=ms2.sxsize)||(ms1.sysize!=ms2.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS1>(nameof(ms1)+" &SetRe("+nameof(ms1)+" &,const "+nameof(ms2)+" &)"));
#endif
		for(int i=0;i<ms1.sysize;i++)
		{
			for(int j=0;j<ms1.sxsize;j++)
				SetRe(ms1.dat[(i+ms1.offset1)*ms1.mxsize+j+ms1.offset2],ms2.dat[(i+ms2.offset1)*ms2.mxsize+j+ms2.offset2]);
		}
		return ms1;
	}
	
	template <class M1,class M2>
	TINLINE M1 &_mmsetim(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=m2.xsize)||(m1.ysize!=m2.ysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M1>(nameof(m1)+" &SetIm("+nameof(m1)+" &,const "+nameof(m2)+" &)"));
#endif
		for(int i=0;i<m1.xsize*m1.ysize;i++)
			SetIm(m1.dat[i],m2.dat[i]);
		return m1;
	}

	template <class MS1,class M2>
	TINLINE MS1 &_msmsetim(MS1 &ms1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((ms1.sxsize!=m2.xsize)||(ms1.sysize!=m2.ysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS1>(nameof(ms1)+" &SetIm("+nameof(ms1)+" &,const "+nameof(m2)+" &)"));
#endif
		for(int i=0;i<ms1.sysize;i++)
		{
			for(int j=0;j<ms1.sxsize;j++)
				SetIm(ms1.dat[(i+ms1.offset1)*ms1.mxsize+j+ms1.offset2],m2.dat[i*m2.xsize+j]);
		}
		return ms1;
	}
	
	template <class M1,class MS2>
	TINLINE M1 &_mmssetim(M1 &m1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=ms2.sxsize)||(m1.ysize!=ms2.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M1>(nameof(m1)+" &SetIm("+nameof(m1)+" &,const "+nameof(ms2)+" &)"));
#endif
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<m1.xsize;j++)
				SetIm(m1.dat[i*m1.xsize+j],ms2.dat[(i+ms2.offset1)*ms2.mxsize+j+ms2.offset2]);
		}
		return m1;
	}
	
	template <class MS1,class MS2>
	TINLINE MS1 &_msmssetim(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((ms1.sxsize!=ms2.sxsize)||(ms1.sysize!=ms2.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS1>(nameof(ms1)+" &SetIm("+nameof(ms1)+" &,const "+nameof(ms2)+" &)"));
#endif
		for(int i=0;i<ms1.sysize;i++)
		{
			for(int j=0;j<ms1.sxsize;j++)
				SetIm(ms1.dat[(i+ms1.offset1)*ms1.mxsize+j+ms1.offset2],ms2.dat[(i+ms2.offset1)*ms2.mxsize+j+ms2.offset2]);
		}
		return ms1;
	}
	
	template <class M,class E>
	TINLINE E _mim(const M &m) throw()
	{
		E r(m.lb1,m.ub1,m.lb2,m.ub2);
		for(int i=0;i<m.xsize*m.ysize;i++)
			r.dat[i]=Im(m.dat[i]);
		return r;
	}

	template <class MS,class E>
	TINLINE E _msim(const MS &ms) throw()
	{
		E r(ms.start1,ms.end1,ms.start2,ms.end2);
		for(int i=0;i<ms.sysize;i++)
		{
			for(int j=0;j<ms.sxsize;j++)
				r.dat[i*ms.sxsize+j]=Im(ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2]);
		}
		return r;
	}
	
	template <class M,class E>
	TINLINE E _mre(const M &m) throw()
	{
		E r(m.lb1,m.ub1,m.lb2,m.ub2);
		for(int i=0;i<m.xsize*m.ysize;i++)
			r.dat[i]=Re(m.dat[i]);
		return r;
	}

	template <class MS,class E>
	TINLINE E _msre(const MS &ms) throw()
	{
		E r(ms.start1,ms.end1,ms.start2,ms.end2);
		for(int i=0;i<ms.sysize;i++)
		{
			for(int j=0;j<ms.sxsize;j++)
				r.dat[i*ms.sxsize+j]=Re(ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2]);
		}
		return r;
	}
	
	template <class M1,class M2,class E>
	TINLINE E _mmsect(const M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=m2.xsize)||(m1.ysize!=m2.ysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M1>(nameof(E())+" operator &(const "+nameof(m1)+" &,const "+nameof(m2)+" &)"));
#endif
		int x=1,y=1;

		if((m1.lb1==m2.lb1)&&(m1.lb2==m2.lb2)) { y=m1.lb1; x=m1.lb2; }
		
		E r(y,y+m1.ysize-1,x,x+m1.xsize-1);
		// assume that the matrices are of the same size, no bounds checking here
		for(int i=0;i<m1.xsize*m1.ysize;i++)
			r.dat[i]=m1.dat[i]&m2.dat[i];
		return r;
	}
	
	template <class M,class MS,class E>
	TINLINE E _mmssect(const M &m,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m.xsize!=ms.sxsize)||(m.ysize!=ms.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(E())+" operator &(const "+nameof(m)+" &,const "+nameof(ms)+" &)"));
#endif
		int x=1,y=1;

		if((m.lb1==ms.start1)&&(m.lb2==ms.start2)) { y=m.lb1; x=m.lb2; }
		
		E r(y,y+m.ysize-1,x,x+m.xsize-1);

		// assume that the matrices are of the same size, no bounds checking here
		for(int i=0;i<m.ysize;i++)
		{
			for(int j=0;j<m.xsize;j++)
				r.dat[i*m.xsize+j]=m.dat[i*m.xsize+j]&ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2];
		}

		return r;
	}

	template <class M1,class M2,class E>
	TINLINE E _mmconv(const M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=m2.xsize)||(m1.ysize!=m2.ysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M1>(nameof(E())+" operator |(const "+nameof(m1)+" &,const "+nameof(m2)+" &)"));
#endif
		int x=1,y=1;

		if((m1.lb1==m2.lb1)&&(m1.lb2==m2.lb2)) { y=m1.lb1; x=m1.lb2; }
		
		E r(y,y+m1.ysize-1,x,x+m1.xsize-1);
		// assume that the matrices are of the same size, no bounds checking here
		for(int i=0;i<m1.xsize*m1.ysize;i++)
			r.dat[i]=m1.dat[i]|m2.dat[i];
		return r;
	}
	
	template <class M,class MS,class E>
	TINLINE E _mmsconv(const M &m,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m.xsize!=ms.sxsize)||(m.ysize!=ms.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(E())+" operator |(const "+nameof(m)+" &,const "+nameof(ms)+" &)"));
#endif
		int x=1,y=1;

		if((m.lb1==ms.start1)&&(m.lb2==ms.start2)) { y=m.lb1; x=m.lb2; }
		
		E r(y,y+m.ysize-1,x,x+m.xsize-1);

		// assume that the matrices are of the same size, no bounds checking here
		for(int i=0;i<m.ysize;i++)
		{
			for(int j=0;j<m.xsize;j++)
				r.dat[i*m.xsize+j]=m.dat[i*m.xsize+j]|ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2];
		}

		return r;
	}

	template <class M1,class M2,class E>
	TINLINE E _mmplus(const M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=m2.xsize)||(m1.ysize!=m2.ysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M1>(nameof(E())+" operator +(const "+nameof(m1)+" &,const "+nameof(m2)+" &)"));
#endif
		int x=1,y=1;

		if((m1.lb1==m2.lb1)&&(m1.lb2==m2.lb2)) { y=m1.lb1; x=m1.lb2; }
		
		E r(y,y+m1.ysize-1,x,x+m1.xsize-1);
		// assume that the matrices are of the same size, no bounds checking here
		for(int i=0;i<m1.xsize*m1.ysize;i++)
			r.dat[i]=m1.dat[i]+m2.dat[i];
		return r;
	}
	
	template <class M,class MS,class E>
	TINLINE E _mmsplus(const M &m,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m.xsize!=ms.sxsize)||(m.ysize!=ms.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(E())+" operator +(const "+nameof(m)+" &,const "+nameof(ms)+" &)"));
#endif
		int x=1,y=1;

		if((m.lb1==ms.start1)&&(m.lb2==ms.start2)) { y=m.lb1; x=m.lb2; }
		
		E r(y,y+m.ysize-1,x,x+m.xsize-1);

		// assume that the matrices are of the same size, no bounds checking here
		for(int i=0;i<m.ysize;i++)
		{
			for(int j=0;j<m.xsize;j++)
				r.dat[i*m.xsize+j]=m.dat[i*m.xsize+j]+ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2];
		}

		return r;
	}

	template <class M>
	TINLINE M _mminus(const M &m) throw()
	{
		M sum(m.lb1,m.ub1,m.lb2,m.ub2);
		for(int i=0;i<m.xsize*m.ysize;i++)
			sum.dat[i]= -m.dat[i];
		return sum;
	}
	
	template <class MS,class E>
	TINLINE E _msminus(const MS &ms) throw()
	{
		E r(ms.start1,ms.end1,ms.start2,ms.end2);
		for(int i=0;i<ms.sysize;i++)
		{
			for(int j=0;j<ms.sxsize;j++)
				r.dat[i*ms.sxsize+j]= -ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2];
		}
		return r;
	}
	
	template <class M1,class M2,class E>
	TINLINE E _mmminus(const M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=m2.xsize)||(m1.ysize!=m2.ysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M1>(nameof(E())+" operator -(const "+nameof(m1)+" &,const "+nameof(m2)+" &)"));
#endif
		int x=1,y=1;
		if((m1.lb1==m2.lb1)&&(m1.lb2==m2.lb2)) { y=m1.lb1; x=m1.lb2; }

		E r(y,y+m1.ysize-1,x,x+m1.xsize-1);
		// assume that the matrices are of the same size, no bounds checking here
		for(int i=0;i<m1.xsize*m1.ysize;i++)
			r.dat[i]=m1.dat[i]-m2.dat[i];

		return r;
	}


	template <class M1,class M2>
	TINLINE M1 &_mmconvassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=m2.xsize)||(m1.ysize!=m2.ysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M1>(nameof(m1)+" &operator |=("+nameof(m1)+" &,const "+nameof(m2)+" &)"));
#endif
		for(int i=0;i<m1.ysize*m1.xsize;i++)
			m1.dat[i]|=m2.dat[i];
		return m1;
	}
		
	template <class M,class MS>
	TINLINE M &_mmsconvassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=ms.sxsize)||(m1.ysize!=ms.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(m1)+" &operator |=("+nameof(m1)+" &,const "+nameof(ms)+" &)"));
#endif
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<m1.xsize;j++)
				m1.dat[i*m1.xsize+j]|=ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2];
		}
		return m1;
	}
		
	template <class MS,class M>
	TINLINE MS &_msmconvassign(MS &ms,const M &m1)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=ms.sxsize)||(m1.ysize!=ms.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS>(nameof(ms)+" &"+nameof(ms)+"::operator |=(const "+nameof(m1)+" &)"));
#endif
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<m1.xsize;j++)
				ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2]|=m1.dat[i*m1.xsize+j];
		}
		return ms;
	}
		
	template <class MS1,class MS2>
	TINLINE MS1 &_msmsconvassign(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((ms1.sxsize!=ms2.sxsize)||(ms1.sysize!=ms2.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS1>(nameof(ms1)+" &"+nameof(ms1)+"::operator |=(const "+nameof(ms2)+" &)"));
#endif
		for(int i=0;i<ms1.sysize;i++)
		{
			for(int j=0;j<ms1.sxsize;j++)
				ms1.dat[(i+ms1.offset1)*ms1.mxsize+j+ms1.offset2]|=ms2.dat[(i+ms2.offset1)*ms2.mxsize+j+ms2.offset2];
		}
		return ms1;
	}
		
	template <class M1,class M2>
	TINLINE M1 &_mmsectassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=m2.xsize)||(m1.ysize!=m2.ysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M1>(nameof(m1)+" &operator &=("+nameof(m1)+" &,const "+nameof(m2)+" &)"));
#endif
		for(int i=0;i<m1.ysize*m1.xsize;i++)
			m1.dat[i]&=m2.dat[i];
		return m1;
	}
		
	template <class M,class MS>
	TINLINE M &_mmssectassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=ms.sxsize)||(m1.ysize!=ms.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(m1)+" &operator &=("+nameof(m1)+" &,const "+nameof(ms)+" &)"));
#endif
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<m1.xsize;j++)
				m1.dat[i*m1.xsize+j]&=ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2];
		}
		return m1;
	}
		
	template <class MS,class M>
	TINLINE MS &_msmsectassign(MS &ms,const M &m1)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=ms.sxsize)||(m1.ysize!=ms.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS>(nameof(ms)+" &"+nameof(ms)+"::operator &=(const "+nameof(m1)+" &)"));
#endif
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<m1.xsize;j++)
				ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2]&=m1.dat[i*m1.xsize+j];
		}
		return ms;
	}
		
	template <class MS1,class MS2>
	TINLINE MS1 &_msmssectassign(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((ms1.sxsize!=ms2.sxsize)||(ms1.sysize!=ms2.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS1>(nameof(ms1)+" &"+nameof(ms1)+"::operator &=(const "+nameof(ms2)+" &)"));
#endif
		for(int i=0;i<ms1.sysize;i++)
		{
			for(int j=0;j<ms1.sxsize;j++)
				ms1.dat[(i+ms1.offset1)*ms1.mxsize+j+ms1.offset2]&=ms2.dat[(i+ms2.offset1)*ms2.mxsize+j+ms2.offset2];
		}
		return ms1;
	}
		
	template <class M1,class M2>
	TINLINE M1 &_mmplusassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=m2.xsize)||(m1.ysize!=m2.ysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M1>(nameof(m1)+" &operator +=("+nameof(m1)+" &,const "+nameof(m2)+" &)"));
#endif
		for(int i=0;i<m1.ysize*m1.xsize;i++)
			m1.dat[i]+=m2.dat[i];
		return m1;
	}
		
	template <class M,class MS>
	TINLINE M &_mmsplusassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=ms.sxsize)||(m1.ysize!=ms.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(m1)+" &operator +=("+nameof(m1)+" &,const "+nameof(ms)+" &)"));
#endif
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<m1.xsize;j++)
				m1.dat[i*m1.xsize+j]+=ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2];
		}
		return m1;
	}
		
	template <class MS,class M>
	TINLINE MS &_msmplusassign(MS &ms,const M &m1)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=ms.sxsize)||(m1.ysize!=ms.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS>(nameof(ms)+" &"+nameof(ms)+"::operator +=(const "+nameof(m1)+" &)"));
#endif
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<m1.xsize;j++)
				ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2]+=m1.dat[i*m1.xsize+j];
		}
		return ms;
	}
		
	template <class MS1,class MS2>
	TINLINE MS1 &_msmsplusassign(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((ms1.sxsize!=ms2.sxsize)||(ms1.sysize!=ms2.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS1>(nameof(ms1)+" &"+nameof(ms1)+"::operator +=(const "+nameof(ms2)+" &)"));
#endif
		for(int i=0;i<ms1.sysize;i++)
		{
			for(int j=0;j<ms1.sxsize;j++)
				ms1.dat[(i+ms1.offset1)*ms1.mxsize+j+ms1.offset2]+=ms2.dat[(i+ms2.offset1)*ms2.mxsize+j+ms2.offset2];
		}
		return ms1;
	}
		
	template <class M,class MS,class E>
	TINLINE E _mmsminus(const M &m,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m.xsize!=ms.sxsize)||(m.ysize!=ms.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(E())+" operator -(const "+nameof(m)+" &,const "+nameof(ms)+" &)"));
#endif
		int x=1,y=1;

		if((m.lb1==ms.start1)&&(m.lb2==ms.start2)) { y=m.lb1; x=m.lb2; }
		
		E r(y,y+m.ysize-1,x,x+m.xsize-1);

		// assume that the matrices are of the same size, no bounds checking here
		for(int i=0;i<m.ysize;i++)
		{
			for(int j=0;j<m.xsize;j++)
				r.dat[i*m.xsize+j]=m.dat[i*m.xsize+j]-ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2];
		}

		return r;
	}
	
	template <class MS,class M,class E>
	TINLINE E _msmminus(const MS &ms,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m.xsize!=ms.sxsize)||(m.ysize!=ms.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(E())+" operator -(const "+nameof(ms)+" &,const "+nameof(m)+" &)"));
#endif
		int x=1,y=1;

		if((m.lb1==ms.start1)&&(m.lb2==ms.start2)) { y=m.lb1; x=m.lb2; }
		
		E r(y,y+m.ysize-1,x,x+m.xsize-1);

		// assume that the matrices are of the same size, no bounds checking here
		for(int i=0;i<m.ysize;i++)
		{
			for(int j=0;j<m.xsize;j++)
				r.dat[i*m.xsize+j]=ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2]-m.dat[i*m.xsize+j];
		}
		return r;
	}
	
	template <class MS1,class MS2,class E>
	TINLINE E _msmsminus(const MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((ms1.sxsize!=ms2.sxsize)||(ms1.sysize!=ms2.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(E())+" operator -(const "+nameof(ms1)+" &,const "+nameof(ms2)+" &)"));
#endif
		int x=1,y=1;
		if((ms1.start1==ms2.start1)&&(ms1.start2==ms2.start2)) { y=ms1.start1; x=ms1.start2; }
		E r(y,y+ms1.sysize-1,x,x+ms1.sxsize-1);

		// assume that the matrices are of the same size, no bounds checking here
		for(int i=0;i<ms1.sysize;i++)
		{
			for(int j=0;j<ms1.sxsize;j++)
				r.dat[i*ms1.sxsize+j]=ms1.dat[(i+ms1.offset1)*ms1.mxsize+j+ms1.offset2]-ms2.dat[(i+ms2.offset1)*ms2.mxsize+j+ms2.offset2];
		}
		return r;
	}
	
	template <class M1,class M2>
	TINLINE M1 &_mmminusassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=m2.xsize)||(m1.ysize!=m2.ysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M1>(nameof(m1)+" &operator -=("+nameof(m1)+" &,const "+nameof(m2)+" &)"));
#endif
		for(int i=0;i<m1.ysize*m1.xsize;i++)
			m1.dat[i]-=m2.dat[i];
		return m1;
	}
		
	template <class M,class MS>
	TINLINE M &_mmsminusassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=ms.sxsize)||(m1.ysize!=ms.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(m1)+" &operator -=("+nameof(m1)+" &,const "+nameof(ms)+" &)"));
#endif
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<m1.xsize;j++)
				m1.dat[i*m1.xsize+j]-=ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2];
		}
		return m1;
	}
		
	template <class MS,class M>
	TINLINE MS &_msmminusassign(MS &ms,const M &m1)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.xsize!=ms.sxsize)||(m1.ysize!=ms.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS>(nameof(ms)+" &"+nameof(ms)+"::operator -=(const "+nameof(m1)+" &)"));
#endif
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<m1.xsize;j++)
				ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2]-=m1.dat[i*m1.xsize+j];
		}
		return ms;
	}
		
	template <class MS1,class MS2>
	TINLINE MS1 &_msmsminusassign(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((ms1.sxsize!=ms2.sxsize)||(ms1.sysize!=ms2.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS1>(nameof(ms1)+" &"+nameof(ms1)+"::operator -=(const "+nameof(ms2)+" &)"));
#endif
		for(int i=0;i<ms1.sysize;i++)
		{
			for(int j=0;j<ms1.sxsize;j++)
				ms1.dat[(i+ms1.offset1)*ms1.mxsize+j+ms1.offset2]-=ms2.dat[(i+ms2.offset1)*ms2.mxsize+j+ms2.offset2];
		}
		return ms1;
	}
		
	template <class M1,class M2,class E>
	TINLINE E _mmmult(const M1 &m1, const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E r(m1.lb1,m1.ub1,m2.lb2,m2.ub2);
#if(CXSC_INDEX_CHECK)
		if(m1.xsize!=m2.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(r)+" operator *(const "+nameof(m1)+" &, const "+nameof(m2)+" &)"));
#endif

		if(opdotprec == 1) {
#ifdef CXSC_USE_BLAS
			blasmatmul(m1,m2,r);
			return r;
#else
			r = 0.0;
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for
#endif
			for(int i=0 ; i<m1.ysize ; i++)
			  for(int k=0 ; k<m1.xsize ; k++)
			    for(int j=0 ; j<m2.xsize ; j++)
			      r.dat[i*m2.xsize+j] += m1.dat[i*m1.xsize+k] * m2.dat[k*m2.xsize+j];
			return r;
#endif
		} else {
			dotprecision dot(0.0);
			dot.set_k(opdotprec);
	
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(dot)
#endif
			for(int i=0;i<m1.ysize;i++)
			{
				for(int j=0;j<m2.xsize;j++)
				{
					dot=0.0;
					accumulate_approx(dot,m1[i+Lb(m1,ROW)],m2[Col(j+Lb(m2,COL))]);
					r.dat[i*m2.xsize+j]=rnd(dot);
				}
			}
	
			return r;
		}
	}

	template <class M1,class M2,class E>
	TINLINE E _mmimult(const M1 &m1, const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E r(m1.lb1,m1.ub1,m2.lb2,m2.ub2);
#if(CXSC_INDEX_CHECK)
		if(m1.xsize!=m2.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(r)+" operator *(const "+nameof(m1)+" &, const "+nameof(m2)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) {
			blasmatmul(m1,m2,r);
			return r;
		}
		else
#endif
		{
			idotprecision idot(0.0);
			idot.set_k(opdotprec);
	
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(idot)
#endif
			for(int i=0;i<m1.ysize;i++)
			{
				for(int j=0;j<m2.xsize;j++)
				{
					idot=0.0;
					accumulate(idot,m1[i+Lb(m1,ROW)],m2[Col(j+Lb(m2,COL))]);
					r.dat[i*m2.xsize+j]=rnd(idot);
				}
			}
	
			return r;
		}
	}

	template <class M1,class M2,class E>
	TINLINE E _mmlmult(const M1 &m1, const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E r(m1.lb1,m1.ub1,m2.lb2,m2.ub2);
#if(CXSC_INDEX_CHECK)
		if(m1.xsize!=m2.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(r)+" operator *(const "+nameof(m1)+" &, const "+nameof(m2)+" &)"));
#endif
		dotprecision dot(0.0);

#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(dot)
#endif
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<m2.xsize;j++)
			{
				dot=0.0;
				accumulate(dot,m1[i+Lb(m1,ROW)],m2[Col(j+Lb(m2,COL))]);
				r.dat[i*m2.xsize+j]=l_real(dot);
			}
		}

		return r;
	}

	template <class M1,class M2,class E>
	TINLINE E _mmlimult(const M1 &m1, const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E r(m1.lb1,m1.ub1,m2.lb2,m2.ub2);
#if(CXSC_INDEX_CHECK)
		if(m1.xsize!=m2.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(r)+" operator *(const "+nameof(m1)+" &, const "+nameof(m2)+" &)"));
#endif
		idotprecision idot(0.0);

#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(idot)
#endif
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<m2.xsize;j++)
			{
				idot=0.0;				
				accumulate(idot,m1[i+Lb(m1,ROW)],m2[Col(j+Lb(m2,COL))]);
				r.dat[i*m2.xsize+j]=l_interval(idot);
			}
		}

		return r;
	}

	template<class Tx, class Ty> inline complex fp_c_mult(const Tx& x, const Ty& y) {
          return x*y;
        }

	template<> inline complex fp_c_mult<complex,complex>(const complex& x, const complex& y) {
          return complex(Re(x)*Re(y)-Im(x)*Im(y), Re(x)*Im(y)+Im(x)*Re(y));
        }
        

	template <class M1,class M2,class E>
	TINLINE E _mmcmult(const M1 &m1, const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E r(m1.lb1,m1.ub1,m2.lb2,m2.ub2);
#if(CXSC_INDEX_CHECK)
		if(m1.xsize!=m2.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(r)+" operator *(const "+nameof(m1)+" &, const "+nameof(m2)+" &)"));
#endif

		if(opdotprec == 1) {
#ifdef CXSC_USE_BLAS
			blasmatmul(m1,m2,r);
			return r;
#else
			r = 0.0;
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for
#endif			
			for(int i=0 ; i<m1.ysize ; i++)
			  for(int k=0 ; k<m1.xsize ; k++)
			    for(int j=0 ; j<m2.xsize ; j++)
			      r.dat[i*m2.xsize+j] += fp_c_mult(m1.dat[i*m1.xsize+k], m2.dat[k*m2.xsize+j]);
			return r;
#endif
		} else {
			cdotprecision cdot(0.0);
			cdot.set_k(opdotprec);
	
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(cdot)
#endif
			for(int i=0;i<m1.ysize;i++)
			{
				for(int j=0;j<m2.xsize;j++)
				{
					cdot=0.0;
					accumulate_approx(cdot,m1[i+Lb(m1,ROW)],m2[Col(j+Lb(m2,COL))]);
					r.dat[i*m2.xsize+j]=rnd(cdot);
				}
			}
	
			return r;
		}
	}

	template <class M1,class M2,class E>
	TINLINE E _mmcimult(const M1 &m1, const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E r(m1.lb1,m1.ub1,m2.lb2,m2.ub2);
#if(CXSC_INDEX_CHECK)
		if(m1.xsize!=m2.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(r)+" operator *(const "+nameof(m1)+" &, const "+nameof(m2)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) {
			blasmatmul(m1,m2,r);
			return r;
		}
		else
#endif
		{
			cidotprecision cidot(0.0);
	
			cidot.set_k(opdotprec);
		
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(cidot)
#endif
			for(int i=0;i<m1.ysize;i++)
			{
				for(int j=0;j<m2.xsize;j++)
				{
					cidot=0.0;
					accumulate(cidot,m1[i+Lb(m1,ROW)],m2[Col(j+Lb(m2,COL))]);
					r.dat[i*m2.xsize+j]=rnd(cidot);
				}
			}
	
			return r;
		}
	}

	template <class M1,class M2,class S>
	TINLINE M1 &_mmmultassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(m1.xsize!=m2.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M1>(nameof(m1)+" &operator *=("+nameof(m1)+" &,const "+nameof(m2)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) 
		{
			M1 tmp(1,ColLen(m1),1,RowLen(m1));
			blasmatmul(m1,m2,tmp);
			m1 = tmp;
			return m1;
		}
		else
#endif
		{
			S *ndat=new S[m1.ysize*m2.xsize];
			dotprecision dot(0.0);
			dot.set_k(opdotprec);
	
			int lbr = Lb(m1,ROW);
			int lbc = Lb(m2,COL);
	
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(dot)
#endif
			for(int i=0;i<m1.ysize;i++)
			{
				for(int j=0;j<m2.xsize;j++)
				{
					dot=0.0;
					accumulate_approx(dot,m1[Row(lbr+i)],m2[Col(lbc+j)]);
					ndat[i*m2.xsize+j]=rnd(dot);
				}
			}
			delete [] m1.dat;
			m1.dat=ndat;
			m1.xsize=m2.xsize;
			m1.lb2=m2.lb2;
			m1.ub2=m2.ub2;
			return m1;
		}
	}

	template <class M1,class M2,class S>
	TINLINE M1 &_mmimultassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(m1.xsize!=m2.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M1>(nameof(m1)+" &operator *=("+nameof(m1)+" &,const "+nameof(m2)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) 
		{
			M1 tmp(1,ColLen(m1),1,RowLen(m1));
			blasmatmul(m1,m2,tmp);
			m1 = tmp;
			return m1;
		}
		else
#endif
		{
			S *ndat=new S[m1.ysize*m2.xsize];
			idotprecision idot(0.0);
			idot.set_k(opdotprec);
	
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(idot)
#endif
			for(int i=0;i<m1.ysize;i++)
			{
				for(int j=0;j<m2.xsize;j++)
				{
					idot=0.0;
					accumulate(idot,m1[i+Lb(m1,ROW)],m2[Col(j+Lb(m2,COL))]);
					ndat[i*m2.xsize+j]=rnd(idot);
				}
			}
			delete [] m1.dat;
			m1.dat=ndat;
			m1.xsize=m2.xsize;
			m1.lb2=m2.lb2;
			m1.ub2=m2.ub2;
			return m1;
		}
	}

	template <class M1,class M2,class S>
	TINLINE M1 &_mmlmultassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(m1.xsize!=m2.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M1>(nameof(m1)+" &operator *=("+nameof(m1)+" &,const "+nameof(m2)+" &)"));
#endif
		S *ndat=new S[m1.ysize*m2.xsize];
		dotprecision dot(0.0);

                int lbr = Lb(m1,ROW);
                int lbc = Lb(m2,COL);

#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(dot)
#endif
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<m2.xsize;j++)
			{
				dot=0.0;
                                accumulate(dot,m1[Row(lbr+i)],m2[Col(lbc+j)]);				ndat[i*m2.xsize+j]=l_real(dot);
			}
		}
		delete [] m1.dat;
		m1.dat=ndat;
		m1.xsize=m2.xsize;
		m1.lb2=m2.lb2;
		m1.ub2=m2.ub2;
		return m1;
	}

	template <class M1,class M2,class S>
	TINLINE M1 &_mmlimultassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(m1.xsize!=m2.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M1>(nameof(m1)+" &operator *=("+nameof(m1)+" &,const "+nameof(m2)+" &)"));
#endif
		S *ndat=new S[m1.ysize*m2.xsize];
		idotprecision idot(0.0);

#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(idot)
#endif
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<m2.xsize;j++)
			{
				idot=0.0;
				for(int k=0;k<m1.xsize;k++)
					accumulate(idot,m1.dat[i*m1.xsize+k],m2.dat[k*m2.xsize+j]);
				ndat[i*m2.xsize+j]=l_interval(idot);
			}
		}
		delete [] m1.dat;
		m1.dat=ndat;
		m1.xsize=m2.xsize;
		m1.lb2=m2.lb2;
		m1.ub2=m2.ub2;
		return m1;
	}

	template <class M1,class M2,class S>
	TINLINE M1 &_mmcmultassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(m1.xsize!=m2.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M1>(nameof(m1)+" &operator *=("+nameof(m1)+" &,const "+nameof(m2)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) 
		{
			M1 tmp(1,ColLen(m1),1,RowLen(m1));
			blasmatmul(m1,m2,tmp);
			m1 = tmp;
			return m1;
		}
		else
#endif
		{
			S *ndat=new S[m1.ysize*m2.xsize];
			cdotprecision cdot(0.0);
			cdot.set_k(opdotprec);
	
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(cdot)
#endif
			for(int i=0;i<m1.ysize;i++)
			{
				for(int j=0;j<m2.xsize;j++)
				{
					cdot=0.0;
					accumulate_approx(cdot,m1[i+Lb(m1,ROW)],m2[Col(j+Lb(m2,COL))]);
					ndat[i*m2.xsize+j]=rnd(cdot);
				}
			}
			delete [] m1.dat;
			m1.dat=ndat;
			m1.xsize=m2.xsize;
			m1.lb2=m2.lb2;
			m1.ub2=m2.ub2;
			return m1;
		}
	}

	template <class M1,class M2,class S>
	TINLINE M1 &_mmcimultassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(m1.xsize!=m2.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M1>(nameof(m1)+" &operator *=("+nameof(m1)+" &,const "+nameof(m2)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) 
		{
			M1 tmp(1,ColLen(m1),1,RowLen(m1));
			blasmatmul(m1,m2,tmp);
			m1 = tmp;
			return m1;
		}
		else
#endif
		{
			S *ndat=new S[m1.ysize*m2.xsize];
			cidotprecision cidot(0.0);
			cidot.set_k(opdotprec);
	
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(cidot)
#endif
			for(int i=0;i<m1.ysize;i++)
			{
				for(int j=0;j<m2.xsize;j++)
				{
					cidot=0.0;
					accumulate(cidot,m1[i+Lb(m1,ROW)],m2[Col(j+Lb(m2,COL))]);
					ndat[i*m2.xsize+j]=rnd(cidot);
				}
			}
			delete [] m1.dat;
			m1.dat=ndat;
			m1.xsize=m2.xsize;
			m1.lb2=m2.lb2;
			m1.ub2=m2.ub2;
			return m1;
		}
	}

	template <class M,class MS,class E>
	TINLINE E _mmsmult(const M &m1, const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E r(m1.lb1,m1.ub1,ms.start2,ms.end2);
#if(CXSC_INDEX_CHECK)
		if(m1.xsize!=ms.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(r)+" operator *(const "+nameof(m1)+" &, const "+nameof(ms)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) {
			blasmatmul(m1,ms,r);
			return r;
		}
		else
#endif
		{
			dotprecision dot(0.0);
			dot.set_k(opdotprec);
	
			int lbr = Lb(m1,ROW);
			int lbc = Lb(ms,COL);
	
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(dot)
#endif
			for(int i=0;i<m1.ysize;i++)
			{
				for(int j=0;j<ms.sxsize;j++)
				{
					dot=0.0;
					accumulate_approx(dot,m1[Row(lbr+i)],ms[Col(lbc+j)]);				r.dat[i*ms.sxsize+j]=rnd(dot);
				}
			}
			return r;
		}
	}

	template <class M,class MS,class E>
	TINLINE E _mmsimult(const M &m1, const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E r(m1.lb1,m1.ub1,ms.start2,ms.end2);
#if(CXSC_INDEX_CHECK)
		if(m1.xsize!=ms.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(r)+" operator *(const "+nameof(m1)+" &, const "+nameof(ms)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) {
			blasmatmul(m1,ms,r);
			return r;
		}
		else
#endif
		{
			idotprecision idot(0.0);
			idot.set_k(opdotprec);
	
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(idot)
#endif
			for(int i=0;i<m1.ysize;i++)
			{
				for(int j=0;j<ms.sxsize;j++)
				{
					idot=0.0;
					accumulate(idot,m1[i+Lb(m1,ROW)],ms[Col(j+Lb(ms,COL))]);
					r.dat[i*ms.sxsize+j]=rnd(idot);
				}
			}
			return r;
		}
	}

	template <class M,class MS,class E>
	TINLINE E _mmslmult(const M &m1, const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E r(m1.lb1,m1.ub1,ms.start2,ms.end2);
#if(CXSC_INDEX_CHECK)
		if(m1.xsize!=ms.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(r)+" operator *(const "+nameof(m1)+" &, const "+nameof(ms)+" &)"));
#endif
		dotprecision dot(0.0);

                int lbr = Lb(m1,ROW);
                int lbc = Lb(ms,COL);

#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(dot)
#endif
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<ms.sxsize;j++)
			{
				dot=0.0;
                                accumulate(dot,m1[Row(lbr+i)],ms[Col(lbc+j)]);			
				r.dat[i*ms.sxsize+j]=l_real(dot);
			}
		}
		return r;
	}

	template <class M,class MS,class E>
	TINLINE E _mmslimult(const M &m1, const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E r(m1.lb1,m1.ub1,ms.start2,ms.end2);
#if(CXSC_INDEX_CHECK)
		if(m1.xsize!=ms.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(r)+" operator *(const "+nameof(m1)+" &, const "+nameof(ms)+" &)"));
#endif
		idotprecision idot(0.0);

#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(idot)
#endif
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<ms.sxsize;j++)
			{
				idot=0.0;
				for(int k=0;k<m1.xsize;k++)
					accumulate(idot,m1.dat[i*m1.xsize+k],ms.dat[(k+ms.offset1)*ms.mxsize+j+ms.offset2]);
				r.dat[i*ms.sxsize+j]=l_interval(idot);
			}
		}
		return r;
	}

	template <class M,class MS,class E>
	TINLINE E _mmscmult(const M &m1, const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E r(m1.lb1,m1.ub1,ms.start2,ms.end2);
#if(CXSC_INDEX_CHECK)
		if(m1.xsize!=ms.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(r)+" operator *(const "+nameof(m1)+" &, const "+nameof(ms)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) {
			blasmatmul(m1,ms,r);
			return r;
		}
		else
#endif
		{			cdotprecision cdot(0.0);
			cdot.set_k(opdotprec);
	
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(cdot)
#endif
			for(int i=0;i<m1.ysize;i++)
			{
				for(int j=0;j<ms.sxsize;j++)
				{
					cdot=0.0;
					accumulate_approx(cdot,m1[i+Lb(m1,ROW)],ms[Col(j+Lb(ms,COL))]);
					r.dat[i*ms.sxsize+j]=rnd(cdot);
				}
			}
			return r;
		}
	}

	template <class M,class MS,class E>
	TINLINE E _mmscimult(const M &m1, const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E r(m1.lb1,m1.ub1,ms.start2,ms.end2);
#if(CXSC_INDEX_CHECK)
		if(m1.xsize!=ms.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(r)+" operator *(const "+nameof(m1)+" &, const "+nameof(ms)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) {
			blasmatmul(m1,ms,r);
			return r;
		}
		else
#endif
		{
			cidotprecision cidot(0.0);
			cidot.set_k(opdotprec);
	
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(cidot)
#endif
			for(int i=0;i<m1.ysize;i++)
			{
				for(int j=0;j<ms.sxsize;j++)
				{
					cidot=0.0;
					accumulate(cidot,m1[i+Lb(m1,ROW)],ms[Col(j+Lb(ms,COL))]);
					r.dat[i*ms.sxsize+j]=rnd(cidot);
				}
			}
			return r;
		}
	}

	template <class MS,class M,class E>
	TINLINE E _msmmult(const MS &ms, const M &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E r(ms.start1,ms.end1,m2.lb2,m2.ub2);
#if(CXSC_INDEX_CHECK)
		if(ms.sxsize!=m2.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(r)+" operator *(const "+nameof(ms)+" &, const "+nameof(m2)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) {
			blasmatmul(ms,m2,r);
			return r;
		}
		else
#endif
		{
			dotprecision dot(0.0);
			dot.set_k(opdotprec);
	
			int lbr = Lb(ms,ROW);
			int lbc = Lb(m2,COL);
	
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(dot)
#endif
			for(int i=0;i<ms.sysize;i++)
			{
				for(int j=0;j<m2.xsize;j++)
				{
					dot=0.0;
					accumulate_approx(dot,ms[Row(lbr+i)],m2[Col(lbc+j)]);			
					r.dat[i*m2.xsize+j]=rnd(dot);
				}
			}
			return r;
		}
	}

	template <class MS,class M,class E>
	TINLINE E _msmimult(const MS &ms, const M &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E r(ms.start1,ms.end1,m2.lb2,m2.ub2);
#if(CXSC_INDEX_CHECK)
		if(ms.sxsize!=m2.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(r)+" operator *(const "+nameof(ms)+" &, const "+nameof(m2)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) {
			blasmatmul(ms,m2,r);
			return r;
		}
		else
#endif
		{
			idotprecision idot(0.0);
			idot.set_k(opdotprec);
	
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(idot)
#endif
			for(int i=0;i<ms.sysize;i++)
			{
				for(int j=0;j<m2.xsize;j++)
				{
					idot=0.0;
					accumulate(idot,ms[i+Lb(ms,ROW)],m2[Col(j+Lb(m2,COL))]);
					r.dat[i*m2.xsize+j]=rnd(idot);
				}
			}
			return r;
		}
	}

	template <class MS,class M,class E>
	TINLINE E _msmlmult(const MS &ms, const M &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E r(ms.start1,ms.end1,m2.lb2,m2.ub2);
#if(CXSC_INDEX_CHECK)
		if(ms.sxsize!=m2.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(r)+" operator *(const "+nameof(ms)+" &, const "+nameof(m2)+" &)"));
#endif

		dotprecision dot(0.0);

                int lbr = Lb(ms,ROW);
                int lbc = Lb(m2,COL);

#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(dot)
#endif
		for(int i=0;i<ms.sysize;i++)
		{
			for(int j=0;j<m2.xsize;j++)
			{
				dot=0.0;
                                accumulate(dot,ms[Row(lbr+i)],m2[Col(lbc+j)]);			
				r.dat[i*m2.xsize+j]=l_real(dot);
			}
		}
		return r;
	}

	template <class MS,class M,class E>
	TINLINE E _msmlimult(const MS &ms, const M &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E r(ms.start1,ms.end1,m2.lb2,m2.ub2);
#if(CXSC_INDEX_CHECK)
		if(ms.sxsize!=m2.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(r)+" operator *(const "+nameof(ms)+" &, const "+nameof(m2)+" &)"));
#endif
		idotprecision idot(0.0);

#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(idot)
#endif
		for(int i=0;i<ms.sysize;i++)
		{
			for(int j=0;j<m2.xsize;j++)
			{
				idot=0.0;
				for(int k=0;k<ms.sxsize;k++)
					accumulate(idot,ms.dat[(i+ms.offset1)*ms.mxsize+k+ms.offset2],m2.dat[k*m2.xsize+j]);
				r.dat[i*m2.xsize+j]=l_interval(idot);
			}
		}
		return r;
	}

	template <class MS,class M,class E>
	TINLINE E _msmcmult(const MS &ms, const M &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E r(ms.start1,ms.end1,m2.lb2,m2.ub2);
#if(CXSC_INDEX_CHECK)
		if(ms.sxsize!=m2.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(r)+" operator *(const "+nameof(ms)+" &, const "+nameof(m2)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) {
			blasmatmul(ms,m2,r);
			return r;
		}
		else
#endif
		{
			cdotprecision cdot(0.0);
			cdot.set_k(opdotprec);
	
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(cdot)
#endif
			for(int i=0;i<ms.sysize;i++)
			{
				for(int j=0;j<m2.xsize;j++)
				{
					cdot=0.0;
					accumulate_approx(cdot,ms[i+Lb(ms,ROW)],m2[Col(j+Lb(m2,COL))]);
					r.dat[i*m2.xsize+j]=rnd(cdot);
				}
			}
			return r;
		}
	}

	template <class MS,class M,class E>
	TINLINE E _msmcimult(const MS &ms, const M &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E r(ms.start1,ms.end1,m2.lb2,m2.ub2);
#if(CXSC_INDEX_CHECK)
		if(ms.sxsize!=m2.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(r)+" operator *(const "+nameof(ms)+" &, const "+nameof(m2)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) {
			blasmatmul(ms,m2,r);
			return r;
		}
		else
#endif
		{
			cidotprecision cidot(0.0);
			cidot.set_k(opdotprec);
	
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(cidot)
#endif
			for(int i=0;i<ms.sysize;i++)
			{
				for(int j=0;j<m2.xsize;j++)
				{
					cidot=0.0;
					accumulate(cidot,ms[i+Lb(ms,ROW)],m2[Col(j+Lb(m2,COL))]);
					r.dat[i*m2.xsize+j]=rnd(cidot);
				}
			}
			return r;
		}
	}

	template <class M,class MS,class S>
	TINLINE M &_mmsmultassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(m1.xsize!=ms.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(m1)+" &operator *=("+nameof(m1)+" &,const "+nameof(ms)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) 
		{
			M tmp(1,ColLen(m1),1,RowLen(m1));
			blasmatmul(m1,ms,tmp);
			m1 = tmp;
			return m1;
		}
		else
#endif
		{
			S *ndat=new S[m1.ysize*ms.sxsize];
			dotprecision dot(0.0);
			dot.set_k(opdotprec);
	
			int lbr = Lb(ms,ROW);
			int lbc = Lb(m1,COL);
	
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(dot)
#endif
			for(int i=0;i<m1.ysize;i++)
			{
				for(int j=0;j<ms.sxsize;j++)
				{
					dot=0.0;
					accumulate_approx(dot,m1[Row(lbr+i)],ms[Col(lbc+j)]);			
					ndat[i*ms.sxsize+j]=rnd(dot);
				}
			}
			delete [] m1.dat;
			m1.dat=ndat;
			m1.xsize=ms.sxsize;
			m1.lb2=ms.start2;
			m1.ub2=ms.end2;
			return m1;
		}
	}

	template <class M,class MS,class S>
	TINLINE M &_mmsimultassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(m1.xsize!=ms.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(m1)+" &operator *=("+nameof(m1)+" &,const "+nameof(ms)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) 
		{
			M tmp(1,ColLen(m1),1,RowLen(m1));
			blasmatmul(m1,ms,tmp);
			m1 = tmp;
			return m1;
		}
		else
#endif
		{
			S *ndat=new S[m1.ysize*ms.sxsize];
			idotprecision idot(0.0);
			idot.set_k(opdotprec);
	
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(idot)
#endif
			for(int i=0;i<m1.ysize;i++)
			{
				for(int j=0;j<ms.sxsize;j++)
				{
					idot=0.0;
					accumulate(idot,m1[i+Lb(m1,ROW)],ms[Col(j+Lb(ms,COL))]);
					ndat[i*ms.sxsize+j]=rnd(idot);
				}
			}
			delete [] m1.dat;
			m1.dat=ndat;
			m1.xsize=ms.sxsize;
			m1.lb2=ms.start2;
			m1.ub2=ms.end2;
			return m1;
		}
	}

	template <class M,class MS,class S>
	TINLINE M &_mmslmultassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(m1.xsize!=ms.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(m1)+" &operator *=("+nameof(m1)+" &,const "+nameof(ms)+" &)"));
#endif
		S *ndat=new S[m1.ysize*ms.sxsize];

		dotprecision dot(0.0);

                int lbr = Lb(m1,ROW);
                int lbc = Lb(ms,COL);

#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(dot)
#endif
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<ms.sxsize;j++)
			{
				dot=0.0;
                                accumulate(dot,m1[Row(lbr+i)],ms[Col(lbc+j)]);			
				ndat[i*ms.sxsize+j]=l_real(dot);
			}
		}
		delete [] m1.dat;
		m1.dat=ndat;
		m1.xsize=ms.sxsize;
		m1.lb2=ms.start2;
		m1.ub2=ms.end2;
		return m1;
	}

	template <class M,class MS,class S>
	TINLINE M &_mmslimultassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(m1.xsize!=ms.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(m1)+" &operator *=("+nameof(m1)+" &,const "+nameof(ms)+" &)"));
#endif
		S *ndat=new S[m1.ysize*ms.sxsize];
		idotprecision idot(0.0);

#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(idot)
#endif
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<ms.sxsize;j++)
			{
				idot=0.0;
				for(int k=0;k<m1.xsize;k++)
					accumulate(idot,m1.dat[i*m1.xsize+k],ms.dat[(k+ms.offset1)*ms.mxsize+j+ms.offset2]);
				ndat[i*ms.sxsize+j]=l_interval(idot);
			}
		}
		delete [] m1.dat;
		m1.dat=ndat;
		m1.xsize=ms.sxsize;
		m1.lb2=ms.start2;
		m1.ub2=ms.end2;
		return m1;
	}

	template <class M,class MS,class S>
	TINLINE M &_mmscmultassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(m1.xsize!=ms.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(m1)+" &operator *=("+nameof(m1)+" &,const "+nameof(ms)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) 
		{
			M tmp(1,ColLen(m1),1,RowLen(m1));
			blasmatmul(m1,ms,tmp);
			m1 = tmp;
			return m1;
		}
		else
#endif
		{
			S *ndat=new S[m1.ysize*ms.sxsize];
			cdotprecision cdot(0.0);
			cdot.set_k(opdotprec);
	
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(cdot)
#endif
			for(int i=0;i<m1.ysize;i++)
			{
				for(int j=0;j<ms.sxsize;j++)
				{
					cdot=0.0;
					accumulate_approx(cdot,m1[i+Lb(m1,ROW)],ms[Col(j+Lb(ms,COL))]);
					ndat[i*ms.sxsize+j]=rnd(cdot);
				}
			}
			delete [] m1.dat;
			m1.dat=ndat;
			m1.xsize=ms.sxsize;
			m1.lb2=ms.start2;
			m1.ub2=ms.end2;
			return m1;
		}
	}

	template <class M,class MS,class S>
	TINLINE M &_mmscimultassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(m1.xsize!=ms.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(m1)+" &operator *=("+nameof(m1)+" &,const "+nameof(ms)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) 
		{
			M tmp(1,ColLen(m1),1,RowLen(m1));
			blasmatmul(m1,ms,tmp);
			m1 = tmp;
			return m1;
		}
		else
#endif
		{
			S *ndat=new S[m1.ysize*ms.sxsize];
			cidotprecision cidot(0.0);
			cidot.set_k(opdotprec);
	
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(cidot)
#endif
			for(int i=0;i<m1.ysize;i++)
			{
				for(int j=0;j<ms.sxsize;j++)
				{
					cidot=0.0;
					accumulate(cidot,m1[i+Lb(m1,ROW)],ms[Col(j+Lb(ms,COL))]);
					ndat[i*ms.sxsize+j]=rnd(cidot);
				}
			}
			delete [] m1.dat;
			m1.dat=ndat;
			m1.xsize=ms.sxsize;
			m1.lb2=ms.start2;
			m1.ub2=ms.end2;
			return m1;
		}
	}

	template <class MS1,class MS2,class E>
	TINLINE E _msmsmult(const MS1 &ms1, const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E r(ms1.start1,ms1.end1,ms2.start2,ms2.end2);
#if(CXSC_INDEX_CHECK)
		if(ms1.sxsize!=ms2.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(r)+" operator *(const "+nameof(ms1)+" &, const "+nameof(ms2)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) {
			blasmatmul(ms1,ms2,r);
			return r;
		}
		else
#endif
		{
			dotprecision dot(0.0);
			dot.set_k(opdotprec);
	
			int lbr = Lb(ms1,ROW);
			int lbc = Lb(ms2,COL);
	
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(dot)
#endif
			for(int i=0;i<ms1.sysize;i++)
			{
				for(int j=0;j<ms2.sxsize;j++)
				{
					dot=0.0;
					accumulate_approx(dot,ms1[Row(lbr+i)],ms2[Col(lbc+j)]);			
					r.dat[i*ms2.sxsize+j]=rnd(dot);
				}
			}
			return r;
		}
	}
	
	template <class MS1,class MS2,class E>
	TINLINE E _msmsimult(const MS1 &ms1, const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E r(ms1.start1,ms1.end1,ms2.start2,ms2.end2);
#if(CXSC_INDEX_CHECK)
		if(ms1.sxsize!=ms2.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(r)+" operator *(const "+nameof(ms1)+" &, const "+nameof(ms2)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) {
			blasmatmul(ms1,ms2,r);
			return r;
		}
		else
#endif
		{
			idotprecision idot(0.0);
			idot.set_k(opdotprec);
	
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(idot)
#endif
			for(int i=0;i<ms1.sysize;i++)
			{
				for(int j=0;j<ms2.sxsize;j++)
				{
					idot=0.0;
					accumulate(idot,ms1[i+Lb(ms1,ROW)],ms2[Col(j+Lb(ms2,COL))]);
					r.dat[i*ms2.sxsize+j]=rnd(idot);
				}
			}
			return r;
		}
	}
	
	template <class MS1,class MS2,class E>
	TINLINE E _msmslmult(const MS1 &ms1, const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E r(ms1.start1,ms1.end1,ms2.start2,ms2.end2);
#if(CXSC_INDEX_CHECK)
		if(ms1.sxsize!=ms2.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(r)+" operator *(const "+nameof(ms1)+" &, const "+nameof(ms2)+" &)"));
#endif
		dotprecision dot(0.0);

                int lbr = Lb(ms1,ROW);
                int lbc = Lb(ms2,COL);

#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(dot)
#endif
		for(int i=0;i<ms1.sysize;i++)
		{
			for(int j=0;j<ms2.sxsize;j++)
			{
				dot=0.0;
                                accumulate(dot,ms1[Row(lbr+i)],ms2[Col(lbc+j)]);			
				r.dat[i*ms2.sxsize+j]=l_real(dot);
			}
		}
		return r;
	}
	
	template <class MS1,class MS2,class E>
	TINLINE E _msmslimult(const MS1 &ms1, const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E r(ms1.start1,ms1.end1,ms2.start2,ms2.end2);
#if(CXSC_INDEX_CHECK)
		if(ms1.sxsize!=ms2.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(r)+" operator *(const "+nameof(ms1)+" &, const "+nameof(ms2)+" &)"));
#endif
		idotprecision idot(0.0);

#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(idot)
#endif
		for(int i=0;i<ms1.sysize;i++)
		{
			for(int j=0;j<ms2.sxsize;j++)
			{
				idot=0.0;
				for(int k=0;k<ms1.sxsize;k++)
					accumulate(idot,ms1.dat[(i+ms1.offset1)*ms1.mxsize+k+ms1.offset2],ms2.dat[(k+ms2.offset1)*ms2.mxsize+j+ms2.offset2]);
				r.dat[i*ms2.sxsize+j]=l_interval(idot);
			}
		}
		return r;
	}
	
	template <class MS1,class MS2,class E>
	TINLINE E _msmscmult(const MS1 &ms1, const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E r(ms1.start1,ms1.end1,ms2.start2,ms2.end2);
#if(CXSC_INDEX_CHECK)
		if(ms1.sxsize!=ms2.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(r)+" operator *(const "+nameof(ms1)+" &, const "+nameof(ms2)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) {
			blasmatmul(ms1,ms2,r);
			return r;
		}
		else
#endif
		{
			cdotprecision cdot(0.0);
			cdot.set_k(opdotprec);
	
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(cdot)
#endif
			for(int i=0;i<ms1.sysize;i++)
			{
				for(int j=0;j<ms2.sxsize;j++)
				{
					cdot=0.0;
					accumulate_approx(cdot,ms1[i+Lb(ms1,ROW)],ms2[Col(j+Lb(ms2,COL))]);
					r.dat[i*ms2.sxsize+j]=rnd(cdot);
				}
			}
			return r;
		}
	}
	
	template <class MS1,class MS2,class E>
	TINLINE E _msmscimult(const MS1 &ms1, const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E r(ms1.start1,ms1.end1,ms2.start2,ms2.end2);
#if(CXSC_INDEX_CHECK)
		if(ms1.sxsize!=ms2.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(r)+" operator *(const "+nameof(ms1)+" &, const "+nameof(ms2)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) {
			blasmatmul(ms1,ms2,r);
			return r;
		}
		else
#endif
		{
			cidotprecision cidot(0.0);
			cidot.set_k(opdotprec);
	
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(cidot)
#endif
			for(int i=0;i<ms1.sysize;i++)
			{
				for(int j=0;j<ms2.sxsize;j++)
				{
					cidot=0.0;
					accumulate(cidot,ms1[Lb(ms1,ROW)+i],ms2[Col(j+Lb(ms2,COL))]);
					r.dat[i*ms2.sxsize+j]=rnd(cidot);
				}
			}
			return r;
		}
	}
	
	template <class S,class M,class E>
	TINLINE E _smmult(const S &c, const M &m) throw()
	{
		E r(m.lb1,m.ub1,m.lb2,m.ub2);

		for(int i=0;i<m.xsize*m.ysize;i++)
			r.dat[i]=c*m.dat[i];
		return r;
	}

	template <class M,class S>
	TINLINE M &_msmultassign(M &m,const S &c) throw()
	{
		for(int i=0;i<m.xsize*m.ysize;i++)
			m.dat[i]*=c;
		return m;
	}
	
	template <class S,class MS,class E>
	TINLINE E _smsmult(const S &c, const MS &ms) throw()
	{
		E r(ms.start1,ms.end1,ms.start2,ms.end2);

		for(int i=0;i<ms.sxsize;i++)
		{
			for(int j=0;j<ms.sysize;j++)
				r.dat[i*ms.sxsize+j]=c*ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2];
		}
		return r;
	}

	template <class MS,class S>
	TINLINE MS &_mssmultassign(MS &ms,const S &c) throw()
	{
		for(int i=0;i<ms.sxsize;i++)
		{
			for(int j=0;j<ms.sysize;j++)
				ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2]*=c;
		}
		return ms;
	}

	template <class M,class V,class E>
	TINLINE E _mvmult(const M &m,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
		E r(m.lb1,m.ub1);
#if(CXSC_INDEX_CHECK)
		if(v.size!=m.xsize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(r)+" operator *(const "+nameof(m)+" &,const "+nameof(v)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) {
			blasmvmul(m,v,r);
			return r;
		}
		else
#endif
		{
			dotprecision dot(0.0);
			dot.set_k(opdotprec);
	
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(dot)
#endif
			for(int i=0;i<m.ysize;i++)
			{
				dot=0.0;
				accumulate_approx(dot,m[Row(Lb(m,ROW)+i)],v);
				r.dat[i]=rnd(dot);
			}
			return r;
		}
	}

	template <class M,class V,class E>
	TINLINE E _mvimult(const M &m,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
		E r(m.lb1,m.ub1);
#if(CXSC_INDEX_CHECK)
		if(v.size!=m.xsize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(r)+" operator *(const "+nameof(m)+" &,const "+nameof(v)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) {
			blasmvmul(m,v,r);
			return r;
		}
		else
#endif
		{
			idotprecision idot(0.0);
			idot.set_k(opdotprec);
	
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(idot)
#endif
			for(int i=0;i<m.ysize;i++)
			{
				idot=0.0;
				accumulate(idot,m[i+Lb(m,ROW)],v);
				r.dat[i]=rnd(idot);
			}
			return r;
		}
	}

	template <class M,class V,class E>
	TINLINE E _mvlmult(const M &m,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
		E r(m.lb1,m.ub1);
#if(CXSC_INDEX_CHECK)
		if(v.size!=m.xsize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(r)+" operator *(const "+nameof(m)+" &,const "+nameof(v)+" &)"));
#endif

		dotprecision dot(0.0);

#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(dot)
#endif
		for(int i=0;i<m.ysize;i++)
		{
			dot=0.0;
                        accumulate(dot,m[Row(Lb(m,ROW)+i)],v);
			r.dat[i]=l_real(dot);
		}
		return r;
	}

	template <class M,class V,class E>
	TINLINE E _mvlimult(const M &m,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
		E r(m.lb1,m.ub1);
#if(CXSC_INDEX_CHECK)
		if(v.size!=m.xsize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(r)+" operator *(const "+nameof(m)+" &,const "+nameof(v)+" &)"));
#endif
		idotprecision idot(0.0);

#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(idot)
#endif
		for(int i=0;i<m.ysize;i++)
		{
			idot=0.0;
			for(int j=0;j<m.xsize;j++)
				accumulate(idot,m.dat[i*m.xsize+j],v.dat[j]);
			r.dat[i]=l_interval(idot);
		}
		return r;
	}

	template <class M,class V,class E>
	TINLINE E _mvcmult(const M &m,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
		E r(m.lb1,m.ub1);
#if(CXSC_INDEX_CHECK)
		if(v.size!=m.xsize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(r)+" operator *(const "+nameof(m)+" &,const "+nameof(v)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) {
			blasmvmul(m,v,r);
			return r;
		}
		else
#endif
		{
			cdotprecision cdot(0.0);
			cdot.set_k(opdotprec);
	
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(cdot)
#endif
			for(int i=0;i<m.ysize;i++)
			{
				cdot=0.0;
				accumulate_approx(cdot,m[i+Lb(m,ROW)],v);
				r.dat[i]=rnd(cdot);
			}
			return r;
		}
	}

	template <class M,class V,class E>
	TINLINE E _mvcimult(const M &m,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
		E r(m.lb1,m.ub1);
#if(CXSC_INDEX_CHECK)
		if(v.size!=m.xsize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(r)+" operator *(const "+nameof(m)+" &,const "+nameof(v)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) {
			blasmvmul(m,v,r);
			return r;
		}
		else
#endif
		{
			cidotprecision cidot(0.0);
			cidot.set_k(opdotprec);
	
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(cidot)
#endif
			for(int i=0;i<m.ysize;i++)
			{
				cidot=0.0;
				accumulate(cidot,m[i+Lb(m,ROW)],v);
				r.dat[i]=rnd(cidot);
			}
			return r;
		}
	}

	template <class V,class M,class E>
	TINLINE E _vmmult(const V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
		E r(m.lb2,m.ub2);
#if(CXSC_INDEX_CHECK)
		if(v.size!=m.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(r)+" operator *(const "+nameof(v)+" &,const "+nameof(m)+" &)"));
#endif
		dotprecision dot(0.0);
                dot.set_k(opdotprec);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(dot)
#endif
		for(int i=0;i<m.xsize;i++)
		{
			dot=0.0;
                        accumulate_approx(dot,v,m[Col(Lb(m,COL)+i)]);
			r.dat[i]=rnd(dot);
		}
		return r;
	}

	template <class V,class M,class E>
	TINLINE E _vmimult(const V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
		E r(m.lb2,m.ub2);
#if(CXSC_INDEX_CHECK)
		if(v.size!=m.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(r)+" operator *(const "+nameof(v)+" &,const "+nameof(m)+" &)"));
#endif
		idotprecision idot(0.0);
		idot.set_k(opdotprec);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(idot)
#endif
		for(int i=0;i<m.xsize;i++)
		{
			idot=0.0;
			accumulate(idot,v,m[Col(Lb(m,COL)+i)]);
			r.dat[i]=rnd(idot);
		}
		return r;
	}

	template <class V,class M,class E>
	TINLINE E _vmcmult(const V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
		E r(m.lb2,m.ub2);
#if(CXSC_INDEX_CHECK)
		if(v.size!=m.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(r)+" operator *(const "+nameof(v)+" &,const "+nameof(m)+" &)"));
#endif
		cdotprecision cdot(0.0);
		cdot.set_k(opdotprec);

#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(cdot)
#endif
		for(int i=0;i<m.xsize;i++)
		{
			cdot=0.0;
			accumulate_approx(cdot,v,m[Col(i+Lb(m,COL))]);
			r.dat[i]=rnd(cdot);
		}
		return r;
	}

	template <class V,class M,class E>
	TINLINE E _vmcimult(const V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
		E r(m.lb2,m.ub2);
#if(CXSC_INDEX_CHECK)
		if(v.size!=m.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(r)+" operator *(const "+nameof(v)+" &,const "+nameof(m)+" &)"));
#endif
		cidotprecision cidot(0.0);
		cidot.set_k(opdotprec);

#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(cidot)
#endif
		for(int i=0;i<m.xsize;i++)
		{
			cidot=0.0;
			accumulate(cidot,v,m[Col(i+Lb(m,COL))]);
			r.dat[i]=rnd(cidot);
		}
		return r;
	}

	template <class V,class M,class E>
	TINLINE E _vmlmult(const V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
		E r(m.lb2,m.ub2);
#if(CXSC_INDEX_CHECK)
		if(v.size!=m.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(r)+" operator *(const "+nameof(v)+" &,const "+nameof(m)+" &)"));
#endif
		dotprecision dot(0.0);

#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(dot)
#endif
		for(int i=0;i<m.xsize;i++)
		{
			dot=0.0;
                        accumulate(dot,v,m[Col(Lb(m,COL)+i)]);
			r.dat[i]=l_real(dot);
		}
		return r;
	}

	template <class V,class M,class E>
	TINLINE E _vmlimult(const V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
		E r(m.lb2,m.ub2);
#if(CXSC_INDEX_CHECK)
		if(v.size!=m.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(r)+" operator *(const "+nameof(v)+" &,const "+nameof(m)+" &)"));
#endif
		idotprecision idot(0.0);

#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(idot)
#endif
		for(int i=0;i<m.xsize;i++)
		{
			idot=0.0;
			for(int j=0;j<m.ysize;j++)
				accumulate(idot,m.dat[i+m.xsize*j],v.dat[j]);
			r.dat[i]=l_interval(idot);
		}
		return r;
	}

	template <class V,class M,class S>
	TINLINE V &_vmmultassign(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(v.size!=m.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(v)+" &operator *=("+nameof(v)+" &,const "+nameof(m)+" &)"));
#endif
		v.size=m.xsize;
		S *ndat=new S[v.size];
		v.l=m.lb2;
		v.u=m.ub2;

		dotprecision dot(0.0);
                dot.set_k(opdotprec);

#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(dot)
#endif
		for(int i=0;i<m.xsize;i++)
		{
			dot=0.0;
			accumulate_approx(dot,v,m[Col(i+Lb(m,COL))]);
			ndat[i]=rnd(dot);
		}
		delete [] v.dat;
		v.dat=ndat;
		return v;
	}

	template <class V,class M,class S>
	TINLINE V &_vmimultassign(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(v.size!=m.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(v)+" &operator *=("+nameof(v)+" &,const "+nameof(m)+" &)"));
#endif
		v.size=m.xsize;
		S *ndat=new S[v.size];
		v.l=m.lb2;
		v.u=m.ub2;
		idotprecision idot(0.0);
		idot.set_k(opdotprec);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(idot)
#endif
		for(int i=0;i<m.xsize;i++)
		{
			idot=0.0;
			accumulate(idot,v,m[Col(i+Lb(m,COL))]);
			ndat[i]=rnd(idot);
		}
		delete [] v.dat;
		v.dat=ndat;
		return v;
	}

	template <class V,class M,class S>
	TINLINE V &_vmlmultassign(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(v.size!=m.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(v)+" &operator *=("+nameof(v)+" &,const "+nameof(m)+" &)"));
#endif
		v.size=m.xsize;
		S *ndat=new S[v.size];
		v.l=m.lb2;
		v.u=m.ub2;
		dotprecision dot(0.0);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(dot)
#endif		
		for(int i=0;i<m.xsize;i++)
		{
			dot=0.0;
			for(int j=0;j<m.ysize;j++)
				accumulate(dot,m.dat[i+m.xsize*j],v.dat[j]);
			ndat[i]=S(dot);
		}
		delete [] v.dat;
		v.dat=ndat;
		return v;
	}

	template <class V,class M,class S>
	TINLINE V &_vmlimultassign(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(v.size!=m.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(v)+" &operator *=("+nameof(v)+" &,const "+nameof(m)+" &)"));
#endif
		v.size=m.xsize;
		S *ndat=new S[v.size];
		v.l=m.lb2;
		v.u=m.ub2;
		idotprecision idot(0.0);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(idot)
#endif		
		for(int i=0;i<m.xsize;i++)
		{
			idot=0.0;
			for(int j=0;j<m.ysize;j++)
				accumulate(idot,m.dat[i+m.xsize*j],v.dat[j]);
			ndat[i]=S(idot);
		}
		delete [] v.dat;
		v.dat=ndat;
		return v;
	}

	template <class V,class M,class S>
	TINLINE V &_vmcmultassign(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(v.size!=m.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(v)+" &operator *=("+nameof(v)+" &,const "+nameof(m)+" &)"));
#endif
		v.size=m.xsize;
		S *ndat=new S[v.size];
		v.l=m.lb2;
		v.u=m.ub2;
		cdotprecision cdot(0.0);
		cdot.set_k(opdotprec);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(cdot)
#endif
		for(int i=0;i<m.xsize;i++)
		{
			cdot=0.0;
			accumulate_approx(cdot,v,m[Col(i+Lb(m,COL))]);
			ndat[i]=rnd(cdot);
		}
		delete [] v.dat;
		v.dat=ndat;
		return v;
	}

	template <class V,class M,class S>
	TINLINE V &_vmcimultassign(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(v.size!=m.ysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(v)+" &operator *=("+nameof(v)+" &,const "+nameof(m)+" &)"));
#endif
		v.size=m.xsize;
		S *ndat=new S[v.size];
		v.l=m.lb2;
		v.u=m.ub2;
		cidotprecision cidot(0.0);
		cidot.set_k(opdotprec);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(cidot)
#endif		
		for(int i=0;i<m.xsize;i++)
		{
			cidot=0.0;
			accumulate(cidot,v,m[Col(i+Lb(m,COL))]);
			ndat[i]=rnd(cidot);
		}
		delete [] v.dat;
		v.dat=ndat;
		return v;
	}

	template <class VS,class M,class S>
	TINLINE VS &_vsmmultassign(VS &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((v.size!=m.ysize)||(v.size!=m.xsize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(v)+" &operator *=("+nameof(v)+" &,const "+nameof(m)+" &)"));
#endif
		S *ndat=new S[v.size];
		dotprecision dot(0.0);
                dot.set_k(opdotprec);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(dot)
#endif		
		for(int i=0;i<m.xsize;i++)
		{
			dot=0.0;
			accumulate_approx(dot,v,m[Col(i+Lb(m,COL))]);
			ndat[i]=rnd(dot);
		}
		for(int i=v.start-v.l,j=0;j<m.xsize;j++,i++)
			v.dat[i]=ndat[j];
		return v;
	}

	template <class VS,class M,class S>
	TINLINE VS &_vsmimultassign(VS &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((v.size!=m.ysize)||(v.size!=m.xsize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(v)+" &operator *=("+nameof(v)+" &,const "+nameof(m)+" &)"));
#endif
		S *ndat=new S[v.size];
		idotprecision idot(0.0);
		idot.set_k(opdotprec);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(idot)
#endif		
		for(int i=0;i<m.xsize;i++)
		{
			idot=0.0;
			accumulate(idot,v,m[Col(i+Lb(m,COL))]);
			ndat[i]=rnd(idot);
		}
		for(int i=v.start-v.l,j=0;j<m.xsize;j++,i++)
			v.dat[i]=ndat[j];
		return v;
	}

	template <class VS,class M,class S>
	TINLINE VS &_vsmlmultassign(VS &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((v.size!=m.ysize)||(v.size!=m.xsize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(v)+" &operator *=("+nameof(v)+" &,const "+nameof(m)+" &)"));
#endif
		S *ndat=new S[v.size];
		dotprecision dot(0.0);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(dot)
#endif
		for(int i=0;i<m.xsize;i++)
		{
			dot=0.0;
			for(int j=0,k=v.start-v.l;j<m.ysize;j++,k++)
				accumulate(dot,m.dat[i+m.xsize*j],v.dat[k]);
			ndat[i]=S(dot);
		}
		for(int i=v.start-v.l,j=0;j<m.xsize;j++,i++)
			v.dat[i]=ndat[j];
		return v;
	}

	template <class VS,class M,class S>
	TINLINE VS &_vsmlimultassign(VS &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((v.size!=m.ysize)||(v.size!=m.xsize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(v)+" &operator *=("+nameof(v)+" &,const "+nameof(m)+" &)"));
#endif
		S *ndat=new S[v.size];
		idotprecision idot(0.0);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(idot)
#endif		
		for(int i=0;i<m.xsize;i++)
		{
			idot=0.0;
			for(int j=0,k=v.start-v.l;j<m.ysize;j++,k++)
				accumulate(idot,m.dat[i+m.xsize*j],v.dat[k]);
			ndat[i]=S(idot);
		}
		for(int i=v.start-v.l,j=0;j<m.xsize;j++,i++)
			v.dat[i]=ndat[j];
		return v;
	}

	template <class VS,class M,class S>
	TINLINE VS &_vsmcmultassign(VS &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((v.size!=m.ysize)||(v.size!=m.xsize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(v)+" &operator *=("+nameof(v)+" &,const "+nameof(m)+" &)"));
#endif
		S *ndat=new S[v.size];
		cdotprecision cdot(0.0);
		cdot.set_k(opdotprec);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(cdot)
#endif
		for(int i=0;i<m.xsize;i++)
		{
			cdot=0.0;
			accumulate_approx(cdot,v,m[Col(i+Lb(m,COL))]);
			ndat[i]=rnd(cdot);
		}
		for(int i=v.start-v.l,j=0;j<m.xsize;j++,i++)
			v.dat[i]=ndat[j];
		return v;
	}

	template <class VS,class M,class S>
	TINLINE VS &_vsmcimultassign(VS &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((v.size!=m.ysize)||(v.size!=m.xsize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<M>(nameof(v)+" &operator *=("+nameof(v)+" &,const "+nameof(m)+" &)"));
#endif
		S *ndat=new S[v.size];
		cidotprecision cidot(0.0);
		cidot.set_k(opdotprec);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(cidot)
#endif		
		for(int i=0;i<m.xsize;i++)
		{
			cidot=0.0;
			accumulate(cidot,v,m[Col(i+Lb(m,COL))]);
			ndat[i]=rnd(cidot);
		}
		for(int i=v.start-v.l,j=0;j<m.xsize;j++,i++)
			v.dat[i]=ndat[j];
		return v;
	}

	template <class MS,class V,class E>
	TINLINE E _msvmult(const MS &ms,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>)
#else
	throw()
#endif
	{
		E r(ms.start1,ms.end1);
#if(CXSC_INDEX_CHECK)
		if(v.size!=ms.sxsize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS>(nameof(r)+" operator *(const "+nameof(ms)+" &,const "+nameof(v)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) {
			blasmvmul(ms,v,r);
			return r;
		}
		else
#endif
		{
			dotprecision dot(0.0);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(dot)
#endif			
			for(int i=0;i<ms.sysize;i++)
			{
				dot=0.0;
				accumulate_approx(dot,ms[i+Lb(ms,ROW)],v);
				r.dat[i]=rnd(dot);
			}
			return r;
		}
	}

	template <class MS,class V,class E>
	TINLINE E _msvimult(const MS &ms,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>)
#else
	throw()
#endif
	{
		E r(ms.start1,ms.end1);
#if(CXSC_INDEX_CHECK)
		if(v.size!=ms.sxsize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS>(nameof(r)+" operator *(const "+nameof(ms)+" &,const "+nameof(v)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) {
			blasmvmul(ms,v,r);
			return r;
		}
		else
#endif
		{
			idotprecision idot(0.0);
			idot.set_k(opdotprec);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(idot)
#endif			
			for(int i=0;i<ms.sysize;i++)
			{
				idot=0.0;
				accumulate(idot,ms[i+Lb(ms,ROW)],v);
				r.dat[i]=rnd(idot);
			}
			return r;
		}
	}

	template <class MS,class V,class E>
	TINLINE E _msvlmult(const MS &ms,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>)
#else
	throw()
#endif
	{
		E r(ms.start1,ms.end1);
#if(CXSC_INDEX_CHECK)
		if(v.size!=ms.sxsize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS>(nameof(r)+" operator *(const "+nameof(ms)+" &,const "+nameof(v)+" &)"));
#endif
		dotprecision dot(0.0);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(dot)
#endif		
		for(int i=0;i<ms.sysize;i++)
		{
			dot=0.0;
			for(int j=0;j<ms.sxsize;j++)
				accumulate(dot,ms.dat[i*ms.mxsize+j],v.dat[j]);
			r.dat[i]=l_real(dot);
		}
		return r;
	}

	template <class MS,class V,class E>
	TINLINE E _msvlimult(const MS &ms,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>)
#else
	throw()
#endif
	{
		E r(ms.start1,ms.end1);
#if(CXSC_INDEX_CHECK)
		if(v.size!=ms.sxsize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS>(nameof(r)+" operator *(const "+nameof(ms)+" &,const "+nameof(v)+" &)"));
#endif
		idotprecision idot(0.0);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(idot)
#endif
		for(int i=0;i<ms.sysize;i++)
		{
			idot=0.0;
			for(int j=0;j<ms.sxsize;j++)
				accumulate(idot,ms.dat[i*ms.mxsize+j],v.dat[j]);
			r.dat[i]=l_interval(idot);
		}
		return r;
	}

	template <class MS,class V,class E>
	TINLINE E _msvcmult(const MS &ms,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>)
#else
	throw()
#endif
	{
		E r(ms.start1,ms.end1);
#if(CXSC_INDEX_CHECK)
		if(v.size!=ms.sxsize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS>(nameof(r)+" operator *(const "+nameof(ms)+" &,const "+nameof(v)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) {
			blasmvmul(ms,v,r);
			return r;
		}
		else
#endif
		{	
			cdotprecision cdot(0.0);
			cdot.set_k(opdotprec);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(cdot)
#endif	
			for(int i=0;i<ms.sysize;i++)
			{
				cdot=0.0;
				accumulate_approx(cdot,ms[i+Lb(ms,ROW)],v);
				r.dat[i]=rnd(cdot);
			}
			return r;
		}
	}

	template <class MS,class V,class E>
	TINLINE E _msvcimult(const MS &ms,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>)
#else
	throw()
#endif
	{
		E r(ms.start1,ms.end1);
#if(CXSC_INDEX_CHECK)
		if(v.size!=ms.sxsize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS>(nameof(r)+" operator *(const "+nameof(ms)+" &,const "+nameof(v)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
		if(opdotprec == 1) {
			blasmvmul(ms,v,r);
			return r;
		}
		else
#endif
		{
			cidotprecision cidot(0.0);
			cidot.set_k(opdotprec);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(cidot)
#endif	
			for(int i=0;i<ms.sysize;i++)
			{
				cidot=0.0;
				accumulate(cidot,ms[i+Lb(ms,ROW)],v);
				r.dat[i]=rnd(cidot);
			}
			return r;
		}
	}

	template <class V,class MS,class E>
	TINLINE E _vmsmult(const V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>)
#else
	throw()
#endif
	{
		E r(ms.start2,ms.end2);
#if(CXSC_INDEX_CHECK)
		if(v.size!=ms.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS>(nameof(r)+" operator *(const "+nameof(v)+" &,const "+nameof(ms)+" &)"));
#endif
		dotprecision dot(0.0);
                dot.set_k(opdotprec);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(dot)
#endif		
		for(int i=0;i<ms.sxsize;i++)
		{
			dot=0.0;
			accumulate_approx(dot,v,ms[Col(i+Lb(ms,COL))]);
			r.dat[i]=rnd(dot);
		}
		return r;
	}

	template <class V,class MS,class E>
	TINLINE E _vmsimult(const V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>)
#else
	throw()
#endif
	{
		E r(ms.start2,ms.end2);
#if(CXSC_INDEX_CHECK)
		if(v.size!=ms.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS>(nameof(r)+" operator *(const "+nameof(v)+" &,const "+nameof(ms)+" &)"));
#endif
		idotprecision idot(0.0);
		idot.set_k(opdotprec);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(idot)
#endif		
		for(int i=0;i<ms.sxsize;i++)
		{
			idot=0.0;
			for(int j=0;j<ms.sysize;j++)
				accumulate(idot,v,ms[Col(i+Lb(ms,COL))]);
			r.dat[i]=rnd(idot);
		}
		return r;
	}

	template <class V,class MS,class E>
	TINLINE E _vmscmult(const V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>)
#else
	throw()
#endif
	{
		E r(ms.start2,ms.end2);
#if(CXSC_INDEX_CHECK)
		if(v.size!=ms.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS>(nameof(r)+" operator *(const "+nameof(v)+" &,const "+nameof(ms)+" &)"));
#endif
		cdotprecision cdot(0.0);
		cdot.set_k(opdotprec);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(cdot)
#endif
		for(int i=0;i<ms.sxsize;i++)
		{
			cdot=0.0;
			accumulate_approx(cdot,v,ms[Col(i+Lb(ms,COL))]);
			r.dat[i]=rnd(cdot);
		}
		return r;
	}

	template <class V,class MS,class E>
	TINLINE E _vmscimult(const V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>)
#else
	throw()
#endif
	{
		E r(ms.start2,ms.end2);
#if(CXSC_INDEX_CHECK)
		if(v.size!=ms.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS>(nameof(r)+" operator *(const "+nameof(v)+" &,const "+nameof(ms)+" &)"));
#endif
		cidotprecision cidot(0.0);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(cidot)
#endif		
		for(int i=0;i<ms.sxsize;i++)
		{
			cidot=0.0;
			accumulate(cidot,v,ms[Col(i+Lb(ms,COL))]);
			r.dat[i]=rnd(cidot);
		}
		return r;
	}

	template <class V,class MS,class E>
	TINLINE E _vmslmult(const V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>)
#else
	throw()
#endif
	{
		E r(ms.start2,ms.end2);
#if(CXSC_INDEX_CHECK)
		if(v.size!=ms.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS>(nameof(r)+" operator *(const "+nameof(v)+" &,const "+nameof(ms)+" &)"));
#endif
		dotprecision dot(0.0);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(dot)
#endif		
		for(int i=0;i<ms.sxsize;i++)
		{
			dot=0.0;
			for(int j=0;j<ms.sysize;j++)
				accumulate(dot,ms.dat[i+ms.mxsize*j],v.dat[j]);
			r.dat[i]=l_real(dot);
		}
		return r;
	}

	template <class V,class MS,class E>
	TINLINE E _vmslimult(const V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>)
#else
	throw()
#endif
	{
		E r(ms.start2,ms.end2);
#if(CXSC_INDEX_CHECK)
		if(v.size!=ms.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS>(nameof(r)+" operator *(const "+nameof(v)+" &,const "+nameof(ms)+" &)"));
#endif
		idotprecision idot(0.0);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(idot)
#endif		
		for(int i=0;i<ms.sxsize;i++)
		{
			idot=0.0;
			for(int j=0;j<ms.sysize;j++)
				accumulate(idot,ms.dat[i+ms.mxsize*j],v.dat[j]);
			r.dat[i]=l_interval(idot);
		}
		return r;
	}

	template <class V,class MS,class S>
	TINLINE V &_vmsmultassign(V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(v.size!=ms.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS>(nameof(v)+" &operator *=("+nameof(v)+" &,const "+nameof(ms)+" &)"));
#endif
		v.size=ms.sxsize;
		S *ndat=new S[v.size];
		v.l=ms.start2;
		v.u=ms.end2;
		dotprecision dot(0.0);
                dot.set_k(opdotprec);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(dot)
#endif		
		for(int i=0;i<ms.sxsize;i++)
		{
			dot=0.0;
			accumulate_approx(dot,v,ms[Col(i+Lb(ms,COL))]);
			ndat[i]=rnd(dot);
		}
		delete [] v.dat;
		v.dat=ndat;
		return v;
	}

	template <class V,class MS,class S>
	TINLINE V &_vmsimultassign(V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(v.size!=ms.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS>(nameof(v)+" &operator *=("+nameof(v)+" &,const "+nameof(ms)+" &)"));
#endif
		v.size=ms.sxsize;
		S *ndat=new S[v.size];
		v.l=ms.start2;
		v.u=ms.end2;
		idotprecision idot(0.0);
		idot.set_k(opdotprec);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(idot)
#endif		
		for(int i=0;i<ms.sxsize;i++)
		{
			idot=0.0;
			accumulate(idot,v,ms[Col(i+Lb(ms,COL))]);
			ndat[i]=rnd(idot);
		}
		delete [] v.dat;
		v.dat=ndat;
		return v;
	}

	template <class V,class MS,class S>
	TINLINE V &_vmslmultassign(V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(v.size!=ms.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS>(nameof(v)+" &operator *=("+nameof(v)+" &,const "+nameof(ms)+" &)"));
#endif
		v.size=ms.sxsize;
		S *ndat=new S[v.size];
		v.l=ms.start2;
		v.u=ms.end2;
		dotprecision dot(0.0);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(dot)
#endif		
		for(int i=0;i<ms.sxsize;i++)
		{
			dot=0.0;
			for(int j=0;j<ms.sysize;j++)
				accumulate(dot,ms.dat[i+ms.mxsize*j],v.dat[j]);
			ndat[i]=S(dot);
		}
		delete [] v.dat;
		v.dat=ndat;
		return v;
	}

	template <class V,class MS,class S>
	TINLINE V &_vmslimultassign(V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(v.size!=ms.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS>(nameof(v)+" &operator *=("+nameof(v)+" &,const "+nameof(ms)+" &)"));
#endif
		v.size=ms.sxsize;
		S *ndat=new S[v.size];
		v.l=ms.start2;
		v.u=ms.end2;
		idotprecision idot(0.0);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(idot)
#endif		
		for(int i=0;i<ms.sxsize;i++)
		{
			idot=0.0;
			for(int j=0;j<ms.sysize;j++)
				accumulate(idot,ms.dat[i+ms.mxsize*j],v.dat[j]);
			ndat[i]=S(idot);
		}
		delete [] v.dat;
		v.dat=ndat;
		return v;
	}

	template <class V,class MS,class S>
	TINLINE V &_vmscmultassign(V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(v.size!=ms.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS>(nameof(v)+" &operator *=("+nameof(v)+" &,const "+nameof(ms)+" &)"));
#endif
		v.size=ms.sxsize;
		S *ndat=new S[v.size];
		v.l=ms.start2;
		v.u=ms.end2;
		cdotprecision cdot(0.0);
		cdot.set_k(opdotprec);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(cdot)
#endif		
		for(int i=0;i<ms.sxsize;i++)
		{
			cdot=0.0;
			accumulate_approx(cdot,v,ms[Col(i+Lb(ms,COL))]);
			ndat[i]=rnd(cdot);
		}
		delete [] v.dat;
		v.dat=ndat;
		return v;
	}

	template <class V,class MS,class S>
	TINLINE V &_vmscimultassign(V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(v.size!=ms.sysize) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MS>(nameof(v)+" &operator *=("+nameof(v)+" &,const "+nameof(ms)+" &)"));
#endif
		v.size=ms.sxsize;
		S *ndat=new S[v.size];
		v.l=ms.start2;
		v.u=ms.end2;
		cidotprecision cidot(0.0);
		cidot.set_k(opdotprec);
#if defined(_OPENMP) && defined(CXSC_USE_OPENMP)
#pragma omp parallel for firstprivate(cidot)
#endif		
		for(int i=0;i<ms.sxsize;i++)
		{
			cidot=0.0;
			accumulate(cidot,v,ms[Col(i+Lb(ms,COL))]);
			ndat[i]=rnd(cidot);
		}
		delete [] v.dat;
		v.dat=ndat;
		return v;
	}

	template <class M,class S,class E>
	TINLINE E _msdiv(const M &m,const S &c) throw()
	{
		E r(m.lb1,m.ub1,m.lb2,m.ub2);

		for(int i=0;i<m.xsize*m.ysize;i++)
			r.dat[i]=m.dat[i]/c;
		return r;
	}

	template <class M,class S>
	TINLINE M &_msdivassign(M &m,const S &c) throw()
	{
		for(int i=0;i<m.xsize*m.ysize;i++)
			m.dat[i]/=c;
		return m;
	}

	template <class MS,class S,class E>
	TINLINE E _mssdiv(const MS &ms, const S &c) throw()
	{
		E r(ms.start1,ms.end1,ms.start2,ms.end2);

		for(int i=0;i<ms.sxsize;i++)
		{
			for(int j=0;j<ms.sysize;j++)
				r.dat[i*ms.sxsize+j]=ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2]/c;
		}
		return r;
	}

	template <class MS,class S>
	TINLINE MS &_mssdivassign(MS &ms,const S &c) throw()
	{
		for(int i=0;i<ms.sxsize;i++)
		{
			for(int j=0;j<ms.sysize;j++)
				ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2]/=c;
		}
		return ms;
	}
	
	template <class MS1,class MS2,class E>
	TINLINE E _msmsconv(const MS1 &m1,const MS2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.sxsize!=m2.sxsize)||(m1.sysize!=m2.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(E())+" operator |(const "+nameof(m1)+" &,const "+nameof(m2)+" &)"));
#endif
		int x=1,y=1;

		if((m1.start2==m2.start2)&&(m1.start1==m2.start1)) { y=m1.start1; x=m1.start2; }
		
		E r(y,y+m1.sysize-1,x,x+m1.sxsize-1);

		// no bounds checking here
		for(int i=0;i<m1.sysize;i++)
		{
			for(int j=0;j<m1.sxsize;j++)
				r.dat[i*m1.sxsize+j]=m1.dat[(i+m1.offset1)*m1.mxsize+m1.offset2+j]|m2.dat[(m2.offset1+i)*m2.mxsize+m2.offset2+j];
		}

		return r;
	}

	template <class MS1,class MS2,class E>
	TINLINE E _msmssect(const MS1 &m1,const MS2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.sxsize!=m2.sxsize)||(m1.sysize!=m2.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(E())+" operator &(const "+nameof(m1)+" &,const "+nameof(m2)+" &)"));
#endif
		int x=1,y=1;

		if((m1.start2==m2.start2)&&(m1.start1==m2.start1)) { y=m1.start1; x=m1.start2; }
		
		E r(y,y+m1.sysize-1,x,x+m1.sxsize-1);

		// no bounds checking here
		for(int i=0;i<m1.sysize;i++)
		{
			for(int j=0;j<m1.sxsize;j++)
				r.dat[i*m1.sxsize+j]=m1.dat[(i+m1.offset1)*m1.mxsize+m1.offset2+j]&m2.dat[(m2.offset1+i)*m2.mxsize+m2.offset2+j];
		}

		return r;
	}

	template <class MS1,class MS2,class E>
	TINLINE E _msmsplus(const MS1 &m1,const MS2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1.sxsize!=m2.sxsize)||(m1.sysize!=m2.sysize)) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(E())+" operator +(const "+nameof(m1)+" &,const "+nameof(m2)+" &)"));
#endif
		int x=1,y=1;

		if((m1.start2==m2.start2)&&(m1.start1==m2.start1)) { y=m1.start1; x=m1.start2; }
		
		E r(y,y+m1.sysize-1,x,x+m1.sxsize-1);

		// no bounds checking here
		for(int i=0;i<m1.sysize;i++)
		{
			for(int j=0;j<m1.sxsize;j++)
				r.dat[i*m1.sxsize+j]=m1.dat[(i+m1.offset1)*m1.mxsize+m1.offset2+j]+m2.dat[(m2.offset1+i)*m2.mxsize+m2.offset2+j];
		}

		return r;
	}

	template <class M>
	TINLINE void *_mvoid(const M &m) throw()
	{
		for(int i=0;i<m.xsize*m.ysize;i++)
		{
			if(!!m.dat[i])
				return (void *)1;
		}
		return (void *)0;
	}

	template <class M>
	TINLINE bool _mnot(const M &m) throw()
	{
		for(int i=0;i<m.xsize*m.ysize;i++)
		{
			if(!!m.dat[i])
				return false;
		}
		return true;
	}

	template <class MS>
	TINLINE void *_msvoid(const MS &ms) throw()
	{
		for(int i=0;i<ms.sysize;i++)
		{
			for(int j=0;j<ms.sysize;j++)
			{
				if(!!ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2])
					return (void *)1;
			}
		}
		return (void *)0;
	}

	template <class MS>
	TINLINE bool _msnot(const MS &ms) throw()
	{
		for(int i=0;i<ms.sysize;i++)
		{
			for(int j=0;j<ms.sxsize;j++)
			{
				if(!!ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2])
					return false;
			}
		}
		return true;
	}

	template <class M1,class M2>
	TINLINE bool _mmeq(const M1 &m1,const M2 &m2) throw()
	{
		if((m1.xsize!=m2.xsize)||(m1.ysize!=m2.ysize)) return false;
		for(int i=0;i<m1.xsize*m1.ysize;i++)
		{
			if(m1.dat[i]!=m2.dat[i])
				return false;
		}
		return true;
	}

	template <class M1,class M2>
	TINLINE bool _mmneq(const M1 &m1,const M2 &m2) throw()
	{
		if((m1.xsize!=m2.xsize)||(m1.ysize!=m2.ysize)) return true;
		for(int i=0;i<m1.xsize*m1.ysize;i++)
		{
			if(m1.dat[i]!=m2.dat[i])
				return true;
		}
		return false;
	}

	template <class M1,class M2>
	TINLINE bool _mmless(const M1 &m1,const M2 &m2) throw()
	{
		if((m1.xsize!=m2.xsize)||(m1.ysize!=m2.ysize)) return false;
		for(int i=0;i<m1.xsize*m1.ysize;i++)
		{
			if(m1.dat[i]>=m2.dat[i])
				return false;
		}
		return true;
	}

	template <class M1,class M2>
	TINLINE bool _mmleq(const M1 &m1,const M2 &m2) throw()
	{
		if((m1.xsize!=m2.xsize)||(m1.ysize!=m2.ysize)) return false;
		for(int i=0;i<m1.xsize*m1.ysize;i++)
		{
			if(m1.dat[i]>m2.dat[i])
				return false;
		}
		return true;
	}

	template <class M,class MS>
	TINLINE bool _mmseq(const M &m1,const MS &ms) throw()
	{
		if((m1.xsize!=ms.sxsize)||(m1.ysize!=ms.sysize)) return false;
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<m1.xsize;j++)
			{
				if(m1.dat[i*m1.xsize+j]!=ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2])
					return false;
			}
		}
		return true;
	}

	template <class M,class MS>
	TINLINE bool _mmsneq(const M &m1,const MS &ms) throw()
	{
		if((m1.xsize!=ms.sxsize)||(m1.ysize!=ms.sysize)) return true;
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<m1.xsize;j++)
			{
				if(m1.dat[i*m1.xsize+j]!=ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2])
					return true;
			}
		}
		return false;
	}

	template <class M,class MS>
	TINLINE bool _mmsless(const M &m1,const MS &ms) throw()
	{
		if((m1.xsize!=ms.sxsize)||(m1.ysize!=ms.sysize)) return false;
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<m1.xsize;j++)
			{
				if(m1.dat[i*m1.xsize+j]>=ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2])
					return false;
			}
		}
		return true;
	}

	template <class M,class MS>
	TINLINE bool _mmsleq(const M &m1,const MS &ms) throw()
	{
		if((m1.xsize!=ms.sxsize)||(m1.ysize!=ms.sysize)) return false;
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<m1.xsize;j++)
			{
				if(m1.dat[i*m1.xsize+j]>ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2])
					return false;
			}
		}
		return true;
	}

	template <class MS,class M>
	TINLINE bool _msmless(const MS &ms,const M &m1) throw()
	{
		if((m1.xsize!=ms.sxsize)||(m1.ysize!=ms.sysize)) return false;
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<m1.xsize;j++)
			{
				if(m1.dat[i*m1.xsize+j]<=ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2])
					return false;
			}
		}
		return true;
	}

	template <class MS,class M>
	TINLINE bool _msmleq(const MS &ms,const M &m1) throw()
	{
		if((m1.xsize!=ms.sxsize)||(m1.ysize!=ms.sysize)) return false;
		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<m1.xsize;j++)
			{
				if(m1.dat[i*m1.xsize+j]<ms.dat[(i+ms.offset1)*ms.mxsize+j+ms.offset2])
					return false;
			}
		}
		return true;
	}

	template <class MS1,class MS2>
	TINLINE bool _msmseq(const MS1 &ms1,const MS2 &ms2) throw()
	{
		if((ms1.sxsize!=ms2.sxsize)||(ms1.sysize!=ms2.sysize)) return false;
		for(int i=0;i<ms1.sysize;i++)
		{
			for(int j=0;j<ms1.sxsize;j++)
			{
				if(ms1.dat[(i+ms1.offset1)*ms1.mxsize+j+ms1.offset2]!=ms2.dat[(i+ms2.offset1)*ms2.mxsize+j+ms2.offset2])
					return false;
			}
		}
		return true;
	}

	template <class MS1,class MS2>
	TINLINE bool _msmsneq(const MS1 &ms1,const MS2 &ms2) throw()
	{
		if((ms1.sxsize!=ms2.sxsize)||(ms1.sysize!=ms2.sysize)) return true;
		for(int i=0;i<ms1.sysize;i++)
		{
			for(int j=0;j<ms1.sxsize;j++)
			{
				if(ms1.dat[(i+ms1.offset1)*ms1.mxsize+j+ms1.offset2]!=ms2.dat[(i+ms2.offset1)*ms2.mxsize+j+ms2.offset2])
					return true;
			}
		}
		return false;
	}

	template <class MS1,class MS2>
	TINLINE bool _msmsless(const MS1 &ms1,const MS2 &ms2) throw()
	{
		if((ms1.sxsize!=ms2.sxsize)||(ms1.sysize!=ms2.sysize)) return false;
		for(int i=0;i<ms1.sysize;i++)
		{
			for(int j=0;j<ms1.sxsize;j++)
			{
				if(ms1.dat[(i+ms1.offset1)*ms1.mxsize+j+ms1.offset2]>=ms2.dat[(i+ms2.offset1)*ms2.mxsize+j+ms2.offset2])
					return false;
			}
		}
		return true;
	}

	template <class MS1,class MS2>
	TINLINE bool _msmsleq(const MS1 &ms1,const MS2 &ms2) throw()
	{
		if((ms1.sxsize!=ms2.sxsize)||(ms1.sysize!=ms2.sysize)) return false;
		for(int i=0;i<ms1.sysize;i++)
		{
			for(int j=0;j<ms1.sxsize;j++)
			{
				if(ms1.dat[(i+ms1.offset1)*ms1.mxsize+j+ms1.offset2]>ms2.dat[(i+ms2.offset1)*ms2.mxsize+j+ms2.offset2])
					return false;
			}
		}
		return true;
	}


//------------------- matrix_subv -----------------------------


template <class V,class MV2,class S>
TINLINE V &_vmvassign(V &v,const MV2 &rv) throw()
{
	delete [] v.dat;
	v.dat=new S[rv.size];
	for(int i=0,j=rv.start;i<rv.size;i++,j+=rv.offset)
		v.dat[i]=S(rv.dat[j]);
	v.l=rv.lb;
	v.u=rv.ub;
	v.size=rv.size;
	return v;
}

template <class MV1,class MV2>
TINLINE MV1 &_mvmvassign(MV1 &v,const MV2 &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV1>)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(v.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MV1>(nameof(v)+" &"+nameof(v)+"::operator =(const "+nameof(rv)+" &)"));
#endif
	for(int i=0,j=v.start,k=rv.start;i<v.size;i++,j+=v.offset,k+=rv.offset)
		v.dat[j]=rv.dat[k];
	return v;
}

template <class MV,class S>
TINLINE MV &_mvsassign(MV &v,const  S &r) throw()
{
	for(int i=v.start,j=0;j<v.size;j++,i+=v.offset)
		v.dat[i]=r;
	return v;
}

template <class MV,class V>
TINLINE MV &_mvvassign(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(v.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MV>(nameof(v)+" &"+nameof(v)+"::operator =(const "+nameof(rv)+" &)"));
#endif
	for(int i=0,j=v.start;i<v.size;i++,j+=v.offset)
		v.dat[j]=rv.dat[i];
	return v;
}

template <class MV,class V>
TINLINE V _mvabs(const MV &mv) throw()
{
	V v(mv.lb,mv.ub);
	for(int i=0,j=mv.start;i<mv.size;i++,j+=mv.offset)
		v.dat[i]=abs(mv.dat[j]);
	return v;
}

template <class MV,class V>
TINLINE V _mvdiam(const MV &mv) throw()
{
	V v(mv.lb,mv.ub);
	for(int i=0,j=mv.start;i<mv.size;i++,j+=mv.offset)
		v.dat[i]=diam(mv.dat[j]);
	return v;
}

template <class MV,class V>
TINLINE V _mvmid(const MV &mv) throw()
{
	V v(mv.lb,mv.ub);
	for(int i=0,j=mv.start;i<mv.size;i++,j+=mv.offset)
		v.dat[i]=mid(mv.dat[j]);
	return v;
}

template <class MV,class V>
TINLINE V _mvinf(const MV &mv) throw()
{
	V v(mv.lb,mv.ub);
	for(int i=0,j=mv.start;i<mv.size;i++,j+=mv.offset)
		v.dat[i]=Inf(mv.dat[j]);
	return v;
}

template <class MV,class V>
TINLINE V _mvsup(const MV &mv) throw()
{
	V v(mv.lb,mv.ub);
	for(int i=0,j=mv.start;i<mv.size;i++,j+=mv.offset)
		v.dat[i]=Sup(mv.dat[j]);
	return v;
}

template <class MV,class V>
TINLINE V _mvim(const MV &mv) throw()
{
	V v(mv.lb,mv.ub);
	for(int i=0,j=mv.start;i<mv.size;i++,j+=mv.offset)
		v.dat[i]=Im(mv.dat[j]);
	return v;
}

template <class MV,class V>
TINLINE V _mvre(const MV &mv) throw()
{
	V v(mv.lb,mv.ub);
	for(int i=0,j=mv.start;i<mv.size;i++,j+=mv.offset)
		v.dat[i]=Re(mv.dat[j]);
	return v;
}

template <class DP,class V,class SV>
	TINLINE void _vmvaccu(DP &dp, const V & rv1, const SV &rv2)
#if(CXSC_INDEX_CHECK)
		throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(OP_WITH_WRONG_DIM("void accumulate("+nameof(dp)+" &, const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		for(int i=0;i<rv1.size;i++)
			accumulate(dp,rv1.dat[i],rv2.dat[rv2.start+i*rv2.offset]);
	}

template <class DP,class MV1,class MV2>
	TINLINE void _mvmvaccu(DP &dp, const MV1 & rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(OP_WITH_WRONG_DIM("void accumulate("+nameof(dp)+" &, const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		for(int i=0;i<rv1.size;i++)
			accumulate(dp,rv1.dat[rv1.start+i*rv1.offset],rv2.dat[rv2.start+i*rv2.offset]);
	}

	template <class MV1,class MV2,class S>
	TINLINE S _mvmvmult(const MV1 & rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MV1>(nameof(S())+" operator *(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		dotprecision dot(0.0);
                dot.set_k(opdotprec);
         	accumulate_approx(dot,rv1,rv2);

		return rnd(dot);
	}

	template <class MV1,class MV2,class S>
	TINLINE S _mvmvimult(const MV1 & rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MV1>(nameof(S())+" operator *(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		idotprecision idot(0.0);
		idot.set_k(opdotprec);
		accumulate(idot,rv1,rv2);

		return rnd(idot);
	}

	template <class MV1,class MV2,class S>
	TINLINE S _mvmvlmult(const MV1 & rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MV1>(nameof(S())+" operator *(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		dotprecision dot(0.0);
		for(int i=0;i<rv1.size;i++)
			accumulate(dot,rv1.dat[rv1.start+i*rv1.offset],rv2.dat[rv2.start+i*rv2.offset]);

		return S(dot);
	}

	template <class MV1,class MV2,class S>
	TINLINE S _mvmvlimult(const MV1 & rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MV1>(nameof(S())+" operator *(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		idotprecision idot(0.0);
		for(int i=0;i<rv1.size;i++)
			accumulate(idot,rv1.dat[rv1.start+i*rv1.offset],rv2.dat[rv2.start+i*rv2.offset]);

		return S(idot);
	}

	template <class MV1,class MV2,class S>
	TINLINE S _mvmvcmult(const MV1 & rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MV1>(nameof(S())+" operator *(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		cdotprecision cdot(0.0);
		cdot.set_k(opdotprec);
		accumulate_approx(cdot,rv1,rv2);

		return rnd(cdot);
	}

	template <class MV1,class MV2,class S>
	TINLINE S _mvmvcimult(const MV1 & rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MV1>(nameof(S())+" operator *(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		cidotprecision cidot(0.0);
		cidot.set_k(opdotprec);
		accumulate(cidot,rv1,rv2);

		return rnd(cidot);
	}

	template <class V,class MV,class S>
	TINLINE S _vmvmult(const V &rv1, const MV &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MV>(nameof(S())+" operator *(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		dotprecision dot(0.0);
                dot.set_k(opdotprec);
		accumulate_approx(dot,rv1,rv2);

		return rnd(dot);
	}

	template <class V,class MV,class S>
	TINLINE S _vmvimult(const V &rv1, const MV &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MV>(nameof(S())+" operator *(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		idotprecision idot(0.0);
		idot.set_k(opdotprec);
		accumulate(idot,rv1,rv2);

		return rnd(idot);
	}

	template <class V,class MV,class S>
	TINLINE S _vmvlmult(const V &rv1, const MV &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MV>(nameof(S())+" operator *(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		dotprecision dot(0.0);
		for(int i=0;i<rv1.size;i++)
			accumulate(dot,rv1.dat[i],rv2.dat[rv2.start+i*rv2.offset]);

		return S(dot);
	}

	template <class V,class MV,class S>
	TINLINE S _vmvlimult(const V &rv1, const MV &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MV>(nameof(S())+" operator *(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		idotprecision idot(0.0);
		for(int i=0;i<rv1.size;i++)
			accumulate(idot,rv1.dat[i],rv2.dat[rv2.start+i*rv2.offset]);

		return S(idot);
	}

	template <class V,class MV,class S>
	TINLINE S _vmvcmult(const V &rv1, const MV &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MV>(nameof(S())+" operator *(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		cdotprecision cdot(0.0);
		cdot.set_k(opdotprec);
		accumulate_approx(cdot,rv1,rv2);

		return rnd(cdot);
	}

	template <class V,class MV,class S>
	TINLINE S _vmvcimult(const V &rv1, const MV &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MV>(nameof(S())+" operator *(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		cidotprecision cidot(0.0);
		cidot.set_k(opdotprec);
		accumulate(cidot,rv1,rv2);

		return rnd(cidot);
	}

	template <class MV,class S,class E>
	TINLINE E _mvsmult(const MV &rv, const S &s) throw()
	{
		E p(rv.lb,rv.ub);

		for(int i=0;i<rv.size;i++)
			p.dat[i]=rv.dat[rv.start+i*rv.offset]*s;
		
		return p;
	}


	template <class MV1,class MV2,class E>
	TINLINE E _mvmvconv(const MV1 &rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E sum(rv1.lb,rv1.ub);

#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(sum)+" operator |(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		for (int i=0;i<rv1.size;i++)
			sum.dat[i]=rv1.dat[rv1.start+i*rv1.offset]|rv2.dat[rv2.start+i*rv2.offset];
		return sum;
	}

	template <class MV,class V,class E>
	TINLINE E _mvvconv(const MV &rv1, const V &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E sum(rv1.lb,rv1.ub);

#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(sum)+" operator |(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		for (int i=0;i<rv1.size;i++)
			sum.dat[i]=rv1.dat[rv1.start+i*rv1.offset]|rv2.dat[i];
		return sum;
	}

	template <class MV1,class MV2,class E>
	TINLINE E _mvmvsect(const MV1 &rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E sum(rv1.lb,rv1.ub);

#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(sum)+" operator &(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		for (int i=0;i<rv1.size;i++)
			sum.dat[i]=rv1.dat[rv1.start+i*rv1.offset]&rv2.dat[rv2.start+i*rv2.offset];
		return sum;
	}

	template <class MV,class V,class E>
	TINLINE E _mvvsect(const MV &rv1, const V &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E sum(rv1.lb,rv1.ub);

#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(sum)+" operator &(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		for (int i=0;i<rv1.size;i++)
			sum.dat[i]=rv1.dat[rv1.start+i*rv1.offset]&rv2.dat[i];
		return sum;
	}

	template <class MV1,class MV2,class E>
	TINLINE E _mvmvplus(const MV1 &rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E sum(rv1.lb,rv1.ub);

#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(sum)+" operator +(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		for (int i=0;i<rv1.size;i++)
			sum.dat[i]=rv1.dat[rv1.start+i*rv1.offset]+rv2.dat[rv2.start+i*rv2.offset];
		return sum;
	}

	template <class MV1,class MV2,class E>
	TINLINE E _mvmvminus(const MV1 &rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E sum(rv1.lb,rv1.ub);

#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(sum)+" operator -(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		for (int i=0;i<rv1.size;i++)
			sum.dat[i]=rv1.dat[rv1.start+i*rv1.offset]-rv2.dat[rv2.start+i*rv2.offset];
		return sum;
	}

	template <class MV,class V,class E>
	TINLINE E _mvvplus(const MV &rv1, const V &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E sum(rv1.lb,rv1.ub);

#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(sum)+" operator +(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		for (int i=0;i<rv1.size;i++)
			sum.dat[i]=rv1.dat[rv1.start+i*rv1.offset]+rv2.dat[i];
		return sum;
	}

	template <class MV,class V,class E>
	TINLINE E _mvvminus(const MV &rv1, const V &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E sum(rv1.lb,rv1.ub);

#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(sum)+" operator -(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		for (int i=0;i<rv1.size;i++)
			sum.dat[i]=rv1.dat[rv1.start+i*rv1.offset]-rv2.dat[i];
		return sum;
	}

	template <class V,class MV,class E>
	TINLINE E _vmvminus(const V &rv1, const MV &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E sum(rv1.l,rv1.u);

#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(sum)+" operator -(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		for (int i=0;i<rv1.size;i++)
			sum.dat[i]=rv1.dat[i]-rv2.dat[rv2.start+i*rv2.offset];
		return sum;
	}

	template <class MV,class S,class E>
	TINLINE E _mvsdiv(const MV &rv, const S &s) throw()
	{
		E p(rv.lb,rv.ub);

		for(int i=0;i<rv.size;i++)
			p.dat[i]=rv.dat[rv.start+i*rv.offset]/s;
		
		return p;
	}

template <class MV,class S>
TINLINE MV &_mvsmultassign(MV &v,const S &r) throw()
{
	for(int i=v.start,j=0;j<v.size;j++,i+=v.offset)
		v.dat[i]*=r;
	return v;
}

template <class MV, class S>
TINLINE MV &_mvsplusassign(MV &v,const S &r) throw()
{
	for(int i=v.start,j=0;j<v.size;j++,i+=v.offset)
		v.dat[i]+=r;
	return v;
}
	
template <class MV,class S>
TINLINE MV &_mvsminusassign(MV &v,const S &r) throw()
{
	for(int i=v.start,j=0;j<v.size;j++,i+=v.offset)
		v.dat[i]-=r;
	return v;
}
	
template <class MV,class S>
TINLINE MV &_mvsdivassign(MV &v,const S &r) throw()
{
	for(int i=v.start,j=0;j<v.size;j++,i+=v.offset)
		v.dat[i]/=r;
	return v;
}
	
template <class MV,class V>
TINLINE MV &_mvvconvassign(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(v.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MV>(nameof(v)+" &"+nameof(v)+"::operator |=(const "+nameof(rv)+" &)"));
#endif
	for(int i=0,j=v.start;i<v.size;i++,j+=v.offset)
		v.dat[j]|=rv.dat[i];
	return v;
}

template <class V,class MV>
TINLINE V &_vmvconvassign(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(v.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V>(nameof(rv)+" &"+nameof(rv)+"::operator |=(const "+nameof(v)+" &)"));
#endif
	for(int i=0,j=v.start;i<v.size;i++,j+=v.offset)
		rv.dat[i]|=v.dat[j];
	return rv;
}

template <class MV,class V>
TINLINE MV &_mvvsectassign(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(v.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MV>(nameof(v)+" &"+nameof(v)+"::operator &=(const "+nameof(rv)+" &)"));
#endif
	for(int i=0,j=v.start;i<v.size;i++,j+=v.offset)
		v.dat[j]&=rv.dat[i];
	return v;
}

template <class V,class MV>
TINLINE V &_vmvsectassign(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(v.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V>(nameof(rv)+" &"+nameof(rv)+"::operator &=(const "+nameof(v)+" &)"));
#endif
	for(int i=0,j=v.start;i<v.size;i++,j+=v.offset)
		rv.dat[i]&=v.dat[j];
	return rv;
}

template <class MV,class V>
TINLINE MV &_mvvplusassign(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(v.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MV>(nameof(v)+" &"+nameof(v)+"::operator +=(const "+nameof(rv)+" &)"));
#endif
	for(int i=0,j=v.start;i<v.size;i++,j+=v.offset)
		v.dat[j]+=rv.dat[i];
	return v;
}

template <class V,class MV>
TINLINE V &_vmvplusassign(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(v.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V>(nameof(rv)+" &"+nameof(rv)+"::operator +=(const "+nameof(v)+" &)"));
#endif
	for(int i=0,j=v.start;i<v.size;i++,j+=v.offset)
		rv.dat[i]+=v.dat[j];
	return rv;
}

template <class MV,class V>
TINLINE MV &_mvvminusassign(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(v.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MV>(nameof(v)+" &"+nameof(v)+"::operator -=(const "+nameof(rv)+" &)"));
#endif
	for(int i=0,j=v.start;i<v.size;i++,j+=v.offset)
		v.dat[j]-=rv.dat[i];
	return v;
}

template <class V,class MV>
TINLINE V &_vmvminusassign(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(v.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V>(nameof(rv)+" &"+nameof(rv)+"::operator -=(const "+nameof(v)+" &)"));
#endif
	for(int i=0,j=v.start;i<v.size;i++,j+=v.offset)
		rv.dat[i]-=v.dat[j];
	return rv;
}

	template <class MV,class S>
	TINLINE MV &_mvsusetsup(MV &mv, const S &s) throw()
	{
		for(int i=mv.start,j=0;j<mv.size;j++,i+=mv.offset)
			UncheckedSetInf(mv.dat[i],s);
		return mv;
	}

	template <class MV,class S>
	TINLINE MV &_mvsusetinf(MV &mv, const S &s) throw()
	{
		for(int i=mv.start,j=0;j<mv.size;j++,i+=mv.offset)
			UncheckedSetInf(mv.dat[i],s);
		return mv;
	}

	template <class MV,class S>
	TINLINE MV &_mvssetinf(MV &mv, const S &s) throw()
	{
		for(int i=mv.start,j=0;j<mv.size;j++,i+=mv.offset)
			SetInf(mv.dat[i],s);
		return mv;
	}

	template <class MV,class S>
	TINLINE MV &_mvssetsup(MV &mv, const S &s) throw()
	{
		for(int i=mv.start,j=0;j<mv.size;j++,i+=mv.offset)
			SetSup(mv.dat[i],s);
		return mv;
	}

	template <class MV,class S>
	TINLINE MV &_mvssetre(MV &mv, const S &s) throw()
	{
		for(int i=mv.start,j=0;j<mv.size;j++,i+=mv.offset)
			SetRe(mv.dat[i],s);
		return mv;
	}

	template <class MV,class S>
	TINLINE MV &_mvssetim(MV &mv, const S &s) throw()
	{
		for(int i=mv.start,j=0;j<mv.size;j++,i+=mv.offset)
			SetIm(mv.dat[i],s);
		return mv;
	}

template <class MV,class V>
TINLINE MV &_mvvsetinf(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(v.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MV>(nameof(v)+" &SetInf("+nameof(v)+" &,const "+nameof(rv)+" &)"));
#endif
	for(int i=0,j=v.start;i<v.size;i++,j+=v.offset)
		SetInf(v.dat[j],rv.dat[i]);
	return v;
}

template <class V,class MV>
TINLINE V &_vmvsetinf(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(v.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V>(nameof(rv)+" &SetInf("+nameof(rv)+" &,const "+nameof(v)+" &)"));
#endif
	for(int i=0,j=v.start;i<v.size;i++,j+=v.offset)
		SetInf(rv.dat[i],v.dat[j]);
	return rv;
}

template <class MV,class V>
TINLINE MV &_mvvusetinf(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(v.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MV>(nameof(v)+" &UncheckedSetInf("+nameof(v)+" &,const "+nameof(rv)+" &)"));
#endif
	for(int i=0,j=v.start;i<v.size;i++,j+=v.offset)
		UncheckedSetInf(v.dat[j],rv.dat[i]);
	return v;
}

template <class V,class MV>
TINLINE V &_vmvusetinf(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(v.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V>(nameof(rv)+" &UncheckedSetInf("+nameof(rv)+" &,const "+nameof(v)+" &)"));
#endif
	for(int i=0,j=v.start;i<v.size;i++,j+=v.offset)
		UncheckedSetInf(rv.dat[i],v.dat[j]);
	return rv;
}

template <class MV,class V>
TINLINE MV &_mvvsetsup(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(v.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MV>(nameof(v)+" &SetSup("+nameof(v)+" &,const "+nameof(rv)+" &)"));
#endif
	for(int i=0,j=v.start;i<v.size;i++,j+=v.offset)
		SetSup(v.dat[j],rv.dat[i]);
	return v;
}

template <class V,class MV>
TINLINE V &_vmvsetsup(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(v.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V>(nameof(rv)+" &SetSup("+nameof(rv)+" &,const "+nameof(v)+" &)"));
#endif
	for(int i=0,j=v.start;i<v.size;i++,j+=v.offset)
		SetSup(rv.dat[i],v.dat[j]);
	return rv;
}

template <class MV,class V>
TINLINE MV &_mvvusetsup(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(v.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MV>(nameof(v)+" &UncheckedSetSup("+nameof(v)+" &,const "+nameof(rv)+" &)"));
#endif
	for(int i=0,j=v.start;i<v.size;i++,j+=v.offset)
		UncheckedSetSup(v.dat[j],rv.dat[i]);
	return v;
}

template <class V,class MV>
TINLINE V &_vmvusetsup(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(v.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V>(nameof(rv)+" &UncheckedSetSup("+nameof(rv)+" &,const "+nameof(v)+" &)"));
#endif
	for(int i=0,j=v.start;i<v.size;i++,j+=v.offset)
		UncheckedSetSup(rv.dat[i],v.dat[j]);
	return rv;
}

template <class MV,class V>
TINLINE MV &_mvvsetim(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(v.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MV>(nameof(v)+" &SetIm("+nameof(v)+" &,const "+nameof(rv)+" &)"));
#endif
	for(int i=0,j=v.start;i<v.size;i++,j+=v.offset)
		SetIm(v.dat[j],rv.dat[i]);
	return v;
}

template <class V,class MV>
TINLINE V &_vmvsetim(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(v.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V>(nameof(rv)+" &SetIm("+nameof(rv)+" &,const "+nameof(v)+" &)"));
#endif
	for(int i=0,j=v.start;i<v.size;i++,j+=v.offset)
		SetIm(rv.dat[i],v.dat[j]);
	return rv;
}

template <class MV,class V>
TINLINE MV &_mvvsetre(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(v.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<MV>(nameof(v)+" &SetRe("+nameof(v)+" &,const "+nameof(rv)+" &)"));
#endif
	for(int i=0,j=v.start;i<v.size;i++,j+=v.offset)
		SetRe(v.dat[j],rv.dat[i]);
	return v;
}

template <class V,class MV>
TINLINE V &_vmvsetre(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(v.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V>(nameof(rv)+" &SetRe("+nameof(rv)+" &,const "+nameof(v)+" &)"));
#endif
	for(int i=0,j=v.start;i<v.size;i++,j+=v.offset)
		SetRe(rv.dat[i],v.dat[j]);
	return rv;
}





/*
template <class S,class S,class U>
rmatrix<U> _ms_mult(const rmatrix<S> &m,const S &s,const U &u=0) throw()
{
	rmatrix<U> r(m.lb1,m.ub1,m.lb2,m.ub2);

	for(int i=0;i<m.ysize*m.xsize;i++)
		r.dat[i]=m.dat[i]*s;

	return r;
}

template <class S,class S,class R>
cxscvector<R> _mv_mult(const rmatrix<S> &m,const rvector &v,R) throw()
{
	cxscvector<R> r(m.lb1,m.ub1);

	for(int i=0;i<m.ysize;i++)
	{
		dot=0.0;
		for(int j=0;j<v.size;j++)
			accumulate(dot,m.dat[i*m.xsize+j],v.dat[j]);
		r.dat[i]=rnd(dot);
	}

	return r;
}

template <class S,class S, class R>
rvector _vm_mult(const rvector &v,const rmatrix<S> &m,R) throw()
{
	cxscvector<R> r(m.lb2,m.ub2);

	for(int i=0;i<m.xsize;i++)
	{
		dot=0.0;
		for(int j=0;j<v.size;j++)
			accumulate(dot,m.dat[j*m.xsize+i],v.dat[j]);
		r.dat[i]=rnd(dot);
	}

	return r;
}
*/

} // namespace cxsc

#endif // _CXSC_MATRIX_INL_INCLUDED
