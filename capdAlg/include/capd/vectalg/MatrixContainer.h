/// @addtogroup vectalg
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file MatrixContainer.h
///
/// This file provides a class MatrixContainer
/// This class inherites form general Container class
/// and provides constructors and methods specific for two dimensional data
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_VECTALG_MATRIXCONTAINER_H_
#define _CAPD_VECTALG_MATRIXCONTAINER_H_

#include "capd/vectalg/Container.h"
#include <utility>
#include "capd/settings/compilerSetting.h"
#include "capd/auxil/Dll.h"

namespace capd{
namespace vectalg{

/// This class inherits form general Container class
/// and provides constructors and methods specific for two dimensional data
template<typename Scalar,__size_type rows, __size_type cols>
class MatrixContainer : public Container<Scalar,rows*cols>
{
  typedef Container<Scalar,rows*cols> BaseContainerType;
public:
  typedef Scalar ScalarType;
  typedef typename Container<Scalar,rows*cols>::iterator iterator;
  typedef typename Container<Scalar,rows*cols>::const_iterator const_iterator;
  typedef MatrixContainer ContainerType;
  typedef Container<Scalar,rows> ColumnContainer;
  typedef Container<Scalar,cols> RowContainer;
  typedef __size_type size_type;
  typedef __difference_type difference_type;
  typedef std::pair<size_type,size_type> Dimension;

  static const size_type ROWS = rows;
  static const size_type COLS = cols;

  inline MatrixContainer() : Container<Scalar,rows*cols>() {}
  inline MatrixContainer(size_type _rows, size_type _cols) : Container<Scalar,rows*cols>(_rows*_cols){}
  inline MatrixContainer(const MatrixContainer& mc) : Container<Scalar,rows*cols>(mc) {}
  inline MatrixContainer(size_type _rows, size_type _cols,bool) : Container<Scalar,rows*cols>(_rows*_cols,true){}
  inline MatrixContainer(const Dimension&) : Container<Scalar,rows*cols>(1){}
  inline MatrixContainer(const Dimension&,bool) : Container<Scalar,rows*cols>(1,true){}
  //  MatrixContainer(MatrixContainer&& v) : ContainerType(std::forward<ContainerType>(v)) {}
//  MatrixContainer & operator=(MatrixContainer && v) {
//     ContainerType::operator= ( std::forward<ContainerType>(v));
//   //  std::cout << "\n v move =";
//    return *this;
//  }

  inline size_type numberOfRows() const {return rows;}
  inline size_type numberOfColumns() const {return cols;}

  inline MatrixContainer& operator=(const MatrixContainer& a)
  {
    Container<Scalar,rows*cols>::operator=(a);
    return *this;
  }

  using Container<Scalar,rows*cols>::begin;
  using Container<Scalar,rows*cols>::end;
  using Container<Scalar,rows*cols>::size;
  static Dimension dimension() { return Dimension(rows,cols); }
  inline void resize(size_type r, size_type c)
  {
    if(r!=rows || c!=cols)
      throw std::range_error("Cannot resize MatrixContainer of constant dimensions");
  }
protected:

  using Container<Scalar,rows*cols>::data;
};

// ---------------------------------------------------------------------

template<typename Scalar>
class MatrixContainer<Scalar,0,0> : public Container<Scalar,0>
{
   typedef Container<Scalar,0> BaseContainerType;
public:
  typedef Scalar ScalarType;
  typedef typename Container<Scalar,0>::iterator iterator;
  typedef typename Container<Scalar,0>::const_iterator const_iterator;
  typedef MatrixContainer ContainerType;
  typedef Container<Scalar,0> ColumnContainer;
  typedef Container<Scalar,0> RowContainer;
  typedef __size_type size_type;
  typedef __difference_type difference_type;
  typedef std::pair<size_type,size_type> Dimension;
  static const size_type ROWS = 0;
  static const size_type COLS = 0;

  inline MatrixContainer()
    : Container<Scalar,0>(0),
      m_rows(0), m_cols(0)
  {}
  inline MatrixContainer(size_type _rows, size_type _cols)
    : Container<Scalar,0>(_rows*_cols), m_rows(_rows), m_cols(_cols)
  {}
  inline MatrixContainer(const MatrixContainer& mc)
    : Container<Scalar,0>(mc), m_rows(mc.m_rows), m_cols(mc.m_cols)
  {}
  inline MatrixContainer(size_type _rows, size_type _cols,bool)
    : Container<Scalar,0>(_rows*_cols,true), m_rows(_rows), m_cols(_cols)
  {}
  inline MatrixContainer(const Dimension& d,bool)
    : Container<Scalar,0>(d.first*d.second,true), m_rows(d.first), m_cols(d.second)
  {}
  inline MatrixContainer(const Dimension& d)
    : Container<Scalar,0>(d.first*d.second), m_rows(d.first), m_cols(d.second)
  {}


  MatrixContainer(MatrixContainer&& v) noexcept
    : BaseContainerType(std::move(v)){
    std::swap(m_rows, v.m_rows);
    std::swap(m_cols, v.m_cols);
  }
  MatrixContainer & operator=(MatrixContainer && v) noexcept {
    BaseContainerType::operator=(std::move(v));
    std::swap(m_rows, v.m_rows);
    std::swap(m_cols, v.m_cols);
    return *this;
  }

  friend void swap(MatrixContainer<Scalar,0,0>& A_m1, MatrixContainer<Scalar,0,0>& A_m2)
  {
    std::swap(*static_cast<Container<Scalar,0>*>(&A_m1),*static_cast<Container<Scalar,0>*>(&A_m2));
    std::swap(A_m1.m_rows,A_m2.m_rows);
    std::swap(A_m1.m_cols,A_m2.m_cols);
  }

  inline size_type numberOfRows() const {return m_rows;}
  inline size_type numberOfColumns() const {return m_cols;}

  MatrixContainer& operator=(const MatrixContainer& a)
  {
    Container<Scalar,0>::operator=(a);
    m_rows = a.m_rows;
    m_cols = a.m_cols;
    return *this;
  }

  using Container<Scalar,0>::begin;
  using Container<Scalar,0>::end;
  using Container<Scalar,0>::size;
  inline Dimension dimension() const {return Dimension(m_rows,m_cols);}
  inline void resize(size_type r, size_type c)
   {
     if(r!=m_rows || c!=m_cols)
     {
       Container<Scalar,0>::resize(r*c);
       m_rows = r;
       m_cols = c;
     }
   }
protected:
  size_type m_rows,m_cols;

  using Container<Scalar,0>::resize;
  using Container<Scalar,0>::data;
};

}} // namespace capd::vectalg

#endif // _CAPD_VECTALG_MATRIXCONTAINER_H_

/// @}
