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

  inline MatrixContainer() : BaseContainerType() {}
  inline MatrixContainer(size_type _rows, size_type _cols) : BaseContainerType(_rows*_cols){}
  inline MatrixContainer(const MatrixContainer& mc) : BaseContainerType(mc) {}
  inline MatrixContainer(size_type _rows, size_type _cols,bool) : BaseContainerType(_rows*_cols,true){}
  inline MatrixContainer(const Dimension&) : BaseContainerType(1){}
  inline MatrixContainer(const Dimension&,bool) : BaseContainerType(1,true){}
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
    BaseContainerType::operator=(a);
    return *this;
  }

  using BaseContainerType::begin;
  using BaseContainerType::end;
  using BaseContainerType::size;
  static Dimension dimension() { return Dimension(rows,cols); }
  inline void resize(size_type r, size_type c)
  {
    if(r!=rows || c!=cols)
      throw std::range_error("Cannot resize MatrixContainer of constant dimensions");
  }
protected:

  using BaseContainerType::data;
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
    : BaseContainerType(1u), m_rows(1), m_cols(1)
  {}
  inline MatrixContainer(size_type _rows, size_type _cols)
    : BaseContainerType(_rows*_cols), m_rows(_rows), m_cols(_cols)
  {}
  inline MatrixContainer(const MatrixContainer& mc)
    : BaseContainerType(mc), m_rows(mc.m_rows), m_cols(mc.m_cols)
  {}
  inline MatrixContainer(size_type _rows, size_type _cols,bool)
    : BaseContainerType(_rows*_cols,true), m_rows(_rows), m_cols(_cols)
  {}
  inline MatrixContainer(const Dimension& d,bool)
    : BaseContainerType(d.first*d.second,true), m_rows(d.first), m_cols(d.second)
  {}
  inline MatrixContainer(const Dimension& d)
    : BaseContainerType(d.first*d.second), m_rows(d.first), m_cols(d.second)
  {}


  MatrixContainer(MatrixContainer&& v) noexcept
    : BaseContainerType(std::move(v)), m_rows(0), m_cols(0){
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
    std::swap(*static_cast<BaseContainerType*>(&A_m1),*static_cast<BaseContainerType*>(&A_m2));
    std::swap(A_m1.m_rows,A_m2.m_rows);
    std::swap(A_m1.m_cols,A_m2.m_cols);
  }

  inline size_type numberOfRows() const {return m_rows;}
  inline size_type numberOfColumns() const {return m_cols;}

  MatrixContainer& operator=(const MatrixContainer& a)
  {
    BaseContainerType::operator=(a);
    m_rows = a.m_rows;
    m_cols = a.m_cols;
    return *this;
  }

  using BaseContainerType::begin;
  using BaseContainerType::end;
  size_type size() const { return BaseContainerType::size(); }
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

  using BaseContainerType::resize;
  using BaseContainerType::data;
};

}} // namespace capd::vectalg

#endif // _CAPD_VECTALG_MATRIXCONTAINER_H_

/// @}
