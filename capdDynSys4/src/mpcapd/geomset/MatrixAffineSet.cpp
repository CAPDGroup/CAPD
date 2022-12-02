#include "capd/vectalg/mplib.h"
#include "capd/geomset/MatrixAffineSet.hpp"
#include "capd/vectalg/iobject.hpp"
#include "capd/vectalg/ColumnVector.hpp"
#include "capd/vectalg/RowVector.hpp"
#include "capd/vectalg/Matrix.hpp"

#ifdef __HAVE_MPFR__
namespace capd{
namespace geomset{
  template class MatrixAffineSet<capd::MpIMatrix >;
}}
#endif
