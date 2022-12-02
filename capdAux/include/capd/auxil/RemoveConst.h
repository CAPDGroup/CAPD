/////////////////////////////////////////////////////////////////////////////
/// @file RemoveConst.h
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2015-11-18
/////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2000-2015 by the CAPD Group.
//
// This file constitutes a part of the CAPD library (capdAux),
// distributed under the terms of the GNU General Public License.
// Consult http://capd.ii.uj.edu.pl and  http://redhom.ii.edu.pl/ for details.
/////////////////////////////////////////////////////////////////////////////

#ifndef CAPD_FILE_CAPDAUX_AUXIL_REMOVECONST_H
#define CAPD_FILE_CAPDAUX_AUXIL_REMOVECONST_H

// #include
#include <capd/auxil/Logger.h>

namespace capd
{
  namespace auxil
  {

    // this should be replace by std::remove_cv from C++11
    template<typename T>
    struct RemoveConst
    {
      typedef T type;
    };

    template<typename T>
    struct RemoveConst<const T>
    {
      typedef T type;
    };

  }
}

#endif // CAPD_FILE_CAPDAUX_AUXIL_REMOVECONST_H
