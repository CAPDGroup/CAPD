/////////////////////////////////////////////////////////////////////////////
/// @file SmithFormFactory.h
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-06-19
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.edu.pl/ for details.

#ifndef CAPD_FILE_SMITHFORMFACTORY_H
#define CAPD_FILE_SMITHFORMFACTORY_H

#include "SmithForm.h"

namespace capd
{

  namespace matrixAlgorithms
  {

    class SmithFormFactory
    {
    public:
      explicit SmithFormFactory(bool usePARI):
	_usePARI(usePARI)
      {}

      template<typename MatrixT>
      SmithForm<MatrixT>* operator()(MatrixT& B, bool computeQ, bool computeQinv, bool computeR, bool computeRinv);


    private:
      bool _usePARI;
    };

  }

}

#endif // CAPD_FILE_SMITHFORMFACTORY_H
