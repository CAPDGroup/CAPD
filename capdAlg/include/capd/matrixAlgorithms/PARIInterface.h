/////////////////////////////////////////////////////////////////////////////
/// @file PARIInterface.h
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-06-15
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.edu.pl/ for details.

#ifndef CAPD_FILE_PARIINTERFACE_H
#define CAPD_FILE_PARIINTERFACE_H

#include <capd/vectalg/Matrix.h>
#include <capd/vectalg/Vector.h>

#include <assert.h>

namespace capd
{

  namespace matrixAlgorithms
  {

    class PARIInterface
    {
      static const bool _enabled;
    public:

      virtual ~PARIInterface();
      static bool enabled() { return _enabled; }

      template<typename Scalar>
      capd::vectalg::Vector<Scalar, 0>
      smithForm(const capd::vectalg::Matrix<Scalar, 0, 0>& capdMatrix);

      template<typename Scalar>
      capd::vectalg::Vector<Scalar, 0>
      smithForm(const capd::vectalg::Matrix<Scalar, 0, 0>& capdMatrix,
		capd::vectalg::Matrix<Scalar, 0, 0>& capdU, capd::vectalg::Matrix<Scalar, 0, 0>& capdV);
      template<typename Scalar>
      capd::vectalg::Vector<Scalar, 0>
      smithForm(const capd::vectalg::Matrix<Scalar, 0, 0>& capdMatrix,
		  capd::vectalg::Matrix<Scalar, 0, 0>& capdQ,
		  capd::vectalg::Matrix<Scalar, 0, 0>& capdQInv,
		  capd::vectalg::Matrix<Scalar, 0, 0>& capdR,
		  capd::vectalg::Matrix<Scalar, 0, 0>& capdRInv);

      static PARIInterface& instance()
      {
        static PARIInterface instance;
        return instance;
      }

    private:
      PARIInterface();

    };
  }
}

#endif // CAPD_FILE_PARIINTERFACE_H
