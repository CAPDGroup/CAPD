/////////////////////////////////////////////////////////////////////////////
/// @file PARIInterface.hpp
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-06-17
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.edu.pl/ for details.

#ifndef CAPD_FILE_PARIINTERFACE_HPP
#define CAPD_FILE_PARIINTERFACE_HPP


#include <capd/matrixAlgorithms/PARIInterface.h>
#include <capd/vectalg/Matrix.hpp>
#include <capd/vectalg/Vector.hpp>

#include <capd/config-capdAlg.h>
#include <stdexcept>

#ifdef HAVE_PARI

#include <pari/pari.h>
#include "PARIConvert.h"


namespace capd
{
  namespace matrixAlgorithms
  {
    template<typename Scalar>
    capd::vectalg::Vector<Scalar, 0>
    PARIInterface::smithForm(const capd::vectalg::Matrix<Scalar, 0, 0>& capdMatrix)
    {
      capd::vectalg::Matrix<Scalar, 0, 0> empty;
      return smithForm(capdMatrix, empty, empty, empty, empty);
    }

    template<typename Scalar>
    capd::vectalg::Vector<Scalar, 0>
    PARIInterface::smithForm(const capd::vectalg::Matrix<Scalar, 0, 0>& capdMatrix,
                             capd::vectalg::Matrix<Scalar, 0, 0>& capdU, capd::vectalg::Matrix<Scalar, 0, 0>& capdV)
    {
      capd::vectalg::Matrix<Scalar, 0, 0> empty;
      return smithForm(capdMatrix,
                       empty,
                       capdU,
                       capdV,
                       empty);
    }

    template<typename Scalar>
    capd::vectalg::Vector<Scalar, 0>
    PARIInterface::smithForm(const capd::vectalg::Matrix<Scalar, 0, 0>& capdMatrix,
                             capd::vectalg::Matrix<Scalar, 0, 0>& capdQ,
                             capd::vectalg::Matrix<Scalar, 0, 0>& capdQinv,
                             capd::vectalg::Matrix<Scalar, 0, 0>& capdR,
                             capd::vectalg::Matrix<Scalar, 0, 0>& capdRinv)
    {
      GEN pariMatrix = PARIConvert<capd::vectalg::Matrix<Scalar, 0, 0> >::to(capdMatrix);
      GEN pariQinv, pariR;

      const bool computeQinv = (!capdQinv.empty() || !capdQ.empty());
      const bool computeR = (!capdRinv.empty() || !capdR.empty());

      GEN divisorsPARI = ZM_snfall_i(pariMatrix, (computeQinv ? &pariQinv : NULL),
                                     (computeR ? &pariR : NULL), 1);

      capd::vectalg::Vector<Scalar, 0> divisors = PARIConvert<capd::vectalg::Vector<Scalar, 0> >::from(divisorsPARI);


      if (!capdMatrix.empty()) { // change returned matrices only if smith form computed. Otherwise we return unexpected 0x0 matrices
        if (!capdQinv.empty()) {
          capdQinv = PARIConvert<capd::vectalg::Matrix<Scalar, 0, 0> >::from(pariQinv);
        }

        if (!capdR.empty()) {
          capdR = PARIConvert<capd::vectalg::Matrix<Scalar, 0, 0> >::from(pariR);
        }

        if (!capdQ.empty()) {
           GEN tmp = PARIConvert<long>::to((long)1);
           GEN pariQ = ZM_inv(pariQinv, &tmp);
          capdQ = PARIConvert<capd::vectalg::Matrix<Scalar, 0, 0> >::from(pariQ);
        }

        if (!capdRinv.empty()) {
           GEN tmp = PARIConvert<long>::to((long)1);
          GEN pariRinv = ZM_inv(pariR, &tmp);
          capdRinv = PARIConvert<capd::vectalg::Matrix<Scalar, 0, 0> >::from(pariRinv);
        }
      }

      return capd::vectalg::Vector<Scalar, 0>(std::find(divisors.rbegin(), divisors.rend(), 0).base(), divisors.end());
    }

  }
}

#else

namespace capd
{
  namespace matrixAlgorithms
  {
    template<typename Scalar>
    capd::vectalg::Vector<Scalar, 0>
    PARIInterface::smithForm(const capd::vectalg::Matrix<Scalar, 0, 0>& /*capdMatrix*/)
    {
      throw std::logic_error("Not implemented");
    }

    template<typename Scalar>
    capd::vectalg::Vector<Scalar, 0>
    PARIInterface::smithForm(const capd::vectalg::Matrix<Scalar, 0, 0>& /*capdMatrix*/,
                             capd::vectalg::Matrix<Scalar, 0, 0>& /*capdU*/, capd::vectalg::Matrix<Scalar, 0, 0>& /*capdV*/)
    {
      throw std::logic_error("Not implemented");
    }

    template<typename Scalar>
    capd::vectalg::Vector<Scalar, 0>
    PARIInterface::smithForm(const capd::vectalg::Matrix<Scalar, 0, 0>& /*capdMatrix */,
                             capd::vectalg::Matrix<Scalar, 0, 0>& /*capdQ */,
                             capd::vectalg::Matrix<Scalar, 0, 0>& /*capdQinv */,
                             capd::vectalg::Matrix<Scalar, 0, 0>& /*capdR */,
                             capd::vectalg::Matrix<Scalar, 0, 0>& /*capdRinv */)
    {
      throw std::logic_error("Not implemented");
    }

  }
}

#endif // HAVE_PARI

#define PARI_INTERFACE_SMITH_FORM_INSTANCE(type)                        \
  template                                                              \
  capd::vectalg::Vector<type, 0>                                        \
  PARIInterface::smithForm<type>(const capd::vectalg::Matrix<type, 0, 0>& capdMatrix); \
  template                                                              \
  capd::vectalg::Vector<type, 0>                                        \
  PARIInterface::smithForm(const capd::vectalg::Matrix<type, 0, 0>& capdMatrix, \
                           capd::vectalg::Matrix<type, 0, 0>& capdU, capd::vectalg::Matrix<type, 0, 0>& capdV); \
  template                                                              \
  capd::vectalg::Vector<type, 0>                                        \
  PARIInterface::smithForm(const capd::vectalg::Matrix<type, 0, 0>& capdMatrix, \
                           capd::vectalg::Matrix<type, 0, 0>& capdQ,    \
                           capd::vectalg::Matrix<type, 0, 0>& capdQinv, \
                           capd::vectalg::Matrix<type, 0, 0>& capdR,    \
                           capd::vectalg::Matrix<type, 0, 0>& capdRinv);

#endif // CAPD_FILE_PARIINTERFACE_HPP
