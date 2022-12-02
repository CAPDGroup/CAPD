/////////////////////////////////////////////////////////////////////////////
/// @file PARISmithForm.hpp
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-06-20
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD::RedHom Group.
//
// This file constitutes a part of the CAPD::RedHom library,
// distributed under the terms of the GNU General Public License.
// Consult  http://redhom.ii.edu.pl/ for details.

#ifndef CAPD_FILE_PARISMITHFORM_HPP
#define CAPD_FILE_PARISMITHFORM_HPP

#include <capd/matrixAlgorithms/PARISmithForm.h>
#include <capd/matrixAlgorithms/PARIInterface.h>
#include <stdexcept>

namespace capd
{
  namespace matrixAlgorithms
  {
    template<typename Scalar>
    bool PARISmithForm<Scalar>::enabled()
    {
      return PARIInterface::enabled();
    }

    template<typename Scalar>
    void PARISmithForm<Scalar>::operator()()
    {
      if (!enabled()) {
        throw std::logic_error("PARISmith form is not enabled");
      }

      CAPD_DEBUG("SmithForm for matrix size: " << _m << "x" << _n);
      const bool trace = (_n <= 100 && _m <= 100);

      if (trace) {
        CAPD_TRACE("PARISmithForm args: " << cppReprezentation(_B, "B", "TYPE"));
      }

      PARIInterface& pari = PARIInterface::instance();

      const capd::vectalg::Vector<Scalar, 0> divisors =
        pari.smithForm(_B, _Q, _Qinv, _R, _Rinv);

      for (size_t i = 0, end = _m;
           i < end / 2; ++i) {
        capd::matrixAlgorithms::rowExchange(_Qinv, end - i, i + 1);
        capd::matrixAlgorithms::columnExchange(_Q, end - i, i + 1);
      }

      for (size_t i = 0, end = _n;
           i < end / 2; ++i) {
        capd::matrixAlgorithms::columnExchange(_R, end - i, i + 1);
        capd::matrixAlgorithms::rowExchange(_Rinv, end - i, i + 1);
      }

      _B = Scalar(0);
      _t = divisors.dimension();
      _s = 0;

      for (size_t i = 0, end = divisors.dimension(); i < end; ++i) {
        _B[i][i] = divisors[end - i - 1];
        if (_B[i][i] == Scalar(1)) {
          ++_s;
        }
      }

      if (trace) {
        CAPD_TRACE("PARISmithForm result: " << cppReprezentation(_B, "B", "TYPE"));
      }
    }

#define PARI_SMITH_FORM(type)                   \
    template                                    \
    void PARISmithForm<type>::operator()();     \
    template                                    \
    bool PARISmithForm<type>::enabled();

  }
}

#endif // CAPD_FILE_PARISMITHFORM_H
