/////////////////////////////////////////////////////////////////////////////
/// @file PARIInterface.cpp
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-06-17
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD::RedHom Group.
//
// This file constitutes a part of the CAPD::RedHom library,
// distributed under the terms of the GNU General Public License.
// Consult  http://redhom.ii.edu.pl/ for details.


#include "../../capd/matrixAlgorithms/PARIInterface.hpp"
#include <capd/multiPrec/MpInt.h>
#include <capd/config-capdAlg.h>

using namespace capd::matrixAlgorithms;

#if HAVE_MPCAPD_ALG

#ifdef HAVE_PARI

namespace
{
  template<>
  struct PARIConvert<capd::multiPrec::MpInt>
  {
    static GEN to(const capd::multiPrec::MpInt& value)
    {
      mpz_srcptr z = value.get_mpz_t();
      // Function from pari-gnump-0.0.1
      const long lz = z->_mp_size;
      const long lx = labs (lz);
      const long lx2 = lx + 2;
      int i;
      GEN x = cgeti (lx2);

      assert (sizeof (long) == sizeof (mp_limb_t));

      x [1] = evalsigne ((lz > 0 ? 1 : (lz < 0 ? -1 : 0))) | evallgefint (lx2);
      for (i = 0; i < lx; i++)
	*int_W (x, i) = (z->_mp_d) [i];

      return x;
    }

    static capd::multiPrec::MpInt from(const GEN& x)
    {
      capd::multiPrec::MpInt result;
      mpz_ptr z = result.get_mpz_t();

      const long lx = lgefint (x) - 2;
      const long sign = signe (x);
      int i;

      assert(sizeof (long) == sizeof (mp_limb_t));
      assert(typ (x) == t_INT);

      if (sign == 0)
	mpz_set_ui (z, 0);
      else {
	mpz_realloc2 (z, lx * BITS_IN_LONG);
	z->_mp_size = sign * lx;
	for (i = 0; i < lx; i++)
	  (z->_mp_d) [i] = *int_W (x, i);
      }

      return result;
    }
  };
}

#endif

    using capd::multiPrec::MpInt;

    PARI_INTERFACE_SMITH_FORM_INSTANCE(MpInt);
#endif
