/////////////////////////////////////////////////////////////////////////////
/// @file SmithFixtures.h
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-06-12
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.edu.pl/ for details.

#ifndef CAPD_FILE_SMITHFIXTURES_H
#define CAPD_FILE_SMITHFIXTURES_H

#include <capd/vectalg/Matrix.hpp>

#include <boost/mpl/list.hpp>

namespace capd
{
  namespace test
  {

    using capd::vectalg::Matrix;

    struct SmithFixture
    {
      typedef Matrix<long long, 0, 0> MatrixType;
      static const bool bigInt = false;

      std::vector<Matrix<long long, 0, 0> > matrices_matrix;
      std::vector<Matrix<long long, 0, 0> > matrices_smithMatrix;
    };

    struct ExampleWeb: SmithFixture
    {
      ExampleWeb();
    };

    struct RandomSage: SmithFixture
    {
      static const bool bigInt = true;
      RandomSage();
    };

    struct Torus0: SmithFixture
    {
      Torus0();
    };

    struct Torus1: SmithFixture
    {
      Torus1();
    };


    struct ChessboardComplex_4x4_0: SmithFixture
    {
      ChessboardComplex_4x4_0();
    };

    struct KleinBottle_1: SmithFixture
    {
      KleinBottle_1();
    };

    struct MatchingComplex_7_0: SmithFixture
    {
      MatchingComplex_7_0();
    };


    struct RealProjectiveSpace_3_0: SmithFixture
    {
      RealProjectiveSpace_3_0();
    };

    struct RandomComplex_10x4_0: SmithFixture
    {
      RandomComplex_10x4_0();
    };

    typedef boost::mpl::list<ExampleWeb, RandomSage, Torus0, Torus1, ChessboardComplex_4x4_0, KleinBottle_1, MatchingComplex_7_0,
			     RealProjectiveSpace_3_0, RandomComplex_10x4_0> SmithFixtures;

  }
}

#endif // CAPD_FILE_SMITHFIXTURES_H
