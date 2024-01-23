/////////////////////////////////////////////////////////////////////////////
/// @file LinearEquationFixtures.h
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-06-18
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.edu.pl/ for details.

#ifndef CAPD_FILE_LINEAREQUATIONFIXTURES_H
#define CAPD_FILE_LINEAREQUATIONFIXTURES_H


#include <capd/vectalg/Matrix.hpp>

#include <boost/mpl/list.hpp>

namespace capd
{
  namespace test
  {

    namespace linearEquation
    {

      using capd::vectalg::Matrix;
      using capd::vectalg::Vector;

      struct LinearEquationFixture
      {
	typedef Matrix<long long, 0, 0> MatrixType;
	typedef Vector<long long, 0> VectorType;
	static const bool bigInt = false;

	std::vector<MatrixType> data_le_matrix;
	std::vector<VectorType> data_le_vector;
	std::vector<VectorType> data_le_solution;
      };

      struct Torus_0: LinearEquationFixture
      {
	Torus_0();
      };

      struct Torus_1: LinearEquationFixture
      {
	Torus_1();
      };

      struct ChessboardComplex_4x4_0: LinearEquationFixture
      {
	ChessboardComplex_4x4_0();
      };

      struct KleinBottle_1: LinearEquationFixture
      {
	KleinBottle_1();
      };

      struct MatchingComplex_7_0: LinearEquationFixture
      {
	MatchingComplex_7_0();
      };

      struct RealProjectiveSpace_3_0: LinearEquationFixture
      {
	RealProjectiveSpace_3_0();
      };

      struct RandomComplex_10x4_0: LinearEquationFixture
      {
	RandomComplex_10x4_0();
      };

      typedef boost::mpl::list<Torus_0, Torus_1, ChessboardComplex_4x4_0, KleinBottle_1, MatchingComplex_7_0,
			       RealProjectiveSpace_3_0, RandomComplex_10x4_0> LinearEquationFixtures;
    }
  }
}
#endif // CAPD_FILE_LINEAREQUATIONFIXTURES_H
