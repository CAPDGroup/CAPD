#ifndef CAPD_TEST_COMPARE_H
#define CAPD_TEST_COMPARE_H

#include "capd/map/lib.h"
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "capd/rounding/DoubleRounding.h"

typedef capd::IMap MapType;
typedef MapType::VectorType VectorType;
typedef MapType::ScalarType ScalarType;
typedef MapType::JetType JetType;

template <typename DVectorT, typename IVectorT>
inline void compareResults(const DVectorT & expected, const IVectorT & result, const std::string & msg, double tolerance = 1.3e-12)
{
  typename IVectorT::const_iterator it = result.begin();
  typename DVectorT::const_iterator ex = expected.begin();
  for(int i=0; ex != expected.end(); ++ex, ++i, ++it){
#ifdef __USE_TIGHT__
    // TIGHT intervals are correctly rounded, hence much narrower than filib ones, while the expected value
    // is a plain double evaluation of a formula (error of a few ulps). Allow a small relative slack.
    const double slack = 1e-14 * (1.0 + std::abs(static_cast<double>(*ex)));
    ScalarType enlarged(it->leftBound() - slack, it->rightBound() + slack);
    BOOST_CHECK_MESSAGE(subset(*ex, enlarged), msg << "  derivative : " << i << ", " << *ex << ", " << *it);
#else
    BOOST_CHECK_MESSAGE(subset(*ex,*it), msg << "  derivative : " << i << ", " << *ex << ", " << *it);
#endif
    BOOST_CHECK_SMALL(it->rightBound()-it->leftBound(), tolerance);
  }
}

// wrapper for Mathematica functions

inline double Power(double x, double n) { return pow(x,n); }
inline double Sin(double x) { return sin(x); }
inline double Cos(double x) { return cos(x); }
inline double ArcSin(double x) { return asin(x); }
inline double ArcCos(double x) { return acos(x); }
inline double ArcTan(double x) { return atan(x); }
inline double Sqrt(double x) { return sqrt(x); }
inline double Log(double x) { return log(x); }

#endif // CAPD_TEST_COMPARE_H
