
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#define BOOST_TEST_MODULE IntervalTest
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "capd/rounding/DoubleRounding.h"
#include "capd/intervals/DoubleInterval.h"
#include "capd/intervals/IntervalError.h"

using namespace std;

using interval = capd::interval;

BOOST_AUTO_TEST_CASE(basics)
{
    interval a {};
    BOOST_CHECK_EQUAL(a, interval(0.0, 0.0));
    BOOST_CHECK_EQUAL(a.leftBound(), 0.0);
    BOOST_CHECK_EQUAL(a.rightBound(), 0.0);
    BOOST_CHECK_EQUAL(a.left(), interval(0.0, 0.0));
    BOOST_CHECK_EQUAL(a.right(), interval(0.0, 0.0));

    interval b(1.0);
    BOOST_CHECK_EQUAL(b, interval(1.0, 1.0));
    BOOST_CHECK_EQUAL(b.leftBound(), 1.0);
    BOOST_CHECK_EQUAL(b.rightBound(), 1.0);
    BOOST_CHECK_EQUAL(b.left(), interval(1.0, 1.0));
    BOOST_CHECK_EQUAL(b.right(), interval(1.0, 1.0));

    interval c(2.0, 3.0);
    BOOST_CHECK_EQUAL(c, interval(2.0, 3.0));
    BOOST_CHECK_EQUAL(c.leftBound(), 2.0);
    BOOST_CHECK_EQUAL(c.rightBound(), 3.0);
    BOOST_CHECK_EQUAL(c.left(), interval(2.0, 2.0));
    BOOST_CHECK_EQUAL(c.right(), interval(3.0, 3.0));

    interval d(c);
    BOOST_CHECK_EQUAL(d, interval(2.0, 3.0));
    BOOST_CHECK_EQUAL(d.leftBound(), 2.0);
    BOOST_CHECK_EQUAL(d.rightBound(), 3.0);
    BOOST_CHECK_EQUAL(d.left(), interval(2.0, 2.0));
    BOOST_CHECK_EQUAL(d.right(), interval(3.0, 3.0));

    interval e(2.5, 3.0);
    BOOST_CHECK_EQUAL(e, interval(2.5, 3.0));
    BOOST_CHECK_EQUAL(e.leftBound(), 2.5);
    BOOST_CHECK_EQUAL(e.rightBound(), 3.0);
    BOOST_CHECK_EQUAL(e.left(), interval(2.5, 2.5));
    BOOST_CHECK_EQUAL(e.right(), interval(3.0, 3.0));

    BOOST_CHECK_EQUAL(left(c), 2.0);
    BOOST_CHECK_EQUAL(right(c), 3.0);

    BOOST_CHECK(a < b);
    BOOST_CHECK(a <= b);
    BOOST_CHECK(!(a > b));
    BOOST_CHECK(!(a >= b));

    BOOST_CHECK(!(b==c));
    BOOST_CHECK(b!=c);

    BOOST_CHECK(b < c);
    BOOST_CHECK(b <= c);
    BOOST_CHECK(!(b > c));
    BOOST_CHECK(!(b >= c));

    BOOST_CHECK(b == 1.0);
    BOOST_CHECK(!(b != 1.0));

    BOOST_CHECK(!(b < 1.0));
    BOOST_CHECK(b <= 1.0);
    BOOST_CHECK(!(b > 1.0));
    BOOST_CHECK(b >= 1.0);

    BOOST_CHECK(c.contains(e));
    BOOST_CHECK(!(c.containsInInterior(e)));
    BOOST_CHECK(c.contains(2.5));
    BOOST_CHECK(d.subset(c));
    BOOST_CHECK(!(d.subsetInterior(c)));

    c.split(c,b);
    BOOST_CHECK_EQUAL(c, interval(2.5, 2.5));
    BOOST_CHECK_EQUAL(b, interval(-0.5, 0.5));
}

BOOST_AUTO_TEST_CASE(parse_from_sstream)
{
    interval a {};
    std::istringstream myStr("[3.21312312, 4.324324324]");
    myStr >> a;
    BOOST_CHECK_SMALL(a.leftBound() - 3.2131312, 8.1e-6);
    BOOST_CHECK_CLOSE(a.rightBound(), 4.324324324, 0.0);
}
 
BOOST_AUTO_TEST_CASE(operators)
{
    interval c(2.0, 3.0);
    interval d = c;
    interval e(-1.0, 4.0);
    interval f(4.0, 5.0);
    
    c += e;
    BOOST_CHECK_EQUAL(c, interval(1.0, 7.0));
    
    c -= e;
    BOOST_CHECK_EQUAL(c, interval(-3.0, 8.0));
    
    c *= e;
    BOOST_CHECK_EQUAL(c, interval(-12.0, 32.0));
    
    c /= f;
    BOOST_CHECK_EQUAL(c, interval(-3.0, 8.0));
    
    interval b(1.0);
    b = (b*(c+(-d))/f - e);
    BOOST_CHECK_EQUAL(b, interval(-5.5, 2.5));

    const interval y = interval(-3, -2);

    BOOST_CHECK_EQUAL(y+3, interval(0.0, 1.0));
    BOOST_CHECK_EQUAL(y-3, interval(-6.0, -5.0));
    BOOST_CHECK_EQUAL(y*3, interval(-9.0, -6.0));
    BOOST_CHECK_EQUAL(y/3, interval(-1.0, -2.0/3));
    BOOST_CHECK_EQUAL(3+y, interval(0.0, 1.0));
    BOOST_CHECK_EQUAL(3-y, interval(5.0, 6.0));
    BOOST_CHECK_EQUAL(3*y, interval(-9.0, -6.0));
    BOOST_CHECK_EQUAL(3/y, interval(-1.5, -1.0));
}

BOOST_AUTO_TEST_CASE(functions)
{
    interval e(-1.0, 4.0);

    constexpr double pi    = 3.14159265358979323;
    interval interval_pi = interval::pi();
    BOOST_CHECK(interval_pi.contains(pi));
    BOOST_CHECK_LE(interval_pi.rightBound() - interval_pi.leftBound(), 5.0e-16);

    constexpr double sqrt2 = 1.41421356237309504;
    interval interval_sqrt2 = sqrt(interval(2.0));
    BOOST_CHECK(interval_sqrt2.contains(sqrt2));
    BOOST_CHECK_LE(interval_sqrt2.rightBound() - interval_sqrt2.leftBound(), 3.0e-16);

    constexpr double sqrt3 = 1.73205080756887729;
    interval interval_sqrt3 = sqrt(interval(3.0));
    BOOST_CHECK(interval_sqrt3.contains(sqrt3));
    BOOST_CHECK_LE(interval_sqrt3.rightBound() - interval_sqrt3.leftBound(), 3.0e-16);

    interval x((interval::pi() / 6).leftBound(),
                (interval::pi() / 4).rightBound());

    auto sin_x = sin(x);
    BOOST_CHECK_SMALL(sin_x.leftBound() - 0.5, 4.5e-16);
    BOOST_CHECK_SMALL(sin_x.rightBound() - sqrt2 / 2, 2.3e-16);

    auto cos_x = cos(x);
    BOOST_CHECK_SMALL(cos_x.leftBound() - sqrt2 / 2, 3.4e-16);
    BOOST_CHECK_SMALL(cos_x.rightBound() - sqrt3 / 2, 4.5e-16);

    auto tan_x = tan(x);
    BOOST_CHECK_SMALL(tan_x.leftBound() - 1 / sqrt3, 8.9e-16);
    BOOST_CHECK_SMALL(tan_x.rightBound() - 1, 8.9e-16);

    auto cot_x = cot(x);
    BOOST_CHECK_SMALL(cot_x.leftBound() - 1, 8.9e-16);
    BOOST_CHECK_SMALL(cot_x.rightBound() - sqrt3, 2.5e-15);

    x = interval(1.0/3, 1.0/2);

    constexpr double atan_1_3 = 0.32175055439664219;
    constexpr double atan_1_2 = 0.46364760900080612;

    auto atan_x = atan(x);
    BOOST_CHECK_SMALL(atan_x.leftBound() - atan_1_3, 0.0);
    BOOST_CHECK_SMALL(atan_x.rightBound() - atan_1_2, 3.4e-16);

    constexpr double asin_1_3 = 0.33983690945412194;

    auto asin_x = asin(x);
    BOOST_CHECK_SMALL(asin_x.leftBound() - asin_1_3, 5.6e-17);
    BOOST_CHECK_SMALL(asin_x.rightBound() - pi / 6, 6.7e-16);

    constexpr double acos_1_3 = 1.23095941734077468;

    auto acos_x = acos(x);
    BOOST_CHECK_SMALL(acos_x.leftBound() - pi / 3, 2.5e-15);
    BOOST_CHECK_SMALL(acos_x.rightBound() - acos_1_3, 2.3e-15);
   
    interval y(1.0, 2.0);

    constexpr double e_constant = 2.71828182845904523536;
    auto exp_y = exp(y);
    BOOST_CHECK_SMALL(exp_y.leftBound() - e_constant, 0.0);
    BOOST_CHECK_SMALL(exp_y.rightBound() - e_constant * e_constant, 5.4e-15);

    constexpr double log_e_2 = 0.6931471805599453094;
    auto ln_y = log(y);
    BOOST_CHECK_SMALL(ln_y.leftBound(), 0.0);
    BOOST_CHECK_SMALL(ln_y.rightBound() - log_e_2, 3.4e-16);

    auto pow_y_2_3 = power(y, interval(2,3));
    BOOST_CHECK_SMALL(pow_y_2_3.leftBound() - 1.0, 0.0);
    BOOST_CHECK_SMALL(pow_y_2_3.rightBound() - 8.0, 1.6e-14);
   
    y = interval(-3, -2);

    auto sqr_y = y^2;
    BOOST_CHECK_EQUAL(sqr_y, interval(4, 9));

    auto cube_y = y^3;
    BOOST_CHECK_EQUAL(cube_y, interval(-27, -8));

    auto sqr_e = e^2;
    BOOST_CHECK_EQUAL(sqr_e, interval(0, 16));

    y = interval(4, 9);
    auto pow_y = power(y, interval(-1.5));
    BOOST_CHECK_SMALL(pow_y.leftBound() - 1.0/27, 1.2e-16);
    BOOST_CHECK_SMALL(pow_y.rightBound() - 1.0/8, 2.8e-16);
}

BOOST_AUTO_TEST_CASE(multiplication)
{
    interval a(1.0, 2.0);

    BOOST_CHECK_EQUAL(a * interval(1.0, 2.0), interval(1.0, 4.0));
    BOOST_CHECK_EQUAL(a * interval(-2.0, -1.0), interval(-4.0, -1.0));
    BOOST_CHECK_EQUAL(a * interval(-3.0, 2.0), interval(-6.0, 4.0));

    a = interval(-2.0, -1.0);
    BOOST_CHECK_EQUAL(a * interval(1.0, 2.0), interval(-4.0, -1.0));
    BOOST_CHECK_EQUAL(a * interval(-2.0, -1.0), interval(1.0, 4.0));
    BOOST_CHECK_EQUAL(a * interval(-3.0, 2.0), interval(-4.0, 6.0));

    a = interval(-1.0, 3.0);
    BOOST_CHECK_EQUAL(a * interval(1.0, 2.0), interval(-2.0, 6.0));
    BOOST_CHECK_EQUAL(a * interval(-2.0, -1.0), interval(-6.0, 2.0));
    BOOST_CHECK_EQUAL(a * interval(-3.0, 2.0), interval(-9.0, 6.0));
}

BOOST_AUTO_TEST_CASE(multiplication_with_assignment)
{
    interval a(1.0, 2.0);

    a *= interval(1.0, 2.0);
    BOOST_CHECK_EQUAL(a, interval(1.0, 4.0));

    a = interval(1.0, 2.0);
    a *= interval(-2.0, -1.0);
    BOOST_CHECK_EQUAL(a, interval(-4.0, -1.0));

    a = interval(1.0, 2.0);
    a *= interval(-3.0, 2.0);
    BOOST_CHECK_EQUAL(a, interval(-6.0, 4.0));

    a = interval(-2.0, -1.0);
    a *= interval(1.0, 2.0);
    BOOST_CHECK_EQUAL(a, interval(-4.0, -1.0));

    a = interval(-2.0, -1.0);
    a *= interval(-2.0, -1.0);
    BOOST_CHECK_EQUAL(a, interval(1.0, 4.0));

    a = interval(-2.0, -1.0);
    a *= interval(-3.0, 2.0);
    BOOST_CHECK_EQUAL(a, interval(-4.0, 6.0));

    a = interval(-1.0, 3.0);
    a *= interval(1.0, 2.0);
    BOOST_CHECK_EQUAL(a, interval(-2.0, 6.0));

    a = interval(-1.0, 3.0);
    a *= interval(-2.0, -1.0);
    BOOST_CHECK_EQUAL(a, interval(-6.0, 2.0));

    a = interval(-1.0, 3.0);
    a *= interval(-3.0, 2.0);
    BOOST_CHECK_EQUAL(a, interval(-9.0, 6.0));
}

BOOST_AUTO_TEST_CASE(scalar)
{
    interval a(1.0, 2.0);

    BOOST_CHECK_EQUAL(a * 2, interval(2.0, 4.0));
    BOOST_CHECK_EQUAL(a * 2, interval(2.0, 4.0));
    BOOST_CHECK_EQUAL(2 * a, interval(2.0, 4.0));

    BOOST_CHECK_EQUAL(a / 2, interval(0.5, 1.0));
    BOOST_CHECK_EQUAL(a / 2, interval(0.5, 1.0));
    BOOST_CHECK_EQUAL(2 / a, interval(1.0, 2.0));

    BOOST_CHECK_EQUAL(a + 2, interval(3.0, 4.0));
    BOOST_CHECK_EQUAL(a + 2, interval(3.0, 4.0));
    BOOST_CHECK_EQUAL(2 + a, interval(3.0, 4.0));

    BOOST_CHECK_EQUAL(a - 2, interval(-1.0, 0.0));
    BOOST_CHECK_EQUAL(a - 2, interval(-1.0, 0.0));
    BOOST_CHECK_EQUAL(2 - a, interval(0.0, 1.0));
}

BOOST_AUTO_TEST_CASE(division)
{
     interval a(1.0, 2.0);
     BOOST_CHECK_EQUAL(a / interval(1.0, 2.0), interval(0.5, 2.0));
     BOOST_CHECK_EQUAL(a / interval(-2.0, -1.0), interval(-2.0, -0.5));

     a = interval(-2.0, -1.0);
     BOOST_CHECK_EQUAL(a / interval(1.0, 2.0), interval(-2.0, -0.5));
     BOOST_CHECK_EQUAL(a / interval(-2.0, -1.0), interval(0.5, 2.0));
     
     a = interval(-1.0, 3.0);
     BOOST_CHECK_EQUAL(a / interval(1.0, 2.0), interval(-1.0, 3.0));
     BOOST_CHECK_EQUAL(a / interval(-2.0, -1.0), interval(-3.0, 1.0));
}

BOOST_AUTO_TEST_CASE(divison_with_assignment)
{
     interval a(1.0, 2.0);
     a /= interval(1.0, 2.0);
     BOOST_CHECK_EQUAL(a, interval(0.5, 2.0));

     a = interval(1.0, 2.0);
     a /= interval(-2.0, -1.0);
     BOOST_CHECK_EQUAL(a, interval(-2.0, -0.5));

     a = interval(-2.0, -1.0);
     a /= interval(1.0, 2.0);
     BOOST_CHECK_EQUAL(a, interval(-2.0, -0.5));

     a = interval(-2.0, -1.0);
     a /= interval(-2.0, -1.0);
     BOOST_CHECK_EQUAL(a, interval(0.5, 2.0));
 
     a = interval(-1.0, 3.0);
     a /= interval(1.0, 2.0);
     BOOST_CHECK_EQUAL(a, interval(-1.0, 3.0));

     a = interval(-1.0, 3.0);
     a /= interval(-2.0, -1.0);
     BOOST_CHECK_EQUAL(a, interval(-3.0, 1.0));
}

BOOST_AUTO_TEST_CASE(irrational_division)
{
    auto div_1_10 = interval(1.0) / interval(10.0);

    BOOST_CHECK_LT(div_1_10.leftBound(), div_1_10.rightBound());
    BOOST_CHECK_SMALL(div_1_10.leftBound() - 1.0/10, 1.4e-17);
    BOOST_CHECK_SMALL(div_1_10.rightBound() - 1.0/10, 0.0);
}

BOOST_AUTO_TEST_CASE(postive_negative_zero_comparison)
{
    interval a = interval(+0.0, +0.0);
    interval b = interval(-0.0, +0.0);
    interval c = interval(-0.0, -0.0);

    BOOST_CHECK_EQUAL(a,b);
    BOOST_CHECK_EQUAL(b,c);
    BOOST_CHECK_EQUAL(c,a);

    BOOST_CHECK(a == b);
    BOOST_CHECK(b == c);
    BOOST_CHECK(c == a);
}
