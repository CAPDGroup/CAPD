
#ifndef CAPD_INTERVALS_INTRA_INTERVAL_ATAN2_H
#define CAPD_INTERVALS_INTRA_INTERVAL_ATAN2_H

#include <iostream>
#include "capd/intervals/IntervalError.h"

namespace capd{
  namespace intervals{
    namespace intra{




// ATAN2 function for point intervals with 4 cases (can be not accurate)
// phi in [-pi, pi+eps]
template <typename interval>
      interval atan2forPointsSimple(interval x, interval y) {
        if (x > 0) { // (-90, 90)
          return atan(y / x);
        }
        if (y > 0) { // [0, 180]
          return -atan(x / y) + 0.5 * interval::pi();
        }
        if (y < 0) {    // (-180, 0)
          return -atan(x / y) - 0.5 * interval::pi();
        }
        if (x < 0) {    // [90, 270]
          return atan(y / x) + interval::pi();
        }
        return interval(-1, 1) * interval::pi();
      }

// ATAN2 function for point intervals (overlapping cases for better approximation)
// phi in [-pi,pi+eps]
template <typename interval>
      interval atan2forPoints(interval x, interval y) {
        if ((x > 0) && (abs(y) <= 2 * x)) {  // [-60, 60]
          return atan(y / x);
        }
        if ((y > 0) && (2 * y > abs(x))) {   // [60, 150]
          return -atan(x / y) + 0.5 * interval::pi();
        }
        if (y < 0) {
          if (-2 * y > abs(x)) {              // [-150, 60]
            return -atan(x / y) - 0.5 * interval::pi();
          }
          if (x < 0) {             // [-180, -120]
            return atan(y / x) - interval::pi();
          }
        } else { // y can still be partially negative              //  [90, 180+eps]
          if (x < 0) {
            return atan(y / x) + interval::pi();
          }
        }
        // in the case we did not fall into one of cases (
        return atan2forPointsSimple(x, y);
      }

// Returns angle in the interval (-\pi, \pi + \eps]
// In case box (x,y) constains semiline {x=0, x<0}
// instead of two intervals (-\pi,b] \cup [a,\pi]
// we return [a, b+ 2 \pi]
template <typename interval>
      interval atan2(interval x, interval y) {
        interval leftAngle{}, rightAngle{};
        if (y > 0) {
          if (x > 0) {
            leftAngle = atan2forPoints(x.right(), y.left());
            rightAngle = atan2forPoints(x.left(), y.right());
          } else if (x < 0) {
            leftAngle = atan2forPoints(x.right(), y.right());
            rightAngle = atan2forPoints(x.left(), y.left());
          } else {// x.contains(0)
            leftAngle = atan2forPoints(x.right(), y.left());
            rightAngle = atan2forPoints(x.left(), y.left());
          }
        } else if (y < 0) {
          if (x > 0) {
            leftAngle = atan2forPoints(x.left(), y.left());
            rightAngle = atan2forPoints(x.right(), y.right());
          } else if (x < 0) {
            leftAngle = atan2forPoints(x.left(), y.right());
            rightAngle = atan2forPoints(x.right(), y.left());
          } else {// x.contains(0)
            leftAngle = atan2forPoints(x.left(), y.right());
            rightAngle = atan2forPoints(x.right(), y.right());
          }
        } else {// y.contains(0)
          if (x > 0) {  // angle in [-pi/4, 2 pi]
            leftAngle = atan2forPoints(x.left(), y.left());
            rightAngle = atan2forPoints(x.left(), y.right());
          } else if (x < 0) {  //  (x,y) contains part of {y=0, x<0}
            leftAngle = atan2forPoints(x.right(), y.right());

            rightAngle = atan2forPoints(x.right(), y.left());
            if (rightAngle < 0) {
              rightAngle += 2.0 * interval::pi();
            }
          } else {// x.contains(0)
            throw capd::intervals::IntervalError<typename interval::BoundType>("atan2: Box (x,y) contains zero.", x.leftBound(),
                                                                      x.rightBound());
          }
          return intervalHull(leftAngle, rightAngle);
        }

        auto result = intervalHull(leftAngle, rightAngle);

        return result;
      }




    }
  }
}
#endif // CAPD_INTERVALS_INTRA_INTERVAL_ATAN2_H
