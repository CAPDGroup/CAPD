/// @file intervalCompare.cpp
///
/// Correctness and speed test of the interval arithmetic selected at compile time
/// (FILIB, TIGHT, CXSC or NATIVE - CAPD_INTERVAL_TYPE), used through the capd::interval interface.
///
/// Correctness. For every operation the result is compared with the exact (MPFR, 300 bits) hull of
/// the set of values of the operation on the argument(s). The computed interval must contain this hull
/// (otherwise the arithmetic is not rigorous - a "violation"). The overestimation is measured in ulps.
///
/// Speed. Time (ns per operation) of a loop over precomputed random arguments.
///
/// Each backend is a separate build of CAPD, hence each build has its own executable. They write
/// a CSV file; compare_backends.py puts the results of several builds side by side.
///
/// Quality vs. diameter of the argument (--width): arguments of diameter 0 (points), 1e-13, 1e-10, 1e-8;
/// the diameter of the result is compared with the diameter of the exact hull.
///
/// Usage:  IntervalCompare [--quick|--full] [--width] [--csv file] [--strict]

#include <algorithm>
#include <chrono>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <mpfr.h>

#include "capd/intervals/lib.h"

#if defined(__USE_TIGHT__)
#define BACKEND "tight"
#elif defined(__USE_FILIB__)
#define BACKEND "filib"
#elif defined(__USE_CXSC__)
#define BACKEND "cxsc"
#else
#define BACKEND "native"
#endif

typedef capd::interval Interval;

static const mpfr_prec_t PREC = 300;

// ---------------------------------------------------------------------------------------------
// small RAII wrapper of mpfr_t

struct MP {
  mpfr_t v;
  MP() { mpfr_init2(v, PREC); }
  explicit MP(double d) { mpfr_init2(v, PREC); mpfr_set_d(v, d, MPFR_RNDN); }
  MP(const MP& o) { mpfr_init2(v, PREC); mpfr_set(v, o.v, MPFR_RNDN); }
  MP& operator=(const MP& o) { mpfr_set(v, o.v, MPFR_RNDN); return *this; }
  ~MP() { mpfr_clear(v); }
};

// hull [lo,hi] of the exact set of values (computed with directed rounding, so it is a rigorous outer hull)
struct Hull {
  MP lo, hi;
  Hull() { mpfr_set_inf(lo.v, +1); mpfr_set_inf(hi.v, -1); }
  void add(const MP& d, const MP& u) {
    if (mpfr_less_p(d.v, lo.v)) mpfr_set(lo.v, d.v, MPFR_RNDD);
    if (mpfr_greater_p(u.v, hi.v)) mpfr_set(hi.v, u.v, MPFR_RNDU);
  }
};

typedef int (*MpFun)(mpfr_t, const mpfr_srcptr, mpfr_rnd_t);

// ---------------------------------------------------------------------------------------------
// deterministic random numbers (all doubles are generated exactly, independently of the FPU rounding mode)

struct Rng {
  uint64_t s;
  explicit Rng(uint64_t seed) : s(seed) {}
  uint64_t next() { s ^= s << 13; s ^= s >> 7; s ^= s << 17; return s; }
  int range(int a, int b) { return a + static_cast<int>(next() % static_cast<uint64_t>(b - a + 1)); }
};

// mantissa in [2^52,2^53), value = m * 2^(e-53)
static double makeDouble(uint64_t m, int e) { return std::ldexp(static_cast<double>(m), e - 53); }

// positive interval [lo,hi]; kind: 0 - point, 1 - narrow, 2 - wide
static void positiveInterval(Rng& r, int minExp, int maxExp, double& lo, double& hi) {
  int kind = r.range(0, 3);
  int e = r.range(minExp, maxExp);
  uint64_t m = (1ULL << 52) | (r.next() & ((1ULL << 52) - 1));
  if (kind == 0) { lo = hi = makeDouble(m, e); return; }
  if (kind == 1) {
    uint64_t d = r.next() & ((1ULL << 20) - 1);
    if (m + d >= (1ULL << 53)) m -= (1ULL << 21);
    lo = makeDouble(m, e); hi = makeDouble(m + d, e); return;
  }
  int e2 = r.range(minExp, maxExp);
  uint64_t m2 = (1ULL << 52) | (r.next() & ((1ULL << 52) - 1));
  lo = makeDouble(m, e); hi = makeDouble(m2, e2);
  if (lo > hi) std::swap(lo, hi);
}

// interval around signed numbers, possibly containing zero (wide intervals use two independent endpoints)
static void anyInterval(Rng& r, int minExp, int maxExp, double& lo, double& hi) {
  double a, b;
  positiveInterval(r, minExp, maxExp, a, b);
  int sign = r.range(0, 2); // 0: positive, 1: negative, 2: mixed (only for wide)
  if (sign == 0) { lo = a; hi = b; }
  else if (sign == 1 || a == b) { lo = -b; hi = -a; }
  else { lo = -a; hi = b; }
}

// ---------------------------------------------------------------------------------------------
// operations

enum Arg { ANY, POSITIVE, NONZERO };

struct Op {
  std::string name;
  int arity;
  Arg domain;                 // domain of the (first/second) argument
  int minExp, maxExp;         // range of binary exponents of arguments
  std::function<Interval(const Interval&, const Interval&)> apply;
  // exact hull of the values of the operation on the argument intervals
  std::function<void(Hull&, const MP&, const MP&, const MP&, const MP&)> reference;
};

static void hullMonotone(Hull& h, MpFun f, const MP& a, const MP& b) {
  MP d, u;
  f(d.v, a.v, MPFR_RNDD); f(u.v, b.v, MPFR_RNDU);
  h.add(d, u);
}

// sin (shift=0) and cos (shift=1) hull on [a,b]. Critical points of sin: pi/2 + m*pi, of cos: m*pi
static void hullTrig(Hull& h, bool isSin, const MP& a, const MP& b) {
  MpFun f = isSin ? static_cast<MpFun>(mpfr_sin) : static_cast<MpFun>(mpfr_cos);
  MP da, ua, db, ub;
  f(da.v, a.v, MPFR_RNDD); f(ua.v, a.v, MPFR_RNDU);
  f(db.v, b.v, MPFR_RNDD); f(ub.v, b.v, MPFR_RNDU);
  h.add(da, ua); h.add(db, ub);
  MP pi, shift, t1, t2;
  mpfr_const_pi(pi.v, MPFR_RNDN);
  mpfr_set_zero(shift.v, +1);
  if (isSin) mpfr_div_2ui(shift.v, pi.v, 1, MPFR_RNDN);
  mpfr_sub(t1.v, a.v, shift.v, MPFR_RNDN); mpfr_div(t1.v, t1.v, pi.v, MPFR_RNDN); mpfr_ceil(t1.v, t1.v);
  mpfr_sub(t2.v, b.v, shift.v, MPFR_RNDN); mpfr_div(t2.v, t2.v, pi.v, MPFR_RNDN); mpfr_floor(t2.v, t2.v);
  if (mpfr_greater_p(t1.v, t2.v)) return;           // no critical point inside
  MP one(1.0), mone(-1.0);
  if (mpfr_cmp(t2.v, t1.v) >= 1) { h.add(mone, mone); h.add(one, one); return; }
  long m = mpfr_get_si(t1.v, MPFR_RNDN);
  if (m % 2 == 0) h.add(one, one); else h.add(mone, mone);
}

static std::vector<Op> makeOps() {
  std::vector<Op> ops;
  ops.push_back({"add", 2, ANY, -30, 30,
    [](const Interval& x, const Interval& y) { return x + y; },
    [](Hull& h, const MP& a, const MP& b, const MP& c, const MP& d) {
      MP lo, hi; mpfr_add(lo.v, a.v, c.v, MPFR_RNDD); mpfr_add(hi.v, b.v, d.v, MPFR_RNDU); h.add(lo, hi); }});
  ops.push_back({"sub", 2, ANY, -30, 30,
    [](const Interval& x, const Interval& y) { return x - y; },
    [](Hull& h, const MP& a, const MP& b, const MP& c, const MP& d) {
      MP lo, hi; mpfr_sub(lo.v, a.v, d.v, MPFR_RNDD); mpfr_sub(hi.v, b.v, c.v, MPFR_RNDU); h.add(lo, hi); }});
  ops.push_back({"mul", 2, ANY, -30, 30,
    [](const Interval& x, const Interval& y) { return x * y; },
    [](Hull& h, const MP& a, const MP& b, const MP& c, const MP& d) {
      const MP* p[2] = {&a, &b}; const MP* q[2] = {&c, &d};
      for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) {
        MP lo, hi; mpfr_mul(lo.v, p[i]->v, q[j]->v, MPFR_RNDD); mpfr_mul(hi.v, p[i]->v, q[j]->v, MPFR_RNDU); h.add(lo, hi); } }});
  ops.push_back({"div", 2, NONZERO, -30, 30,
    [](const Interval& x, const Interval& y) { return x / y; },
    [](Hull& h, const MP& a, const MP& b, const MP& c, const MP& d) {
      const MP* p[2] = {&a, &b}; const MP* q[2] = {&c, &d};
      for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) {
        MP lo, hi; mpfr_div(lo.v, p[i]->v, q[j]->v, MPFR_RNDD); mpfr_div(hi.v, p[i]->v, q[j]->v, MPFR_RNDU); h.add(lo, hi); } }});
  ops.push_back({"sqr", 1, ANY, -30, 30,
    [](const Interval& x, const Interval&) { return sqr(x); },
    [](Hull& h, const MP& a, const MP& b, const MP&, const MP&) {
      MP lo, hi, l1, h1, l2, h2;
      mpfr_mul(l1.v, a.v, a.v, MPFR_RNDD); mpfr_mul(h1.v, a.v, a.v, MPFR_RNDU);
      mpfr_mul(l2.v, b.v, b.v, MPFR_RNDD); mpfr_mul(h2.v, b.v, b.v, MPFR_RNDU);
      h.add(l1, h1); h.add(l2, h2);
      if (mpfr_sgn(a.v) <= 0 && mpfr_sgn(b.v) >= 0) { MP z(0.0); h.add(z, z); } }});
  ops.push_back({"sqrt", 1, POSITIVE, -60, 60,
    [](const Interval& x, const Interval&) { return sqrt(x); },
    [](Hull& h, const MP& a, const MP& b, const MP&, const MP&) { hullMonotone(h, mpfr_sqrt, a, b); }});
  ops.push_back({"exp", 1, ANY, -30, 5,
    [](const Interval& x, const Interval&) { return exp(x); },
    [](Hull& h, const MP& a, const MP& b, const MP&, const MP&) { hullMonotone(h, mpfr_exp, a, b); }});
  ops.push_back({"log", 1, POSITIVE, -60, 60,
    [](const Interval& x, const Interval&) { return log(x); },
    [](Hull& h, const MP& a, const MP& b, const MP&, const MP&) { hullMonotone(h, mpfr_log, a, b); }});
  ops.push_back({"atan", 1, ANY, -30, 30,
    [](const Interval& x, const Interval&) { return atan(x); },
    [](Hull& h, const MP& a, const MP& b, const MP&, const MP&) { hullMonotone(h, mpfr_atan, a, b); }});
  ops.push_back({"sin", 1, ANY, -30, 4,
    [](const Interval& x, const Interval&) { return sin(x); },
    [](Hull& h, const MP& a, const MP& b, const MP&, const MP&) { hullTrig(h, true, a, b); }});
  ops.push_back({"cos", 1, ANY, -30, 4,
    [](const Interval& x, const Interval&) { return cos(x); },
    [](Hull& h, const MP& a, const MP& b, const MP&, const MP&) { hullTrig(h, false, a, b); }});
  return ops;
}

static void generate(Rng& r, const Op& op, double& lo, double& hi) {
  switch (op.domain) {
    case POSITIVE: positiveInterval(r, op.minExp, op.maxExp, lo, hi); break;
    case NONZERO: {
      positiveInterval(r, op.minExp, op.maxExp, lo, hi);
      if (r.range(0, 1)) { double t = lo; lo = -hi; hi = -t; }
      break; }
    default: anyInterval(r, op.minExp, op.maxExp, lo, hi);
  }
}

// ---------------------------------------------------------------------------------------------

struct Result {
  std::string name;
  long count = 0, violations = 0, errors = 0;
  double sumExcess = 0, maxExcess = 0, nsPerOp = 0;
};

static double ulpOf(double x) {
  x = std::fabs(x);
  if (x == 0.0) return std::numeric_limits<double>::denorm_min();
  return std::nextafter(x, std::numeric_limits<double>::infinity()) - x;
}

// excess of the computed interval [l,r] over the hull, in ulps; sets violated if the hull is not contained
static double excessUlps(double l, double r, const Hull& h, bool& violated) {
  MP ml(l), mr(r), dl, dr;
  violated = false;
  if (mpfr_greater_p(ml.v, h.lo.v) || mpfr_less_p(mr.v, h.hi.v)) violated = true;
  mpfr_sub(dl.v, h.lo.v, ml.v, MPFR_RNDN);
  mpfr_sub(dr.v, mr.v, h.hi.v, MPFR_RNDN);
  double scale = std::max(std::fabs(l), std::fabs(r));
  double hi = mpfr_get_d(h.hi.v, MPFR_RNDN), lo = mpfr_get_d(h.lo.v, MPFR_RNDN);
  scale = std::min(scale, std::max(std::fabs(lo), std::fabs(hi)) == 0 ? scale : std::max(std::fabs(lo), std::fabs(hi)));
  double u = ulpOf(scale);
  double ex = (std::max(mpfr_get_d(dl.v, MPFR_RNDN), 0.0) + std::max(mpfr_get_d(dr.v, MPFR_RNDN), 0.0)) / u;
  return ex;
}

// ---------------------------------------------------------------------------------------------
// Quality of enclosures as a function of the diameter of the argument (--width mode).
//
// Arguments have a prescribed absolute diameter w (0 means a point interval) and moderate magnitude
// (|x| in [0.25,4)), so that the diameter is meaningful for all operations. For each argument we compare
// the diameter of the computed result with the diameter of the exact hull.

// interval [x,x+w] (or [-x-w,-x] if neg) with x in [2^(eMin-1), 2^eMax); x+w is rounded by MPFR, so the
// arguments do not depend on the FPU rounding mode (i.e. they are identical for all backends)
static void fixedWidthArgument(Rng& r, bool allowNegative, double w, double& lo, double& hi) {
  int e = r.range(-1, 2);
  uint64_t m = (1ULL << 52) | (r.next() & ((1ULL << 52) - 1));
  double x = makeDouble(m, e);
  bool neg = allowNegative && r.range(0, 1);
  double xw = x;
  if (w > 0) { MP t(x); mpfr_add_d(t.v, t.v, w, MPFR_RNDN); xw = mpfr_get_d(t.v, MPFR_RNDN); }
  if (neg) { lo = -xw; hi = -x; } else { lo = x; hi = xw; }
}

// Statistics and text output must not run in the rounding mode left by the interval library (TIGHT keeps
// the FPU in round-up: sums would be biased and printf/iostream would round the printed digits upwards).
struct NearestRounding {
  int old;
  NearestRounding() : old(std::fegetround()) { std::fesetround(FE_TONEAREST); }
  ~NearestRounding() { std::fesetround(old); }
};

static double mpDiff(const MP& hi, const MP& lo) {
  MP d; mpfr_sub(d.v, hi.v, lo.v, MPFR_RNDN); return mpfr_get_d(d.v, MPFR_RNDN);
}

static int runWidthTest(bool full, bool strict, const std::string& csv) {
  const int M = full ? 20000 : 2000;
  const double widths[] = {0.0, 1.0e-13, 1.0e-10, 1.0e-8};
  std::vector<Op> ops = makeOps();
  std::ofstream out(csv);
  out << std::setprecision(12);   // ratio of diameters differs from 1 only in the 8th digit for width 1e-8
  out << "backend,op,width,count,violations,errors,mean_diam,mean_hull_diam,median_ratio,max_ratio,mean_excess_ulp,max_excess_ulp\n";
  std::printf("backend: %s   mode: width, %s   arguments per operation and width: %d\n", BACKEND, full ? "full" : "quick", M);
  std::printf("%-6s %9s %6s %12s %12s %12s %12s %12s\n", "op", "width", "viol", "mean diam", "hull diam", "median ratio", "mean exc/ulp", "max exc/ulp");
  long bad = 0;
  for (const Op& op : ops) {
    for (double w : widths) {
      Rng rng(0xC0FFEE1234567ULL ^ std::hash<std::string>()(op.name));
      long viol = 0, err = 0, n = 0;
      double sumDiam = 0, sumHull = 0, sumExc = 0, maxExc = 0, maxRatio = 0;
      std::vector<double> ratios;
      for (int i = 0; i < M; ++i) {
        double l1, h1, l2 = 1.0, h2 = 1.0;
        fixedWidthArgument(rng, op.domain != POSITIVE, w, l1, h1);
        if (op.arity == 2) fixedWidthArgument(rng, op.domain != POSITIVE, w, l2, h2);
        try {
          Interval r = op.apply(Interval(l1, h1), Interval(l2, h2));
          NearestRounding nearest;
          Hull h;
          op.reference(h, MP(l1), MP(h1), MP(l2), MP(h2));
          bool violated;
          double ex = excessUlps(r.leftBound(), r.rightBound(), h, violated);
          if (violated) ++viol;
          double dr = mpDiff(MP(r.rightBound()), MP(r.leftBound()));
          double dh = mpDiff(h.hi, h.lo);
          sumDiam += dr; sumHull += dh; sumExc += ex; maxExc = std::max(maxExc, ex);
          if (w > 0 && dh > 0) { double q = dr / dh; ratios.push_back(q); maxRatio = std::max(maxRatio, q); }
          ++n;
        } catch (...) { ++err; }
      }
      NearestRounding nearest;
      double median = 0;
      if (!ratios.empty()) {
        std::nth_element(ratios.begin(), ratios.begin() + ratios.size() / 2, ratios.end());
        median = ratios[ratios.size() / 2];
      }
      double cnt = n > 0 ? static_cast<double>(n) : 1.0;
      std::printf("%-6s %9.0e %6ld %12.3e %12.3e %12.6f %12.3f %12.1f\n", op.name.c_str(), w, viol,
                  sumDiam / cnt, sumHull / cnt, median, sumExc / cnt, maxExc);
      out << BACKEND << ',' << op.name << ',' << w << ',' << n << ',' << viol << ',' << err << ',' << sumDiam / cnt << ','
          << sumHull / cnt << ',' << median << ',' << maxRatio << ',' << sumExc / cnt << ',' << maxExc << '\n';
      bad += viol + err;
    }
  }
  std::printf("results written to %s\n", csv.c_str());
  return (strict && bad > 0) ? 1 : 0;
}

int main(int argc, char** argv) {
  bool full = false, strict = false, widthMode = false, csvGiven = false;
  std::string csv = std::string("intervals_") + BACKEND + ".csv";
  for (int i = 1; i < argc; ++i) {
    if (!std::strcmp(argv[i], "--full")) full = true;
    else if (!std::strcmp(argv[i], "--quick")) full = false;
    else if (!std::strcmp(argv[i], "--strict")) strict = true;
    else if (!std::strcmp(argv[i], "--width")) widthMode = true;
    else if (!std::strcmp(argv[i], "--csv") && i + 1 < argc) { csv = argv[++i]; csvGiven = true; }
    else { std::cerr << "usage: " << argv[0] << " [--quick|--full] [--width] [--csv file] [--strict]\n"; return 2; }
  }
  if (widthMode) {
    if (!csvGiven) csv = std::string("intervals_width_") + BACKEND + ".csv";
    return runWidthTest(full, strict, csv);
  }
  const int M = full ? 20000 : 2000;      // arguments per operation
  const int reps = full ? 200 : 20;       // repetitions in the speed test

  std::vector<Op> ops = makeOps();
  std::vector<Result> results;
  long totalViolations = 0, totalErrors = 0;

  for (const Op& op : ops) {
    Rng rng(0x9E3779B97F4A7C15ULL ^ std::hash<std::string>()(op.name));
    std::vector<double> lo1(M), hi1(M), lo2(M), hi2(M);
    for (int i = 0; i < M; ++i) {
      generate(rng, op, lo1[i], hi1[i]);
      if (op.arity == 2) {
        // second argument: for div it must not contain zero, for others any
        generate(rng, op, lo2[i], hi2[i]);
      } else { lo2[i] = hi2[i] = 1.0; }
    }
    std::vector<Interval> A(M), B(M);
    for (int i = 0; i < M; ++i) { A[i] = Interval(lo1[i], hi1[i]); B[i] = Interval(lo2[i], hi2[i]); }

    Result res; res.name = op.name; res.count = M;
    // ---- correctness
    for (int i = 0; i < M; ++i) {
      try {
        Interval r = op.apply(A[i], B[i]);
        Hull h;
        op.reference(h, MP(lo1[i]), MP(hi1[i]), MP(lo2[i]), MP(hi2[i]));
        bool violated;
        double ex = excessUlps(r.leftBound(), r.rightBound(), h, violated);
        if (violated) {
          ++res.violations;
          if (res.violations <= 3)
            std::fprintf(stderr, "[%s] %s violation: arg=[%.17g,%.17g]%s result=[%.17g,%.17g]\n", BACKEND, op.name.c_str(),
                         lo1[i], hi1[i], op.arity == 2 ? " (+second)" : "", r.leftBound(), r.rightBound());
        }
        res.sumExcess += ex; res.maxExcess = std::max(res.maxExcess, ex);
      } catch (...) { ++res.errors; }
    }
    // ---- speed
    volatile double sink = 0;
    try {
      auto t0 = std::chrono::steady_clock::now();
      double acc = 0;
      for (int k = 0; k < reps; ++k)
        for (int i = 0; i < M; ++i) acc += op.apply(A[i], B[i]).rightBound();
      auto t1 = std::chrono::steady_clock::now();
      sink = acc;
      res.nsPerOp = std::chrono::duration<double, std::nano>(t1 - t0).count() / (static_cast<double>(reps) * M);
    } catch (...) { res.nsPerOp = -1; ++res.errors; }
    (void)sink;
    totalViolations += res.violations; totalErrors += res.errors;
    results.push_back(res);
  }

  std::printf("backend: %s   mode: %s   arguments per operation: %d\n", BACKEND, full ? "full" : "quick", M);
  std::printf("%-6s %10s %10s %14s %14s %12s\n", "op", "violations", "errors", "mean excess/ulp", "max excess/ulp", "ns/op");
  for (const Result& r : results)
    std::printf("%-6s %10ld %10ld %14.3f %14.1f %12.1f\n", r.name.c_str(), r.violations, r.errors,
                r.sumExcess / r.count, r.maxExcess, r.nsPerOp);

  std::ofstream out(csv);
  out << "backend,op,count,violations,errors,mean_excess_ulp,max_excess_ulp,ns_per_op\n";
  for (const Result& r : results)
    out << BACKEND << ',' << r.name << ',' << r.count << ',' << r.violations << ',' << r.errors << ','
        << r.sumExcess / r.count << ',' << r.maxExcess << ',' << r.nsPerOp << '\n';
  std::printf("results written to %s\n", csv.c_str());

  return (strict && (totalViolations > 0 || totalErrors > 0)) ? 1 : 0;
}
