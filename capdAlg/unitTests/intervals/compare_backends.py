#!/usr/bin/env python3
"""Puts side by side the results of IntervalCompare from several CAPD builds
(one build directory per interval backend: FILIB, TIGHT, CXSC, NATIVE).

Usage:
    compare_backends.py [--run quick|full] build_filib/.../intervals_FILIB.csv ...
    compare_backends.py --run quick DIR1 DIR2 ...        # DIR = CAPD build directory, the executable is started in it

With --run the IntervalCompare executable is run in every given build directory first.
Without --run the given CSV files (or build directories containing intervals_*.csv) are only read.
Exit status is 1 if any backend returned an interval that does not contain the exact result.
"""
import csv
import glob
import os
import subprocess
import sys

EXE = os.path.join("capdAlg", "unitTests", "intervals", "IntervalCompare")


def find_csvs(path):
    if os.path.isdir(path):
        return glob.glob(os.path.join(path, "**", "intervals_*.csv"), recursive=True)
    return [path]


def main(argv):
    mode = None
    paths = []
    it = iter(argv)
    for a in it:
        if a == "--run":
            mode = next(it)
        else:
            paths.append(a)
    if not paths:
        print(__doc__)
        return 2

    csvs = []
    for p in paths:
        if mode:
            exe = os.path.join(p, EXE)
            out = os.path.join(p, "intervals_run.csv")
            subprocess.run([exe, "--" + mode, "--csv", out], check=True, stdout=subprocess.DEVNULL)
            csvs.append(out)
        else:
            csvs += find_csvs(p)

    data = {}   # backend -> op -> row
    ops = []
    for f in csvs:
        with open(f) as fh:
            for row in csv.DictReader(fh):
                data.setdefault(row["backend"], {})[row["op"]] = row
                if row["op"] not in ops:
                    ops.append(row["op"])
    backends = sorted(data)

    def table(title, key, fmt):
        print("\n" + title)
        print("%-6s" % "op" + "".join("%14s" % b for b in backends))
        for op in ops:
            line = "%-6s" % op
            for b in backends:
                r = data[b].get(op)
                line += "%14s" % (fmt(float(r[key])) if r else "-")
            print(line)

    table("violations (result does not contain the exact hull; 0 = rigorous)", "violations", lambda v: "%d" % v)
    table("mean overestimation [ulp]", "mean_excess_ulp", lambda v: "%.2f" % v)
    table("max overestimation [ulp]", "max_excess_ulp", lambda v: "%.1f" % v)
    table("time [ns/op]", "ns_per_op", lambda v: "%.1f" % v)

    if "filib" in data:
        print("\nspeed relative to filib (<1 - faster than filib)")
        print("%-6s" % "op" + "".join("%14s" % b for b in backends))
        for op in ops:
            base = float(data["filib"][op]["ns_per_op"]) if op in data["filib"] else 0
            line = "%-6s" % op
            for b in backends:
                r = data[b].get(op)
                line += "%14s" % ("%.2f" % (float(r["ns_per_op"]) / base) if r and base > 0 else "-")
            print(line)

    bad = sum(int(r["violations"]) + int(r["errors"]) for b in data.values() for r in b.values())
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
