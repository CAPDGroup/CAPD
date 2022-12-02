//////////////////////////////////////////////////////////////////////////////
///
///  @file TexWriterExample.cpp
///  
///  Example how to use TexWriter class.
///
///  Compile and run:
///    - make
///    - ./TexWriterExample2  < example.txt
///    - pdflatex out.tex
///
///  @author kapela  @date   Mar 30, 2011
//////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2011 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include "capd/basicalg/TexWriter.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include "capd/capdAlglib.h"
using namespace std;
using namespace capd;

int main(int, char**) {

  // TexWriter that will write to screen
  cout.precision(16);
  cout << fixed;
  TexWriter out(cout);
  
  // Second TexWriter will write to file 'out.tex' that can be next compiled with LaTeX.
  char filename[] = "out.tex";
  std::ofstream file(filename);
  TexWriter tex(file);
  tex.setEquationSymbol(TexWriter::InlineEquation);

  tex.writeDocumentHeader("\\textheight=24cm \\textwidth=16cm \\oddsidemargin=0.0cm\\evensidemargin=0.0cm \\def\\ii{\\mbox{i\\,}}");
  tex << "\\textbf{Examples of the intervals styles}\\\\\n";
  tex <<
      "\n\\begin{tabular}{|l|l|l|l|}\\hline\n"
      "   FloatSci & FloatSci & FloatFix & floating point style \\\\\n"
      "      2     &    10    &     5    & precision  \\\\\n"
      "      C     &   Math   &     -    & base style \\\\\n"
      "      +     &          &          & plus style \\\\\\hline \n";

  const int numberOfIntervals = 4;
  DInterval intervals[numberOfIntervals] = { 
    DInterval(1.000000001,1.000000042),
    DInterval(-1.23456789123,-1.23456778912),
    DInterval(1.23456789123e5,1.23456778912e5),
    DInterval(1.2345678912345678e12,1.234567789123456e12)
  };
  for (int k = 0; k < numberOfIntervals; ++k) {

    DInterval x = intervals[k];
    
    cout << x << "\n";
    out << x <<"\n----\n";

    tex << " & & & \\\\\n";
    tex.setFloatStyle(TexWriter::FloatSci).precision(2).setBaseStyle(0).setPlusSymbol(0);
    tex << "$" << x << "$ &";                            // equivalent to tex.write(x) <<" &";
    tex.precision(10).setBaseStyle(1).setPlusSymbol(1);
    tex.write(x) << " & ";
    tex.setFloatStyle(TexWriter::FloatFix).precision(5).write(x);
    tex<< " & "<< printToString(x, 16) << " \\\\\\hline\n";
  }

  tex << "\n\\end{tabular}\\\\\n\n\n";
  
  std::complex<DInterval> c(DInterval(14354.32435,14354.3343543543), DInterval(34253.11111,34253.12312));

  tex.setEquationSymbol(2) << "\\textbf{Complex numbers style} ";
  tex.write(c)<<"\n\n";
  out << c <<"\n===\n";

  tex << "\n\n \\textbf{Examples of the vector styles}";
  double dd[] = { 1.131312, 3.0, 4.53254543345};
  IVector v(3);
  std::copy(dd, dd+3, v.begin());
  v += DInterval(1e-13, 1e-10);
  out << v<<"\n===\n";
  tex.write(v);
  tex.setVectorStyle(1) << "$$" << v << "$$\n";
  tex.setVectorStyle(2);
  tex.write(v);
  tex.writeDocumentFooter();

  std::cout << "\nTeX output was writen to file \"" << filename << "\"."
      << "\nYou can compile it with\n  pdflatex " << filename << "\n";
  return 0;
}
