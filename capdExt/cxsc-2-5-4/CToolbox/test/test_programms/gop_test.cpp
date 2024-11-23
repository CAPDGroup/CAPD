//============================================================================
//
//                              Program/Module
//                                   from
//                 C++ TOOLBOX FOR VERIFIED COMPUTING I
//                         Basic Numerical Problems
//
//      Copyright (c) 1995   Rolf Hammer, Matthias Hocks, Dietmar Ratz
//
// For details on theory, algorithms, and programs, see the book
//
//  R. Hammer, M. Hocks, U. Kulisch, D. Ratz:  C++ Toolbox for
//  Verified Computing I - Basic Numerical Problems. Springer-Verlag,
//  Heidelberg, New York, 1995.
//
//============================================================================
//----------------------------------------------------------------------------
// Example: Global optimization
// This program uses module 'gop' to compute the global optimizers of the
// function of Branin
//
//    fB(x) = sqr(5/pi*x[1] - 51/(40*sqr(pi))*sqr(x[1]) + x[2] - 6)
//
//            + 10*(1-1/8/pi)*cos(x[1]) + 10
//
// and the function of Levy
//             5                         5
//    fL(x) = sum i cos((i-1)x[1] + i)  sum j cos((j+1)x[2] + j)
//            i=1                       j=1
//
//            + sqr(x[1] + 1.42513) + sqr(x[2] + 0.80032)
//
// A starting interval and a tolerance must be entered.
//----------------------------------------------------------------------------
#include <gop.hpp>         // Global optimization
#include <stacksz.hpp>     // To increase stack size for some
                           // special C++ compilers
#include "comparator.hpp"

using namespace cxsc;
using namespace std;

stringstream output;

HessType fBranin ( const HTvector& x )
{
  HessType   hh(2);
  interval   pi = Pi();         // enclosure of pi

  hh = sqr(5.0/pi*x[1] - 51.0/(40.0*sqr(pi))*sqr(x[1]) + x[2] - 6.0)
       + 10.0*(1.0-1.0/8.0/pi)*cos(x[1]) + 10.0;
  return hh;
}

HessType fLevy ( const HTvector& x )
{
  HessType   isum(2), jsum(2), hh(2);
  int        i;

  isum = 0.0; jsum = 0.0;
  for (i = 1; i <= 5; i++) {
    isum = isum + double(i)*cos(double(i-1)*x[1] + double(i));
    jsum = jsum + double(i)*cos(double(i+1)*x[2] + double(i));
  }
  hh = isum * jsum +
       sqr(x[1] + interval(142513.0)/100000.0) +    // Avoid real con-
       sqr(x[2] + interval(80032.0)/100000.0);      // version error
  return hh;
}

//----------------------------------------------------------------------------
// Function for printing and reading informations to call the function
// 'AllGOp()'. This function must be called with the function 'f', a string
// 'Name' containing textual description of that function, and an integer
// 'dim' specifying the dimension of the problem.
//----------------------------------------------------------------------------
void compute ( HTscalar_FctPtr f, char* Name, int dim, istream &input )
{
  ivector   SearchInterval(dim);
  interval  Minimum;
  real      Tolerance;
  imatrix   Opti;
  intvector Unique;
  int       NumberOfOptis, i, Error;

  output << SetPrecision(23,15) << Scientific;   // Output format

  output << "Computing all global minimizers of the " << Name << endl;
  output << "Search interval      : "; input >> SearchInterval;
  output << "Tolerance (relative) : "; input >> Tolerance;
  output << endl;

  AllGOp(f,SearchInterval,Tolerance,
           Opti,Unique,NumberOfOptis,Minimum,Error);

  for (i = 1; i <= NumberOfOptis; i++) {
    output << Opti[i];
    if (Unique[i])
      output << "encloses a locally unique candidate for a global minimizer!";
    else
      output << "may contain a local or global minimizer!";
    output << endl << endl;
  }


  if (NumberOfOptis != 0)
    output << Minimum << endl
         << "encloses the global minimum value!" << endl << endl;
  output << NumberOfOptis << " interval enclosure(s)" << endl;

  if (Error)
    output << endl << AllGOpErrMsg(Error) << endl;
  else if (NumberOfOptis == 1 && Unique[1])
    output << endl << "We have validated that there is "
                    "a unique global optimizer!" << endl;
}

// ARGV-VERSION
int main(int argc, char *argv[]) { 
  int result = 1;
  if (argc > 3) {
    ifstream input(argv[1]); 
    ifstream compinput(argv[2]); 
    ifstream expected(argv[3]);
    
    compute(fBranin,(char*)"Function of Branin",2,input);
    output << endl << endl;
    compute(fLevy,(char*)"Function of Levy",2,input);
  
    ComparatorCTests cmp;
    cmp.debugmode(false);
    cmp.compare(output,compinput);   
    result = cmp.getResult(expected);
    cout << "return code = " << result << endl; 
  }
  return result;
}

