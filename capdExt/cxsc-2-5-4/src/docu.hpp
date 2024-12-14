/*
**  CXSC is a C++ library for eXtended Scientific Computing (V 2.5.4)
**
**  Copyright (C) 1990-2000 Institut fuer Angewandte Mathematik,
**                          Universitaet Karlsruhe, Germany
**            (C) 2000-2014 Wiss. Rechnen/Softwaretechnologie
**                          Universitaet Wuppertal, Germany   
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Library General Public
**  License as published by the Free Software Foundation; either
**  version 2 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Library General Public License for more details.
**
**  You should have received a copy of the GNU Library General Public
**  License along with this library; if not, write to the Free
**  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/* CVS $Id: docu.hpp,v 1.25 2014/01/30 17:23:44 cxsc Exp $ */

/*!
\mainpage

- \ref mainpage_sec_overview
- \ref mainpage_sec_examples
- \ref mainpage_sec_toolbox
- \ref mainpage_sec_thread
- \ref mainpage_sec_credits

\section mainpage_sec_overview Overview

The speed of digital Computers is ever increasing. While the emphasis in computing was traditionally on speed, more emphasis can now be put on accuracy and reliability of results. Numerical mathematics has devised algorithms which deliver highly accurate and automatically verified results by applying mathematical fixed-point theorems. This means that these computation carry their own accuracy control. However, their implementation requires suitable arithmetic support and powerful programming tools which are not generally available.

Different hardware solutions are available for Personal Computers, Workstations, Mainframes and Super Computers. In particular a vector arithmetic coprocessor for the PC has been developed in VLSI-technology. Language support is available on the basis of FORTRAN, PASCAL, and C (ACRITH-XSC, Fortran-XSC, PASCAL-XSC, and C- XSC). Problem-solving routines with automatic result verification have been developed for many standard problems of numerical analysis as for linear or nonlinear systems of equations, for differential and integral equations, etc. as well as for a large number of applications in the engineering and natural sciences.

Language eXtensions for Scientific Computation provide all features indispensable for modern numerical software development, such as
- Operator concept (user-defined operators)
- Overloading concept
- Module concept
- Dynamic arrays
- Controlled rounding
- Predefined arithmetic data types real, (extended real), complex, interval, complex interval, and corresponding vector and matrix types
- Predefined arithmetic operators of highest accuracy for the arithmetic data types
- Predefined elementary functions of highest accuracy for the arithmetic data types
- Data type dotprecision for the exact representation of dot products
- Library of mathematical problem-solving routines with automatic result verification and high accuracy

\section mainpage_sec_examples Simple programming examples

To learn how to use the C-XSC class library you can read the page \ref cxscexamples and follow the examples described there.

For more experienced users the \ref toolbox is the more interesting place. There you find much more sophisticated examples on how to use the C-XSC class library.

\section mainpage_sec_toolbox C++ Toolbox for Verified Computing

The \ref toolbox is the C++ edition of the <em>Numerical Toolbox for Verified Computing</em>. The programs of the original edition
were written in PASCAL-XSC, a PASCAL eXtension for Scientific Computation.

The methods presented here are practical, reliable, and elegant. They are provided in theory, algorithmic descriptions, and
implementations to solve a number of basic numerical problems in a reliable way.

\section mainpage_sec_prec Selecting the precision of dot product computations

By default, C-XSC uses a software emulation of a long fix-point accumulator for dot product computations, which delivers results with maximum accuracy. This means that all operations by default are only one rounding away from the exact result. However, due to the software emulation these operations are slow. 

The user can change the default behaviour for the overloaded operators by setting the global integer variable opdotprec to a value K. All dot products computed with overloaded operators (for example in a matrix-matrix product), are then carried out in (simulated) K-fold double precision. For K=1, pure floating point operations are used. For K>=2, the DotK algorithm, based on so-called error free transformations is used. Accuracy of the results and computing costs increase with K. For K=0, the default setting, the long accumulator is used as before.

When computing dot product expressions using the dotprecision classes, the precision used in each call of the accumulate-function can be set via the member function set_dotprec of the dotprecision classes. The meaning of the setting is the same as for the operators.

For more information, see:
Zimmer, M.; Krämer, W.; Bohlender, G.; Hofschuster, W.:
Extension of the C-XSC Library with Scalar Products with Selectable Accuracy
Published in: Serdica Journal of Computing, Vol. 4, No. 3, p. 349-370, 2010

\section mainpage_sec_blas Using a BLAS-library

For best performance, it is suggested to make use of the BLAS support of C-XSC. BLAS support can be used for the overloaded operator* of the dense vector and matrix types. To activate it, the global variable opdotprec must be set to 1 and the flag -DCXSC_USE_BLAS must be set during compilation. The program must then be linked to an optimized BLAS and CBLAS library. When activating BLAS support, and intermediate midpoint-radius-format is used for intervals, which might result in some additional overestimation.

Please note that the BLAS library used must allow the changing of the rounding mode. Suitable BLAS-librarier in our tests were ATLAS and the Intel MKL. Not supported at the moment is GoTo BLAS.

For more information, see:

Zimmer, M.; Krämer, W.; Hofschuster, W.:
Using C-XSC in High Performance Computing
Preprint 2009/5, Universität Wuppertal, 2009
Revised version submitted for publication (PARA2010):
Walter Krämer, Michael Zimmer, Werner Hofschuster:
Using C-XSC for High Performance Verified Computing
PARA 2010, Reykjavik, Iceland, Part II, LNCS 7134, Springer-Verlag, pp. 168-178, 2012

\section mainpage_sec_openmp Activating OpenMP support

Some of the C-XSC operators can be parallelized using OpenMP. This is especially helpful when using higher precision computations or if no BLAS-library is available. To activate OpenMP support, set the appropriate OpenMP flag of your compiler and additionally the flag -DCXSC_USE_OPENMP during the compilation of your program.


\section mainpage_sec_thread Using C-XSC in a Multi-Threaded Environment

Since C-XSC is now nearly twenty years old, thread-safety has not been a focus in its
development for a long time. In recent years, a lot of work has been invested
into making C-XSC fit for high performance computing and the threadsafety
of C-XSC has been vastly improved in the process.

More information can be found in:

Zimmer, M.:
"Using C-XSC in a Multi-Threaded Environment"
Preprint BUW-WRSWT 2011/2, Universität Wuppertal, 2011, http://www2.math.uni-wuppertal.de/wrswt/preprints/prep_11_2.pdf


\section mainpage_sec_credits Credits

The work on C-XSC started in 1990 at the Institute for Applied Mathematics (Prof. Kulisch), 
University of Karlsruhe. Many colleagues and scientists have directly and
indirectly contributed to the realization of C-XSC. The authors would like to thank 
each of them for his or her cooperation. Special thanks go to 
U. Allend&ouml;rfer, C. Baumhof, H. Berlejung, H. Bleher, H. B&ouml;hm,
B. Bohl, G. Bohlender, F. Blomquist, K. Braune, H.H. Chen, D. Cordes, 
A. Davidenkoff, H.-C. Fischer, M. Grimmer, K. Gr&uuml;ner, R. Hammer, 
M. Hinz, M. Hocks, B. H&ouml;ffgen, W. Hofschuster, P. Januscke, E. Kaucher, 
R. Kelch, R. Kirchner, R. Klatte, W. Klein,
W. Kr&auml;mer, U. Kulisch, C. Lawo, M. Metzger, W.L. Miranker, M. Neaga, 
M. Neher, D. Ratz, M. Rauch, S. Ritterbusch, S.M. Rump, R. Saier,
D. Schiriaev, L. Schmidt, G. Schumacher, U. Storck, J. Suckf&uuml;ll, 
F. Toussaint, C. Ullrich, W. Walter, S. Wedner, G. Werheit, A. Wiethoff, 
H.W. Wippermann, J. Wolff von Gudenberg and M. Zimmer.

C-XSC is an outcome of an ongoing collaboration of the Institute for Applied Mathematics (Prof. Kulisch), 
University of Karlsruhe and the Institute for Scientific 
Computing/Software Engineering (Prof. Kr&auml;mer), University of Wuppertal. 
For the latest news and up to date software contact http://www.math.uni-wuppertal.de/~xsc/ .

Thanks to the referees for valuable comments and suggestions.

\date February 2014

*/

/*!
\page toolbox C++ Toolbox for Verified Computing

- \ref toolbox_sec_resultverification
- \ref toolbox_sec_history
- \ref toolbox_sec_existingSolutions

\section toolbox_sec_resultverification Why Numerical Result Verification?

Floating-point arithmetic is the fast way to perform scientific and engineering calculations. Today, individual floating-point operations on most computers are of
maximum accuracy in the sense that the rounded result differs at most by 1 unit in the last place from the exact result. However, after two or more operations, the
result can be completely wrong. Computers now carry out up to \f$ 10^{11} \f$ floating-point operations per second. Thus, particular attention must be paid to the reliability of
the computed results. In recent years, techniques have been developed in numerical analysis which make it possible for the computer itself to verify the correctness of
computed results for numerous problems and applications. Moreover, the computer frequently establishes the existence and uniqueness of the solution in this way. For
example, a verified solution of a system of ordinary differential equations is just as valid as a solution obtained by a computer algebra system, which still requires a
valid formula evaluation. Furthermore, the numerical routine remains applicable even if the problem does not have a closed-form solution.

\section toolbox_sec_history A Brief History of Computing

Methods of calclulation have put their stamp on mathematics throughout history. Progress in connection with numbers has meant progress in mathematics also. The
invention of mechanical calculating machines by B. Pascal, G. W. Leibniz, and others, of logarithms and the slide rule, and of many other mechanical and analog
computational aids, represent significant milestones in scientific computing.

A great leap in development took place about 50 years ago with the development of the first electronic computers. This technology iommediately made it possible to
perform arithmetic operations faster than its mechanical or electromechanical predeccors by a factor of about 1000. The great technological gains of this century
would have been impossible without modern computation. Today's automobile, airplane, space travel, modern radio and television, and not least the rapid further
development of the computer itself wre enabled by enormous computational capacity. On the other hand, advances in computer hardware gave massive stimulation
to further development of algorithmic and numerical mathematics.

Further development of circuit technology, which again would not have been possible without the computer, has increased the computational capacity of a processor
by a factor of about \f$ 10^9 \f$, compared to the first electronic computers of the 1950`s. Comparison of the numbers \f$ 10^3 \f$ and \f$ 10^9 \f$ shows that the real computer revolution
took place <em>after</em> the development of the first electronic computers. Remarkably, there was nothing equivalent to this development on the arithmetical-mathematical
side. The enormous advances in computer technology really suggested an attempt to make the computer also more powerful arithmetically. On this, however, mathematicians
exerted hardly any influence. In the area of scientific and engineering applications of cencern here, the computer is used today in especially the same way
as it was in the middle 1950`s. The same four floating-point operations of addition, subtraction, multiplication, and division are still being performed, albeit much faster.

As more and more computations can be done faster and faster, larger scientific problems can be addressed, and it becomes increasingly critical to make scientific
computation more trustworthy than the exclusive use of ordinary floating-point arithmetic. This development began about 25 years ago and has gained a considerable
following.
However, methods for validated computation have yet to achieve universal acceptance.

\section toolbox_sec_existingSolutions Existing solutions for standard Algorithms

- \ref evalofpoly
- Automatic Differentiation
- Nonlinear Equations in One Variable
- Global Optimization (one-dimensional)
- Evaluation of Arithmetic Expressions
- Zeros of Complex Polynomials
- (Fast) Linear Systems of Equations
- Linear Optimization
- Automatic Differentiation for Gradients, Jacobians and Hessians
- Nonlinear Systems of Equations
- Global Optimization (multi-dimensional)
*/

/*!
\page evalofpoly Evaluation of Polynomials

- \ref evalofpoly_sec_overview
- \ref evalofpoly_sec_algodescr
- \ref evalofpoly_sec_implementation

\section evalofpoly_sec_overview Overview

On this page, we consider the evaluation of a polynomial function of a single variable. We usually compute the value of an arithmetic function by replacing
each arithmetic operation by its corresponding floating-point machine operation. Roundoff errors and cancellations sometimes cause the calculated
result to be drastically wrong. For similar reasons, a naive interval evaluation of a polynomial may lead to intervals so large as to be practically useless. Roundoff and
cancellation errors are especially dangerous if we are evaluating a function close to a root, as we will see when we compute verified enclosures of zeroes of polynomials.

We present an algorithm to evaluate a real polynomial \f$ p : R \to R \f$ defined as
\f[
p(t) = \sum \limits_{i=0}^n p_i t^i \;,\; p_i , t \in R \;,\; i = 0 , ... , n \;,\; p_n \not= 0
\f]

We assume the coefficients \f$ p_i \f$ to be representable in the floating-point number system of the host computer. The algorithm achieves maximum accuracy, even in the
neighborhood of a root where cancellation dooms an ordinary floating-point evaluation.

\section evalofpoly_sec_algodescr Algorithmic Description

We present the algorithm <em>RPolyEval</em> for the evaluation of a real polynomial \f$ p(t) = \sum \limits_{i=0}^n p_i t^i \f$ with maximum accuracy. Except for the special cases \f$ n = 0 \f$ and \f$ n = 1 \f$,
which can be calculated directly, an iterative solution method is used. We first compute a floating-point approximation of \f$ p(t) \f$. We then carry out a residual iteration by solving a linear system of equations.
The new solution interval determined in the next step is checked for being of maximum accuracy, i.e. for being exact to one unit in the last place of the mantissa (1 ulp).

\image html "algorpoly.png" "The RPolyEval Algorithm"

\section evalofpoly_sec_implementation Implementation and Examples

We list the C++ program code for the evaluation of a real polynomial with maximum accuracy. Interval variables are named with double characters, e.g. \f$ rr[i] \f$ denotes the interval \f$ [r]_i \f$.

\subsection evalofpoly_subsec_rpoly Module rpoly

The header file of module <em>rpoly</em> supplies the definition of the class <em>RPolynomial</em> representing a real polynomial \f$ p(t) = \sum \limits_{i=0}^n p_i t^i \;,\; p_i , t \in R \f$. Since arithmetic operations
for real polynomials are not required by the algorithm, they are not provided by module <em>rpoly</em>.

\includelineno Modules/rpoly.hpp

The input operator >> echos the indices on the standard output device to give an information on the index of the coefficient actually to be entered.

\includelineno Modules/rpoly.cpp

\subsection evalofpoly_subsec_rpeval Module rpeval

The header file of module <em>rpeval</em> supplies the interface of the function <em>RPolyEval()</em>, for evaluating a real polynomial with maximum accuracy and the interface of the
function <em>RPolyEvalErrMsg()</em> which can be used to get an error message for the error code returned by <em>RPolyEval()</em>.

\includelineno Modules/rpeval.hpp

The function <em>RPolyEval()</em> is an implementation of the algorithm. It computes an approximation and an enclosure of the value of the real polynomial \f$ p(t) = \sum \limits_{i=0}^n p_i t^i \f$ with \f$ p_i \in R \f$ at a point \f$ t \in R \f$.

\includelineno Modules/rpeval.cpp

\subsection evalofpoly_subsec_examples Examples

We consider the polynomial \f$ p(t) = t^4 - 8t^3 + 24t^2 - 32t + 16 \f$ to be evaluated in the neighborhood of the real value \f$ t=2.0001 \f$ and the polynomial \f$ q(t) = -t^3 + 3t^2 - 3t + 1 \f$
to be evaluated in the neighborhood of the real value \f$ t = 1.000005 \f$. To make sure that the arguments are representable on the computer, we use the machine numbers
\f$ t = (2.0001) \f$ and \f$ t = (1.000005) \f$, respectively.

\image html "tbpolyt.png" "Polynomial p(t) in interval [1.9995 , 2.0005]"

\image html "tbpolyq.png" "Polynomial q(t) in interval [0.99999 , 1.00001]"

To illustrate the difficulties that may occur with the calculation of a polynomial value when using floating-point arithmetic, these two polynomials have been evaluated
in floating-point arithmetic for 100 values in the neighborhood of their roots. The corresponding plots are given in the above two figures. These two examples are
solved reliably using the following sample program to call <em>RPolyEval()</em>.

\includelineno Programs/rpe_ex.cpp

Our implementation of the algorithm produces the output listed below. In both cases the floating-point approximation, naively using Horner's nested multiplication
form, does not <em>lie</em> within the verified enclosure of \f$ p(t) \f$ and \f$ q(t)\f$. An even worse side effect of rounding errors occuring during the evaluation using Horner's
scheme is that the standard algorithm returns not only wrong digits but also an incorrect sign for the values of \f$ p( 2.0001) \f$ and \f$ q( 1.000005) \f$. To avoid an
incorrect interpretation of the resulting interval, we stress that this interval is numerically proved to be an enclosure of the exact value of the given polynomials for
the machine numbers \f$ t = (2.0001) \f$ and \f$ t = (1.000005) \f$.
*/

/*!
\page cxscexamples Simple programming examples

On this page we will try to show you, how to use the C-XSC class library with some short and simple examples.

- \ref cxscexamples_sec_ex1
- \ref cxscexamples_sec_ex2
- \ref cxscexamples_sec_ex3
- \ref cxscexamples_sec_ex4
- \ref cxscexamples_sec_ex5
- \ref cxscexamples_sec_ex6
- \ref cxscexamples_sec_ex7
- \ref cxscexamples_sec_ex8

\section cxscexamples_sec_ex1 Example 1 - An Introduction

We will start with the following simple program showing you how to use intervals, compute one operation and print out the result:

\includelineno example.cpp

Let's start investigating the interesting lines. With the line

\skipline interval.hpp

you include the functionality of the interval class of C-XSC in your program. After that you have to inform the compiler about the namespace cxsc - where all classes and methods of C-XSC are stored in - to use C-XSC
without fully qualified identifiers.

\skipline using namespace cxsc

Then you declare your variables and assign adequate values in the following lines.

\skipline interval a
\until b =

Finally you want to print out the result for your desired computation.

\skipline cout

So it's just that easy to use C-XSC.

\section cxscexamples_sec_ex2 Example 2 - Input / Output

\includelineno io.cpp

\section cxscexamples_sec_ex3 Example 3 - Compute all zeros of a function

\includelineno allzeros.cpp

\section cxscexamples_sec_ex4 Example 4 - Interval Newton method

\includelineno inewton.cpp

\section cxscexamples_sec_ex5 Example 5 - 

\includelineno lexample.cpp

\section cxscexamples_sec_ex6 Example 6 - 

\includelineno linewton.cpp

\section cxscexamples_sec_ex7 Example 7 - 

\includelineno rungekutta.cpp

\section cxscexamples_sec_ex8 Example 8 - 

\includelineno trace.cpp


*/
