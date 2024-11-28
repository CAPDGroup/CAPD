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

/* CVS $Id: expreval.cpp,v 1.15 2014/01/30 17:49:26 cxsc Exp $ */

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
// File: expreval (implementation)
// Purpose: Computation of an enclosure of the value of a real arithmetic
//    expression composed of the operations +, -, *, /, and ^, where ^
//    denotes exponentiation by an integer.
// Method: Iterative refinement using a defect correction mechanism. The
//    successively computed corrections for the result are stored in a
//    staggered correction format.
// Class Staggered (staggered data representation):
//    Staggered()         : constructors
//    operators =         : assignment of a real argument
//    operators +, -, *, /: both operands of type 'Staggered' or one of
//                          type 'Staggered' and one of type 'real'
//    Power()             : argument of type 'Staggered', exponent of type
//                          'integer' (Exponentiation by an integer)
//    Eval()              : main function for evaluation of an expression
//    EvalErrMsg()        : to get an error message text
// Class StaggArray (arrays of type 'Staggered'):
//    StaggArray()        : constructors
//    ~StaggArray()       : destructor
//    operator []         : component access
//---------------------------------------------------------------------------
#include <stdlib.h>         // For function 'exit()'
#include <string.h>         // String handling
#include <idot.hpp>         // Interval dotprecision type
#include <i_util.hpp>       // Interval utility functions
#include <expreval.hpp>

using namespace cxsc;
using namespace std;

const int
  MaxStaggPrec = 10;      // Maximum number of staggered corrections which can
                          // be stored in a variable of type 'Staggered'.
                          // Recompile this module if its value is changed.
                          //--------------------------------------------------

typedef struct IR         // Type for an intermediate result
{                         //--------------------------------
  Staggered  Entry;       // Staggered entry
  IR        *Next;        // Next entry
} IntRes, *IntResPtr;

//----------------------------------------------------------------------------
// Static variables globally used within other functions this module
//----------------------------------------------------------------------------
static int
  DivByZero,      // Error flag
  InitFlag,       // Signals the initialization process
  ActStaggPrec;   // Actual length of the staggered format

static IntResPtr
  HeadPtr,        // Head of the list of intermediate results
  ActPtr,         // Pointer to the actual intermediate result
  FreePtr;        // Pointer to freed intermediate results

//----------------------------------------------------------------------------
// Error codes used in this module.
//----------------------------------------------------------------------------
const int
  NoError   = 0,   // No error occurred
  ItFailed  = 1,   // Maximal number of defect corrections exceeded
  ZeroDiv   = 2;   // Division by zero that could not be
                   // removed by iterative refinement

//----------------------------------------------------------------------------
// Error messages depending on the error code.
//----------------------------------------------------------------------------
char* EvalErrMsg ( int Err )
{
  static char Msg[80] = "";

  if (Err != NoError) {
    char Hlp[60];

    switch (Err) {
      case ItFailed:
        sprintf(Hlp,"Maximal number of defect corrections (=%1d) exceeded",
                MaxStaggPrec); break;
      case ZeroDiv:
        strcpy(Hlp,"Division by zero occurred"); break;
      default:
        strcpy(Hlp,"Code not defined");
    }
    sprintf(Msg,"Error: %s!",Hlp);
  }
  return(Msg);
} // EvalErrMsg


//----------------------------------------------------------------------------
// Constructors, destructors, component access for the classes 'Staggered'
// and 'StaggArray
//----------------------------------------------------------------------------
Staggered::Staggered ( )                                        // Constructor
{                                                               //------------
  Resize(Val,0,MaxStaggPrec);
  for (int i = 0; i <= MaxStaggPrec; i++) Val[i] = 0.0;
  Err = 0.0;
}

Staggered::Staggered ( const Staggered& x )                // Copy constructor
{                                                          //-----------------
  Resize(Val,0,MaxStaggPrec);
  for (int i = 0; i <= MaxStaggPrec; i++) Val[i] = x.Val[i];
  Err = x.Err;
}

StaggArray::StaggArray ( )                                      // Constructor
{                                                               //------------
  Dim = 0;  SA = NULL;
}

StaggArray::StaggArray ( int n )                                // Constructor
{                                                               //------------
  if (n < 1) {
    cerr << "Lower bound is < 1 in 'StaggArray'-constructor!" << endl;
    exit(-1);
  }
  if ( !(SA = new Staggered[n]) ) {
    cerr << "Not enough memory for variable of type 'StaggArray'!" << endl;
    exit(-1);
  }

  Dim = n;
}

StaggArray::StaggArray ( const StaggArray& x )                   // Copy constructor
{                                                          //-----------------
  if ( (Dim < 1) || (x.SA == NULL) ) {
    cerr << "Illegal input parameter in 'StaggArray(StaggArray&)'!" << endl;
    exit(-1);
  }
  if ( !(SA = new Staggered[x.Dim]) ) {
    cerr << "Not enough memory for variable of type 'StaggArray'!" << endl;
    exit(-1);
  }

  Dim = x.Dim;
  for (int i = 0; i < Dim; i++) SA[i] = x.SA[i];
}

StaggArray::~StaggArray ( )                                      // Destructor
{                                                                //-----------
  Dim = 0;
  delete [] SA;
  SA = NULL;
}

Staggered& StaggArray::operator[] ( int i )                // Component access
{                                                          //-----------------
  if ( (i < 1) || (i > Dim) ) {
    cerr << "Illegal index access on variable of type 'StaggArray'!" << endl;
    exit(-1);
  }
  return( SA[--i] );
}

//----------------------------------------------------------------------------
// Since the intermediate results of operations for staggered data must be
// iteratively updated, it is necessary to store these data. For this
// purpose, a linear linked list is used. Its head is accessed via 'HeadPtr',
// whereas its actual entry is accessed via 'ActPtr'. 'FreePtr' is a pointer
// to a list of already allocated but actually unused entries. To prevent
// creation of garbage in memory, the entries of this list are used first
// when allocating a new intermediate result. All functions for the handling
// of the list of intermediate results are locally defined.
//----------------------------------------------------------------------------
static void InitList ( )            // Initialize list of intermediate results
{                                   //----------------------------------------
  DivByZero = 0;
  if (HeadPtr != NULL)       // Use list of freed entries
    FreePtr = HeadPtr;
  else                       // A list was not yet allocated
    FreePtr = NULL;
  HeadPtr = NULL;
  ActPtr = NULL;
}

static void ResetList ( )               // Reset error flag and set the actual
{                                       // pointer to the head of the list.
  DivByZero = 0;                        //------------------------------------
  ActPtr = HeadPtr;
}

static void AllocEntry ( IntResPtr& p )         // Allocate memory for a new
{                                               // entry. Use a previously
  if (FreePtr != NULL) {                        // created entry from the list
    p = FreePtr;                                // of freed entries if any.
    FreePtr = FreePtr->Next;                    //----------------------------
  }
  else {
    p = new IntRes;
    p->Next = NULL;
  }
}

void InitEntry ( real Approx )            // Get a new list entry and initial-
{                                         // ize its first staggered component
  IntResPtr p;                            // with the value 'Approx'.
                                          //----------------------------------
  AllocEntry(p);
  p->Entry.Val[0] = Approx;
  p->Next = NULL;
  if (HeadPtr == NULL)
    { HeadPtr = p;  ActPtr = p; }
  else
    { ActPtr->Next = p;  ActPtr = p; }
}

void UpdateError ( interval Error )              // Update the error component
{                                                // of the actual entry.
  ActPtr->Entry.Err = Error;                     //---------------------------
  ActPtr = ActPtr->Next;
}

void UpdateStaggComp ( int i )                    // Update the i-th staggered
{                                                 // component of all list
  ActPtr = HeadPtr;                               // entries by the midpoint
  while (ActPtr != NULL) {                        // of the error interval.
    ActPtr->Entry.Val[i] = mid(ActPtr->Entry.Err);//--------------------------
    ActPtr->Entry.Err = 0.0;
    ActPtr = ActPtr->Next;
  }
}

//----------------------------------------------------------------------------
// Assingment operators used to convert real input data to the staggered type
// and to assign a variables of type 'Staggered'. A real operand is assumed
// to be exact! It must not have been rounded. In particular, 'r' cannot be
//    - a real value such as 0.1 which is not exactly representable in the
//      internal binary format,
//    - rounded by conversion to the internal binary data format during
//      input,
//    - any arithmetic expression which cannot be evaluated exactly.
//----------------------------------------------------------------------------
Staggered& Staggered::operator= ( const real& r )
{
  int i;

  Val[0] = r;
  for (i = 1; i <= MaxStaggPrec; i++) Val[i] = 0.0;
  Err = 0.0;
  return *this;
}

Staggered& Staggered::operator= ( const Staggered& x )
{
  int i;

  for (i = 0; i <= MaxStaggPrec; i++) Val[i] = x.Val[i];
  Err = x.Err;
  return *this;
}

//----------------------------------------------------------------------------
// Arithmetic operators +, -, *, and / for operands of type 'Staggered'. A
// division by an interval containing zero will be avoided and is marked by
// setting the 'DivByZero' flag. If 'DivByZero' is set, all succeeding opera-
// tions are not executed. The evaluation is stopped at this point but may be
// restarted after updating a new staggered component by the midpoints of the
// error enclosures computed so far.
//----------------------------------------------------------------------------
// In the comments below the notations #*(...) and ##(...) are used to indi-
// cate the evaluation of the expression specified in round brackets using
// the exact scalar product. I.e. each component of the result is computed
// with only one final rounding. The symbol #* holds for rounding to nearest
// whereas ## holds for rounding to the smallest enclosing interval. An exact
// scalar product may be implemented using dotprecision accumulators.
//----------------------------------------------------------------------------
Staggered operator+ ( const Staggered& x , const Staggered& y )
{
  int           i;
  Staggered     z;
  idotprecision IAccu;

  if (!DivByZero) {
    if (InitFlag) {          // Initialize first staggered component
      z.Val[0] = x.Val[0] + y.Val[0];
      InitEntry(z.Val[0]);
    }
    else {
      z = ActPtr->Entry;   // Get actual values of z

      // Error:  dz := ##( x + y - z + dx + dy )
      //----------------------------------------
      IAccu  = x.Err;
      IAccu += y.Err;
      for (i = 0; i <= ActStaggPrec; i++) {
        IAccu += x.Val[i];
        IAccu += y.Val[i];
        IAccu -= z.Val[i];
      }
      rnd(IAccu,z.Err);
      UpdateError(z.Err);
    }
  } // if (!DivByZero)
  return(z);
} // operator+

Staggered operator- ( const Staggered& x , const Staggered& y )
{
  int           i;
  Staggered     z;
  idotprecision IAccu;

  if (!DivByZero) {
    if (InitFlag) {          // Initialize first staggered component
      z.Val[0] = x.Val[0] - y.Val[0];
      InitEntry(z.Val[0]);
    }
    else {
      z = ActPtr->Entry;   // Get actual values of z

      // Error:  dz := ##( x - y - z + dx - dy )
      //----------------------------------------
      IAccu  = x.Err;
      IAccu -= y.Err;
      for (i = 0; i <= ActStaggPrec; i++) {
        IAccu += x.Val[i];
        IAccu -= y.Val[i];
        IAccu -= z.Val[i];
      }
      rnd(IAccu,z.Err);
      UpdateError(z.Err);
    }
  } // if (!DivByZero)
  return(z);
} // operator-

Staggered operator* ( const Staggered& x , const Staggered& y )
{
  int           i, j;
  Staggered     z;
  idotprecision IAccu;

  if (!DivByZero) {
    if (InitFlag) {          // Initialize first staggered component
      z.Val[0] = x.Val[0] * y.Val[0];
      InitEntry(z.Val[0]);
    }
    else {
      z = ActPtr->Entry;   // Get actual values of z

      // Error:  dz := ##( x*y - z + y*dx + x*dy + dx*dy )
      //--------------------------------------------------
      IAccu = 0.0;
      accumulate(IAccu,x.Err,y.Err);
      for (i = 0; i <= ActStaggPrec; i++) {
        IAccu -= z.Val[i];
        for (j = 0; j <= ActStaggPrec; j++)
          accumulate(IAccu,x.Val[i],y.Val[j]);
      }
      for (i = 0; i <= ActStaggPrec; i++) {
        accumulate(IAccu,y.Val[i],x.Err);
        accumulate(IAccu,x.Val[i],y.Err);
      }
      rnd(IAccu,z.Err);
      UpdateError(z.Err);
    }
  } // if (!DivByZero)
  return(z);
} // operator*

Staggered operator/ ( const Staggered& x , const Staggered& y )
{
  int           i, j;
  Staggered     z;
  interval      num, denom;
  idotprecision IAccu;

  if (!DivByZero) {
    if (InitFlag)            // Initialize first staggered component
      if (y.Val[0] != 0.0) {
        z.Val[0] = x.Val[0] / y.Val[0];
        InitEntry(z.Val[0]);
      }
      else {
        DivByZero = 1;
        InitEntry(0.0);
      }
    else {
      z = ActPtr->Entry;   // Get actual values of z

      // Error:  dz := ##( x - z*y + dx - z*dy )  /  ##( y + dy )
      //---------------------------------------------------------
      IAccu = x.Err;
      for (i = 0; i <= ActStaggPrec; i++) {
        IAccu += x.Val[i];
        for (j = 0; j <= ActStaggPrec; j++)
          accumulate(IAccu,-z.Val[i],y.Val[j]);
      }
      for (i = 0; i <= ActStaggPrec; i++)
        accumulate(IAccu,-z.Val[i],y.Err);
      rnd(IAccu,num);
      IAccu = y.Err;
      for (i = 0; i <= ActStaggPrec; i++)
        IAccu += y.Val[i];
      rnd(IAccu,denom);

      if ( in(0.0,denom) )
        { z.Err = 0.0; DivByZero = 1; }
      else
        z.Err = num / denom;
      UpdateError(z.Err);
    }
  } // if (!DivByZero)
  return(z);
} // operator/

//----------------------------------------------------------------------------
// Arithmetic operators for different operands the one of type 'real' and the
// other one of type 'Staggered'. All these operators are implemented by
// first coercing both operands to 'Staggered' type and then calling the
// corresponding operators for the type 'Staggered'.
//----------------------------------------------------------------------------
Staggered operator+ ( const real& x, const Staggered& y )
{
  Staggered z;
  z = x; return(z+y);
}

Staggered operator+ ( const Staggered& x , const real& y )
{
  Staggered z;
  z = y; return(x+z);
}

Staggered operator- ( const real& x, const Staggered& y )
{
  Staggered z;
  z = x; return(z-y);
}

Staggered operator- ( const Staggered& x , const real& y )
{
  Staggered z;
  z = y; return(x-z);
}

Staggered operator* ( const real& x, const Staggered& y )
{
  Staggered z;
  z = x; return(z*y);
}

Staggered operator* ( const Staggered& x , const real& y )
{
  Staggered z;
  z = y; return(x*z);
}

Staggered operator/ ( const real& x, const Staggered& y )
{
  Staggered z;
  z = x; return(z/y);
}

Staggered operator/ ( const Staggered& x , const real& y )
{
  Staggered z;
  z = y; return(x/z);
}

//----------------------------------------------------------------------------
// Power function for integer exponents using the binary shift method. Thus
// the number of successive multiplications is reduced from n to log(2,n).
// Note: Since x^n is considered to be a monomial we define x^0 := 1.
//----------------------------------------------------------------------------
Staggered Power ( const Staggered& x, const int n )
{
  int        m;
  Staggered  p, z;

  if (!DivByZero) {
    p = 1.0;
    if (n != 0) {
      m = (n > 0) ? n : -n;

      p = 1.0; z = x;              // Binary shift method
      while (m > 0) {              //---------------------
        if (m % 2 == 1) p = p * z;
        m = m / 2;                 // integer division!
        if (m > 0) z = z * z;
      }

      if (n < 0) p = 1.0 / p;
    }
  } // if (!DivByZero)
  return p;                // Note: 'p' is undefined if 'DivByZero' is set
} // Power

//----------------------------------------------------------------------------
// Purpose: The function 'Eval()' may be used for the computation of an
//    enclosure of the value of a real arithmetic expression composed of the
//    operations +, -, *, /, and ^, where ^ denotes exponentiation by an
//    integer.
// Parameters:
//    In : 'f'        : a function of type 'Staggered' whose arguments are
//                      passed in an array of type 'Staggered'
//         'Arg'      : the real-valued arguments of 'f' stored as
//                      components of a real vector
//         'Eps'      : desired accuracy
//    Out: 'Approx'   : result computed with standard floating-point
//                      arithmetic
//         'Encl'     : verified enclosure of the result
//         'StaggPrec': number of corrections needed
//         'Err'      : error code
// Description:
//    The expression 'f' is evaluated for the real arguments which are stored
//    sequentially in the vector 'Arg'. Initially, the real arguments are
//    converted to arguments of type 'Staggered'. When 'f' is evaluated for
//    the first time using the special staggered arithmetic, the list of
//    intermediate results is initialized. Each time 'f' is evaluated again,
//    the error of every intermediate result is enclosed. The midpoints of
//    these enclosures are used to update the intermediate results. The
//    iteration is finished if the error of the last intermediate result
//    (= value of 'f') is less than the desired accuracy. Otherwise, the
//    iteration is halted after 'MaxStaggPrec' steps.
//----------------------------------------------------------------------------
void Eval ( Stagg_FctPtr f,
            rvector      Arg,
            real         Eps,
            real&        Approx,
            interval&    Encl,
            int&         StaggPrec,
            int&         Err )
{
  int           n = Ub(Arg)-Lb(Arg)+1;
  int           i, Success;
  StaggArray    x(n);
  Staggered     StaggRes;
  idotprecision IAccu;

  for (i = 1; i <= n; i++) x[i] = Arg[i];    // Initialize arguments

  InitList();                    // Initialize list for intermediate results.
  InitFlag = 1;                  // So far, no staggered corrections computed.
  StaggRes = f(x);               // Compute first staggered component.

  // In general, the first component of 'StaggRes' will now hold the usual
  // floating-point approximation, except when a division by zero occurred.
  // In the latter case, a signaling NaN (= Not a Number, for more details
  // see the C-XSC reference) is returned.
  //-----------------------------------------------------------------------
  if (DivByZero)
    Approx = SignalingNaN;
  else
    Approx = StaggRes.Val[0];
                                                 // Initial approximations are
  InitFlag = 0; ActStaggPrec = 0;                // already computed.
  do {
    ResetList();               // Compute new enclosures of the absolute error
    StaggRes = f(x);           // of any intermediate result.

    // Compute an enclosure of f(x) by ##( StaggRes + StaggRes.Err)
    //-------------------------------------------------------------
    IAccu = StaggRes.Err;
    for (i = 0; i <= ActStaggPrec; i++) IAccu += StaggRes.Val[i];
    rnd(IAccu,Encl);

    Success = ( !DivByZero && (RelDiam(Encl) <= Eps) );

    // Increment actual staggered precision and store the next component of
    // the intermediate results by updating the list of intermediate results.
    //-----------------------------------------------------------------------
    if ( !Success && (ActStaggPrec <= MaxStaggPrec) ) {
      ActStaggPrec++;
      UpdateStaggComp(ActStaggPrec);
    }

  } while ( !(Success || (ActStaggPrec == MaxStaggPrec)) );
  StaggPrec = ActStaggPrec;

  if (Success)             // Set error code
    Err = NoError;
  else if (DivByZero)
    Err = ZeroDiv;
  else
    Err = ItFailed;
} // Eval




