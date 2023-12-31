/// @addtogroup vecttst
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file floatMatrixAlgorithmsTest.cpp
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details. 

#include <fstream>
#include <sstream>
#include <iomanip>
#include "capd/intervals/minmax_interval.h"
#include "capd/vectalg/Norm.hpp"
#include "capd/vectalg/lib.h"
#include "capd/rounding/DoubleRounding.h"
#include "capd/vectalg/Matrix.hpp"

using namespace capd;
using namespace capd::matrixAlgorithms;

void test_gauss(void)
{
   std::ifstream inp("DMatrix.txt");

   std::cout << "\n\nGauss test\n";
   std::cout << "==========\n";
   if(!inp) std::cout << "Failed to open inp\n";

   DVector v;
   DMatrix a,b;

   inp >> v;
   inp >> a >> b;

   IVector iv(v);
   IMatrix ia(a);

   DVector r = gauss(a,v);

   capd::rounding::DoubleRounding::roundUp();
   IVector ir = gauss(ia,iv);
   capd::rounding::DoubleRounding::roundNearest();

   std::cout << "v=" << v << " iv=" << iv << std::endl;
   std::cout << "a=" << a << " ia=" << ia << std::endl;
   std::cout << "r=" << r << std::endl;
   std::cout << "ir=" << ir << std::endl;
   for(int i=0;i<(int)ir.dimension();i++)
   {
      std::cout << "     " << diam(ir[i]) << std::endl;
   }
}

bool assertion_failed(const char * message = ""){
   std::cerr << "Assertion FAILED : " << message << std::endl;
   exit(1);
}

bool assert_true(bool condition, const char * message = ""){
  if(condition)
    return true; 
  return assertion_failed(message);
}

bool assert_subset(const interval & a, const interval & b, double eps = -1.0, const char * message = ""){
  if(a.subset(b) and ((eps < 0) or (abs(b-a) <= eps)))
    return true;
  std::cout << "\n abs(b-a) : " << abs(b-a) << "  eps : " << eps << std::endl;
  std::cout << std::setprecision(17) << "\n a = " << a << "\n b = " << b << std::endl;
  return assertion_failed(message);  
}
bool assert_vectorSubset(
   const IVector & a, const IVector & b, 
   double eps = -1.0, const char *  = ""
){
  if(a.dimension() != b.dimension())
    return assertion_failed("Not equal dimensions of vectors.");
  IVector::const_iterator ia = a.begin(),  ib = b.begin(), stop = a.end();
  bool res= true;
  while(ia != stop){
    //std::cout << *ia <<"  " << *ib << std::endl;
    res = res && assert_subset(*ia++, *ib++, eps, "Vector subset.");
  }
  return res;
}
bool assert_matrixSubset(
   const IMatrix & a, const IMatrix & b, 
   double eps = -1.0, const char *  = ""
){
  if(a.numberOfColumns() != b.numberOfColumns() or a.numberOfRows()!=b.numberOfRows())
    return assertion_failed("Not equal dimensions of matrices.");
  IMatrix::const_iterator ia = a.begin(),  ib = b.begin(), stop = a.end();
  bool res= true;
  while(ia != stop){
  //  std::cout << *ia <<"  " << *ib << std::endl;
    res = res && assert_subset(*ia++, *ib++, eps, "Matrix subset.");
  }
  return res;
}

void QRdecomp_test(void)
{
   std::cout << "\n\nQRdecomp test\n";
   std::cout << "=============\n";

   interval a[] = {
     7,63,0, 
     2,18,10,
     3,30,0
   };
   double tolerance = 1.0e-13; // tolerance
   IMatrix ia(3,3,a);
   IMatrix iq(ia.numberOfRows(),ia.numberOfColumns()), ir(ia.numberOfRows(),ia.numberOfColumns());

   std::cout << "ia=" << ia << "\n";

   QR_decompose(ia,iq,ir);
   std::cout << "iq=" << iq << "\n";
   std::cout << "ir=" << ir << "\n";
   std::cout << "iq*ir=" << iq*ir << "\n\n";
   assert_matrixSubset(ia, iq*ir, 1.e-12, "QR decomposition");
   interval zero(0.0), one(1.0);
   std::cout << "\n orthogonality test \n" ;
   for(int i=0;i<(int)iq.numberOfRows();i++)
   {
      for(int j=0;j<(int)iq.numberOfRows();j++)
      {
         interval sp=iq.row(i)*iq.row(j);
         std::cout  << sp << " ";
         assert_subset(i==j?one:zero, sp, tolerance);
      }
      std::cout << std::endl;
   }
}

void QRdecomp_test2(void)
{
   std::cout << "\n\nQRdecomp test\n";
   std::cout << "=============\n";

   interval a[] = {
     7,63,0, 
     2,18,10,
     3,30,0
   };
   
   // TODO: check tolerance changes
   double tolerance = 1.0e-12; // tolerance
   IMatrix ia(3,3,a), A(ia), B(ia);
   IMatrix iq(ia.numberOfRows(),ia.numberOfColumns()),
           ir(ia.numberOfRows(),ia.numberOfColumns());
   IVector v(3);
   v[0]=5; v[1]=4; v[2]=3;
   IVector sizes(3);
   ZVector perm(3);
   std::cout << "ia=" << ia << "\n";

   QRdecomposeWithPivoting(A,v,iq,ir,sizes,perm);
   std::cout << "iq=" << iq << "\n";
   std::cout << "ir=" << ir << "\n";
   std::cout << "iq*ir=" << iq*ir << "\n\n";
   assert_matrixSubset(ia, iq*ir, 3.e-12, "QR decomposition");
   interval zero(0.0), one(1.0);
   std::cout << "\n orthogonality test \n" ;
   for(int i=0;i<(int)iq.numberOfRows();i++)
   {
      for(int j=0;j<(int)iq.numberOfRows();j++)
      {
         interval sp=iq.row(i)*iq.row(j);
         std::cout  << sp << " ";
         assert_subset(i==j?one:zero, sp, tolerance);
      }
      std::cout << std::endl;
   }
  orthonormalize(B, v);
  std::cout << "\n B = " << B;
  std::cout << "\n orthogonality test \n" ;
     for(int i=0;i<(int)B.numberOfRows();i++)
   {
      for(int j=0;j<(int)B.numberOfRows();j++)
      {
         interval sp=B.row(i)*B.row(j);
         std::cout  << sp << " ";
         assert_subset(i==j?one:zero, sp, tolerance);
      }
      std::cout << std::endl;
   }
  
}

  void diagonalize_test(void)
{
   std::cout << "\n\ndiagonalize test\n";
   std::cout << "================\n";
  // std::ifstream inp("diagdane.txt", std::ios::in);
  // if(!inp) std::cout << "Failed to open inp\n";
   double data[] = {7,63,3,  63,18,10,   3,10,0};
   DMatrix a(3, 3, data);
  
   IMatrix ia(a), id(a.numberOfRows(),a.numberOfColumns());

   std::cout << "a=" << a << "\n";
   std::cout << "ia=" << ia << "\n";

   symMatrixDiagonalize(ia, id,  1.0e-12);
   for(int i=0;i<(int)id.numberOfRows();i++)
   {
      std::cout  << "id[" << i << "]=" << id.row(i) << "\n";
   }
   std::cout  << "Spectral radius ia = " << spectralRadiusOfSymMatrix(ia) << "\n";
   std::cout  << "Spectral radius id = " << spectralRadiusOfSymMatrix(id) << "\n";
}
void norm_test(void)
{
   std::cout << "\n\nnorm test\n";
   std::cout << "============================\n";

   double data[] = {7,63,3,   12,-18,10,    11,-8,0};
   //{{170.5, -169.5, 330.}, {-169.5, 170.5, -330.}, {330., -330., 670.}} 
   DMatrix a(3, 3, data);
   IMatrix ia(a);

   IEuclNorm en;
   IMaxNorm mn;
   ISumNorm sn;
   IEuclLNorm eln;
   IMaxLNorm mln;
   ISumLNorm sln;
   DVector re(3), im(3);
 
   std::cout << "Euclidean norm = " << std::setprecision(16) << en(ia) << " \n";
                       // 66.043028616955155 alglib
   assert_subset(interval(66.043028616955157), en(ia), 1.0e-10);
   std::cout << "Maximum norm = " << mn(ia) << " \n";
   assert_subset(interval(73), mn(ia), 0.0);
   std::cout << "Sum norm = " << sn(ia) << " \n";
   assert_subset(interval(89), sn(ia), 0.0);
   std::cout << "Euclidean logarithmic norm = " << eln(ia) << " \n";
   assert_subset(eln(ia), interval(35.148869407, 35.148869408), 1.0e-9);
   std::cout << "Maximum logarithmic norm = " << mln(ia) << " \n";
   assert_subset(interval(73), mln(ia), 0.0);
   std::cout << "Sum logarithmic norm = " << sln(ia) << " \n";
   assert_subset(interval(53), sln(ia), 0.0);
   
}
void norm_test2(void)
{
   std::cout << "\n\nnorm test 2\n";
   std::cout << "============================\n";

   double data[] =
   //{170.5, -169.5, 330., -169.5, 170.5, -330., 330., -330., 670.};
   {-977, 983, -2020, 983, -977, 2020, -2020, 2020, -3980}; 
   DMatrix a(3, 3, data);
   IMatrix ia(a);
   
   IEuclNorm en;
   IMaxNorm mn;
   ISumNorm sn;
   IEuclLNorm eln;
   IMaxLNorm mln;
   ISumLNorm sln;
  
   std::cout << "Euclidean norm = " << std::setprecision(16) << en(ia) << " \n";
   assert_subset(6000, en(ia), 1.0e-9);
   std::cout << "Maximum norm = " << mn(ia) << " \n";
   assert_subset(interval(8020), mn(ia), 0.0);
   std::cout << "Sum norm = " << sn(ia) << " \n";
   assert_subset(interval(8020), sn(ia), 0.0);
   std::cout << "Euclidean logarithmic norm = " << eln(ia) << " \n";
   assert_subset( interval(60), eln(ia), 1.0e-9);
   std::cout << "Maximum logarithmic norm = " << mln(ia) << " \n";
   assert_subset(interval(2026), mln(ia), 0.0);
   std::cout << "Sum logarithmic norm = " << sln(ia) << " \n";
   assert_subset(interval(2026), sln(ia), 0.0);
   
}
void norm_test3(void)
{
   std::cout << "\n\nnorm test 3\n";
   std::cout << "============================\n";

   double data[] =
   {170.5, -169.5, 330., -169.5, 170.5, -330., 330., -330., 670.};
   //{-977, 983, -2020, 983, -977, 2020, -2020, 2020, -3980}; 
   DMatrix a(3, 3, data);
   IMatrix ia(a);
   IEuclNorm en;
   IMaxNorm mn;
   ISumNorm sn;
   IEuclLNorm eln;
   IMaxLNorm mln;
   ISumLNorm sln;
  
   std::cout << "Euclidean norm = " << std::setprecision(16) << en(ia) << " \n";
   assert_subset(1000, en(ia), 1.0e-10);
   std::cout << "Maximum norm = " << mn(ia) << " \n";
   assert_subset(interval(1330), mn(ia), 0.0);
   std::cout << "Sum norm = " << sn(ia) << " \n";
   assert_subset(interval(1330), sn(ia), 0.0);
   std::cout << "Euclidean logarithmic norm = " << eln(ia) << " \n";
   assert_subset(interval(1000), eln(ia),  1.0e-9);
   std::cout << "Maximum logarithmic norm = " << mln(ia) << " \n";
   assert_subset(interval(1330), mln(ia), 0.0);
   std::cout << "Sum logarithmic norm = " << sln(ia) << " \n";
   assert_subset(interval(1330), sln(ia), 0.0);
   
}
void ztest_gauss(void)
{
   std::cout << "\n\n(z) Gauss test\n";
   std::cout << "==============\n";

   IMatrix ia(3,3);
   IVector zero(3);
   
   ia[0][0]=-0.237304;
   ia[0][1]= -1.18286;
   ia[0][2]=8.49793e-322;

   ia[1][0]=1.18286;
   ia[1][1]=-0.000731345;
   ia[1][2]=-5.7114e-321;

   ia[2][0]=-6.81317e-321;
   ia[2][1]=-1.34386e-321;
   ia[2][2]=1.0;
  { 
   IVector iv(3);
   iv[0]=interval(-7.45037e-16,-5.4116e-16);
   iv[1]=interval(-2.12588e-15,-9.94006e-16);
   iv[2]=interval(-4.44089e-21,4.44089e-21);

   IVector r = gauss(ia,iv);
   IVector res = ia*r - iv;

   std::cout << "ia*r - iv = " << res << "\n";
   std::cout << "iv = " << iv << "\n";
   std::cout << "r = " << r << "\n";
   
   assert_vectorSubset(zero, res, 1.0e-14);
  } 
  {
   IVector iv(3);
   iv[0]=5.0;
   iv[1]=0.0;
   iv[2]=0.0;
   
   IVector r = gauss(ia,iv);
   IVector res = ia*r - iv;

   std::cout << "ia*r - iv = " << res << "\n";
   std::cout << "iv = " << iv << "\n";
   std::cout << "r = " << r << "\n";
   assert_vectorSubset(zero, res, 1.0e-14);
  } 
}

int main(int , char* [])
{
   //init_fpunit();

   try{
      if((CAPD_DEFAULT_DIMENSION !=0) && (CAPD_DEFAULT_DIMENSION!=3)) throw std::runtime_error(
         " Wrong dimension!  \n dimension in file vectalg/dimension.h should be 3.\n Change and recompile!"
         );
        QRdecomp_test2();
        diagonalize_test();
        norm_test();
        norm_test2();
        norm_test3();
//      matrix_test();
        ztest_gauss();
      // test_gauss();
   }catch(std::exception& e)
   {
      std::cout << e.what();
   }
   catch(...)
   {
      std::cout << "Unknown exception!\n";
   }
   capd::rounding::DoubleRounding::roundNearest();
   return 0;
}

/// @}
