const char cpp_file[] ="\
#include \"%s.hpp\"\n\
\n\
template class %s::%s<double>;\n\
template class %s::%s<capd::interval>;\n\
\n\
\n";

// #####################################################################

const char hpp_header[] ="\
#ifndef __%s_%s_HPP__\n\
#define __%s_%s_HPP__\n\n\
#include \"%s.h\"\n\
#include \"capd/autodiff/eval.hpp\"\n\
#include \"capd/vectalg/Vector.hpp\"\n\
#include \"capd/vectalg/Matrix.hpp\"\n\
#include \"capd/vectalg/Norm.hpp\"\n\
#include \"capd/matrixAlgorithms/floatMatrixAlgorithms.hpp\"\n\
#include \"capd/map/Function.hpp\"\n\
#include \"capd/dynsys/TaylorHOE.hpp\"\n\
#include \"capd/poincare/TimeMap.hpp\"\n\
#include \"capd/poincare/PoincareMap.hpp\"\n\
\n\
#define CAPD__TEMPLATE_MAP_CLASSNAME__ %s\n\
\n\
namespace %s{\n\
\n";

// #####################################################################


const char hpp_endheader[] ="\n\
} /* namespace */\n\
\n\
#undef CAPD__TEMPLATE_MAP_CLASSNAME__\n\
\n\
#endif\n\n\
";

// #####################################################################

const char header[] ="\
#ifndef __%s_%s_H__\n\
#define __%s_%s_H__\n\
\n\
#include <algorithm>\n\
#include \"capd/capdlib.h\"\n\
\n\
#define CAPD__TEMPLATE_MAP_CLASSNAME__ %s\n\
\n\
namespace %s{\n\
\n\
template<class Scalar>\n\
class CAPD__TEMPLATE_MAP_CLASSNAME__ {\n\
public:\n\
  const static int dimIn = %d;\n\
  const static int dimOut = %d;\n\
  const static int jetSize = %d;\n\
  const static int timeJetSize = %d;\n\
  const static int maxOrder = %d;\n\
  const static int dataSize = %d;\n\
\n\
  typedef capd::vectalg::Matrix<Scalar,dimOut,dimIn> MatrixType;\n\
  typedef typename MatrixType::RowVectorType VectorType;\n\
  typedef typename MatrixType::ColumnVectorType ImageVectorType;\n\
  typedef typename VectorType::ScalarType ScalarType;\n\
  typedef capd::map::Function<VectorType> FunctionType;\n\
\n\
  CAPD__TEMPLATE_MAP_CLASSNAME__ (){\n\
    std::fill(data,data+dataSize,capd::TypeTraits<ScalarType>::zero());\n\
    data[dimIn*timeJetSize+1]=capd::TypeTraits<ScalarType>::one(); // set dt/dt=1\n\
    setConstants();\n\
  }\n\
  void setConstants();\n\
\n\
  ImageVectorType operator()(const VectorType& u) const;\n\
  ImageVectorType operator()(ScalarType t, const VectorType& u) const{\n\
    this->setCurrentTime(t);\n\
    return (*this)(u);\n\
  }\n\
\n\
  ImageVectorType operator()(const VectorType& u, MatrixType& out_derivative) const;\n\
  ImageVectorType operator()(ScalarType t, const VectorType& u, MatrixType& out_derivative) const{\n\
    this->setCurrentTime(t);\n\
    return (*this)(u,out_derivative);\n\
  }\n\
\n\
  MatrixType operator[](const VectorType& u) const{\n\
    MatrixType der;\n\
    (*this)(u,der);\n\
    return der;\n\
  }\n\
\n\
  MatrixType derivative(ScalarType t, const VectorType& u) const{\n\
    this->setCurrentTime(t);\n\
    return (*this)[u];\n\
  }\n\
\n\
  void computeODECoefficients(VectorType coeffs[], unsigned order) const;\n\
  void computeODECoefficients(VectorType coeffs[], MatrixType dCoeffs[], unsigned order) const;\n\
\n\
  unsigned degree() const { return 1; }\n\
  unsigned dimension() const {return dimIn;}\n\
  unsigned imageDimension() const {return dimOut; }\n\
\n\
  void setOrder(unsigned r){ if(r>=maxOrder) throw std::logic_error(\"Cannot change order in generated class.\");}\n\
  unsigned getOrder() const { throw \"getOrder\"; }\n\
\n\
  void setCurrentTime(const ScalarType& t) const { data[dimIn*timeJetSize]=t; }\n\
  const ScalarType& getCurrentTime() const { return data[dimIn*timeJetSize]; }\n\
  void differentiateTime() const {}\n\
\n\
  void setParameter(unsigned d, const ScalarType& value) { data[(dimIn+d+1)*timeJetSize] = value; }\n\
  void setParameters(const VectorType& values){\n\
     for(int i=0;i<values.dimension();++i) setParameter(i,values[i]);\n\
  }\n\
\n\
private:\n\
  mutable ScalarType data[dataSize];\n\
  const static int pos[];\n\
};\n\
\n\
typedef CAPD__TEMPLATE_MAP_CLASSNAME__ <double>::VectorType DVector;\n\
typedef CAPD__TEMPLATE_MAP_CLASSNAME__ <double>::MatrixType DMatrix;\n\
typedef CAPD__TEMPLATE_MAP_CLASSNAME__ <double> D%s;\n\
typedef CAPD__TEMPLATE_MAP_CLASSNAME__ <capd::interval>::VectorType IVector;\n\
typedef CAPD__TEMPLATE_MAP_CLASSNAME__ <capd::interval>::MatrixType IMatrix;\n\
typedef CAPD__TEMPLATE_MAP_CLASSNAME__ <capd::interval> I%s;\n\
\n\
typedef capd::map::Function<DVector> DFunction;\n\
typedef capd::map::Function<IVector> IFunction;\n\
\n\
typedef capd::dynsys::BasicTaylor< CAPD__TEMPLATE_MAP_CLASSNAME__ <double> > DTaylor;\n\
typedef capd::dynsys::TaylorHOE< CAPD__TEMPLATE_MAP_CLASSNAME__ <capd::interval> > ITaylor;\n\
\n\
typedef capd::poincare::TimeMap< DTaylor > DTimeMap;\n\
typedef capd::poincare::TimeMap< ITaylor > ITimeMap;\n\
typedef capd::poincare::PoincareMap< DTaylor > DPoincareMap;\n\
typedef capd::poincare::PoincareMap< ITaylor > IPoincareMap;\n\
\n\
\n\
} /* namespace */\n\
\n\
#undef CAPD__TEMPLATE_MAP_CLASSNAME__\n\
\n\
#endif\n\n\
";

// #####################################################################

const char evalOperatorBegin[] ="\
template<class Scalar>\n\
typename CAPD__TEMPLATE_MAP_CLASSNAME__ <Scalar>::ImageVectorType\nCAPD__TEMPLATE_MAP_CLASSNAME__ <Scalar>::operator()(const VectorType& u) const{\n\
  ScalarType* p = data;\n\
  for(typename VectorType::const_iterator b=u.begin();b!=u.end();++b,p+=timeJetSize)\n\
    *p=*b;\n\
  \n\
";

// #####################################################################

const char evalOperatorEnd[] ="\n\
  ImageVectorType result;\n\
  for(int i=0;i<result.dimension();++i)\n\
    result[i] = data[pos[i]];\n\
  return result;\n\
}\n\
\n\
";

// #####################################################################

const char evalOperatorDerBegin[] ="\n\
template<class Scalar>\n\
typename CAPD__TEMPLATE_MAP_CLASSNAME__ <Scalar>::ImageVectorType\nCAPD__TEMPLATE_MAP_CLASSNAME__ <Scalar>::operator()(const VectorType& u, MatrixType& der) const{\n\
  ScalarType* p = data;\n\
  int i,j;\n\
  for(i=0;i<u.dimension();++i,p+=timeJetSize){\n\
    *p=u[i];\n\
    ScalarType* q=p+maxOrder;\n\
    for(j=0;j<u.dimension();++j,q+=maxOrder)\n\
      *q = (i==j) ? capd::TypeTraits<ScalarType>::one() : capd::TypeTraits<ScalarType>::zero();\n\
  }\n\
  \n\
";

// #####################################################################

const char evalOperatorDerEnd[] ="\n\
  ImageVectorType result(this->imageDimension());\n\
  for(i=0;i<this->imageDimension();++i){\n\
    ScalarType* p = data+pos[i];\n\
    result[i] = *p;\n\
    p+=maxOrder;\n\
    for(j=0;j<this->dimension();++j,p+=maxOrder)\n\
      der(i+1,j+1) = *p;\n\
  }\n\
  return result;\n\
}\n\
\n\
";

// #####################################################################

const char evalODEBegin[] ="\
template<class Scalar>\n\
void CAPD__TEMPLATE_MAP_CLASSNAME__ <Scalar>::computeODECoefficients(VectorType u[], unsigned order) const{\n\
  ScalarType* p = data;\n\
  int i;\n\
  for(typename VectorType::const_iterator b=u[0].begin();b!=u[0].end();++b,p+=timeJetSize)\n\
    *p=*b;\n\
\n\
  for(unsigned k=0;k<order;++k){\n\
";

// #####################################################################

const char evalODEEnd[] ="\n\
    for(i=0, p=data+k+1;i<dimIn;++i, p+=timeJetSize)\n\
      u[k+1][i] = *p = data[pos[i]+k]/(k+1);\n\
  }\n\
}\n\
";

// #####################################################################

const char evalODEDerBegin[] ="\
template<class Scalar>\n\
void CAPD__TEMPLATE_MAP_CLASSNAME__ <Scalar>::computeODECoefficients(VectorType u[], MatrixType der[], unsigned order) const{\n\
  ScalarType* p = data;\n\
  typename VectorType::const_iterator b=u[0].begin();\n\
  typename MatrixType::const_iterator i=der[0].begin(), e = der[0].begin() + dimIn;\n\
  for(;b!=u[0].end();++b,p+=timeJetSize){\n\
    *p=*b;\n\
    ScalarType* q = p+maxOrder;\n\
    while(i!=e) { *q=*i; q+=maxOrder; ++i; }\n\
    e += dimIn;\n\
  }\n\
\n\
  for(unsigned k=0;k<order;++k){\n\
";

// #####################################################################

const char evalODEDerEnd[] ="\n\
    int j;\n\
    for(j=0, p=data+k+1;j<dimIn;++j, p+=timeJetSize){\n\
      u[k+1][j] = *p = data[pos[j]+k]/(k+1);\n\
      for(int c=0, r=maxOrder; c<dimIn; ++c, r+=maxOrder)\n\
        der[k+1](j+1,c+1) = *(p+r) = data[pos[j]+r+k]/(k+1);\n\
    }\n\
  }\n\
}\n\
";
