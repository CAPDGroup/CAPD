/// @addtogroup map
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file Map.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_MAP_MAP_CODETRANSLATION_HPP_
#define _CAPD_MAP_MAP_CODETRANSLATION_HPP_

#include "capd/map/Map.h"
#include "capd/map/Map_headerTemplate.h"
#include <cstdio>
#include <sstream>

namespace capd{
namespace map{

template<typename MatrixT>
void Map<MatrixT>::codeTranslation(const char className[], const char userNamespace[], const char relativePath[], size_type maxTimeOrder) const{
  using namespace capd::autodiff;
  generateHeaderFile(className,userNamespace,relativePath,maxTimeOrder);
  generateHppFile(className,userNamespace,relativePath,maxTimeOrder);
  generateCppFile(className,userNamespace,relativePath,maxTimeOrder);
}

// #####################################################################

template<typename MatrixT>
void Map<MatrixT>::generateHeaderFile(const char className[], const char userNamespace[], const char relativePath[], size_type maxTimeOrder) const{
  std::ostringstream fullPath;

  fullPath << relativePath << "/" << className << ".h";
  FILE* file = fopen(fullPath.str().c_str(),"w");
  fprintf(file,header,
              userNamespace,className,
              userNamespace,className,
              className,userNamespace,
              (int)this->dimension(),
              (int)this->imageDimension(),
              (int)(this->m_dag.jetSize()),
              (int)(maxTimeOrder*this->m_dag.jetSize()),
              (int)maxTimeOrder,
              (int)(maxTimeOrder*this->m_dag.numberOfNodes()*this->m_dag.timeJetSize()),
              className, className
         );
  fclose(file);
}

// #####################################################################

template<typename MatrixT>
void Map<MatrixT>::generateCppFile(const char className[], const char userNamespace[], const char relativePath[], size_type /*maxTimeOrder*/) const{
  std::ostringstream fullPath;

  fullPath << relativePath << "/" << className << ".cpp";
  FILE* file = fopen(fullPath.str().c_str(),"w");
  fprintf(file,cpp_file,
              className,
              userNamespace,className,
              userNamespace,className
         );
  fclose(file);
}

// #####################################################################

template<typename MatrixT>
void Map<MatrixT>::generateHppFile(const char className[], const char userNamespace[], const char relativePath[], size_type maxTimeOrder) const{
  using namespace capd::autodiff;

  std::ostringstream fullPath;
  fullPath << relativePath << "/" << className << ".hpp";
  FILE* file = fopen(fullPath.str().c_str(),"w");
  fprintf(file,hpp_header,
      className,userNamespace,
      className,userNamespace,
      className,className,
      userNamespace
      );

  int i;
  int timeJetSize = (int)(maxTimeOrder*this->m_dag.jetSize());
  fprintf(file,"template<class Scalar>\nconst int CAPD__TEMPLATE_MAP_CLASSNAME__ <Scalar>::pos[] ={%d",(int)this->m_pos[0]*timeJetSize);
  for(i=1;i<(int)this->imageDimension();++i)
    fprintf(file,",%d",(int)this->m_pos[i]*timeJetSize);
  fprintf(file,"};\n\n");

  // ############# setConstants() #########################
  fprintf(file,"template<class Scalar>\nvoid CAPD__TEMPLATE_MAP_CLASSNAME__ <Scalar>::setConstants(){\n");
  for(i=0;i<(int)this->m_fullGraph.size();++i){
    if(this->m_fullGraph[i].isConst)
      fprintf(file,"  data[%d]=%.17f;\n",(int)i*timeJetSize,toDouble(leftBound(this->m_dag(VarNo(i),CoeffNo(0)))));
  }
  fprintf(file,"}\n\n");

  // ############# operator(const VectorType& u) #########################
  fprintf(file,"%s",evalOperatorBegin);
  for(i=0;i<(int)this->m_evalPath.size();++i){
    fprintf(file,"  capd::autodiff::%s::evalC0HomogenousPolynomial(data+%d,data+%d,data+%d);\n",
        this->m_nodes[i]->name(),
        this->m_evalPath[i].left*timeJetSize,
        this->m_evalPath[i].right*timeJetSize,
        this->m_evalPath[i].result*timeJetSize
        );
  }
  fprintf(file,"%s",evalOperatorEnd);

  // ####### operator(const VectorType& u, MatrixType& out_der) #############
  fprintf(file,"%s",evalOperatorDerBegin);
  for(i=0;i<(int)this->m_evalPath.size();++i){
    fprintf(file,"  capd::autodiff::%s::evalC0HomogenousPolynomial(data+%d,data+%d,data+%d);\n",
        this->m_nodes[i]->name(),
        this->m_evalPath[i].left*timeJetSize,
        this->m_evalPath[i].right*timeJetSize,
        this->m_evalPath[i].result*timeJetSize
        );
    fprintf(file,"  capd::autodiff::%s::evalHomogenousPolynomial(1,data+%d,data+%d,data+%d,dimIn,maxOrder);\n",
        this->m_nodes[i]->name(),
        this->m_evalPath[i].left*timeJetSize,
        this->m_evalPath[i].right*timeJetSize,
        this->m_evalPath[i].result*timeJetSize
        );
  }
  fprintf(file,"%s",evalOperatorDerEnd);

  // ####### computeODECoefficients(VectorType coeffs[], unsigned order) #############
  fprintf(file,"%s",evalODEBegin);
  for(i=0;i<(int)this->m_evalPath.size();++i)
    fprintf(file,"    capd::autodiff::%s::evalC0(data+%d,data+%d,data+%d,k);\n",
        this->m_nodes[i]->name(),
        this->m_evalPath[i].left*timeJetSize,
        this->m_evalPath[i].right*timeJetSize,
        this->m_evalPath[i].result*timeJetSize
        );
  fprintf(file,"%s",evalODEEnd);

  // ####### computeODECoefficients(VectorType c[], MatrixType v[], unsigned order) #############
  fprintf(file,"%s",evalODEDerBegin);
  for(i=0;i<(int)this->m_evalPath.size();++i)
    fprintf(file,"    capd::autodiff::%s::eval(1,data+%d,data+%d,data+%d,dimIn,maxOrder,k);\n",
        this->m_nodes[i]->name(),
        this->m_evalPath[i].left*timeJetSize,
        this->m_evalPath[i].right*timeJetSize,
        this->m_evalPath[i].result*timeJetSize
        );
  fprintf(file,"%s",evalODEDerEnd);
  fprintf(file,"%s",hpp_endheader);
  fclose(file);
}

// #####################################################################

}} // namespace capd::map

#endif // _CAPD_MAP_MAP_HPP_

/// @}
