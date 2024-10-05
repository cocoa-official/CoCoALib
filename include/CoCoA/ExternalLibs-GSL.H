#ifndef CoCoA_gsl_H
#define CoCoA_gsl_H

//   Copyright (c)  2011  Bruno Simoes

//   This file is part of the source of CoCoALib, the CoCoA Library.
//
//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.

#include "CoCoA/PREPROCESSOR_DEFNS.H"

#ifdef CoCoA_WITH_GSL

#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>

namespace CoCoA
{

  class matrix;          // fwd decl -- defined in matrix.H
  class ring;          // fwd decl -- defined in ring.H
  class RingElem;          // fwd decl -- defined in ring.H
  class ConstMatrixView; // fwd decl -- defined in MatrixView.H

  matrix NewMatrix(const ring& R, gsl_matrix* m);
  matrix NewMatrix(const ring& R, gsl_matrix* m, bool transpose);
  
  matrix NewMatrix(const ring& R, gsl_vector* m);
  gsl_matrix* NewMatrixGSL(ConstMatrixView M);
  gsl_matrix* NewMatrixGSL(matrix& M);
  void UpdateMatrix(matrix& A, gsl_matrix* m);
  
  std::vector<RingElem> NewVector(const ring& R, gsl_vector* m);
  gsl_vector* NewVectorGSL(std::vector<RingElem>& v);

  // Print tools
  void PrintVector(gsl_vector *m);
  void PrintMatrix(gsl_matrix *m);
  
  // SVD
  std::vector<matrix> GslSVD(ConstMatrixView A);
  std::vector<RingElem> GslSingularValues(ConstMatrixView A);
  std::vector<RingElem> GslSolveSVD(ConstMatrixView A, std::vector<RingElem>& b);
  
  // LU
  matrix GslLU(ConstMatrixView A);
  std::vector<RingElem> GslSolverLU(ConstMatrixView A, std::vector<RingElem>& b);
  
  // Cholesky Decomposition
  matrix GslCholeskyDecomposition(ConstMatrixView A);

  // QR
  std::vector<matrix> GslQR(ConstMatrixView A);
  std::vector<matrix> GslQRPT(ConstMatrixView A);
  std::vector<RingElem> GslSolverQR(ConstMatrixView A, std::vector<RingElem>& b);

  void GslQRDecomposition(ConstMatrixView A, matrix& qr, std::vector<RingElem>& tau);
    
  // B-Decomp
  void GslBidiagDecomposition(ConstMatrixView A, matrix& D, std::vector<RingElem>& tau_U, std::vector<RingElem>& tau_V);
  void GslBidiagDecompUnpack(ConstMatrixView A, matrix& U, matrix& V, std::vector<RingElem>& diagonal, std::vector<RingElem>& superdiagonal);
  
}
#endif

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/ExternalLibs-GSL.H,v 1.9 2022/02/18 14:11:45 abbott Exp $
// $Log: ExternalLibs-GSL.H,v $
// Revision 1.9  2022/02/18 14:11:45  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.8  2017/03/13 12:17:08  abbott
// Summary: Include PREPROCESSOR_DEFNS.H (in header files which need the info)
//
// Revision 1.7  2014/09/02 15:40:09  bigatti
// -- added "Gsl" prefix
//
// Revision 1.6  2014/04/17 13:37:07  bigatti
// -- MatrixViews --> MatrixView
//
// Revision 1.5  2013/07/23 12:32:13  bigatti
// -- removed indlude library.H
//
// Revision 1.4  2011/04/27 09:41:04  bigatti
// -- minor reshape
//

#endif
