//   Copyright (c) 2012  John Abbott,  Anna M Bigatti
//   Orig authors: 2012 Anna M. Bigatti, Bruno Simoes
//
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

#include "BuiltInFunctions.H"
#include "BuiltInOneLiners.H"

using namespace std;
using namespace boost;
using namespace boost::iostreams;
using namespace CoCoA::AST;
using namespace CoCoA::LexerNS;
using namespace CoCoA::ParserNS;

namespace CoCoA {
namespace InterpreterNS {

//  extern std::vector<NameFunPair> builtIns; // declared in BuiltInFunctions.C

#ifndef CoCoA_WITH_GSL

  DECLARE_MISSING_EXTLIB(GslLU, "GSL")
  DECLARE_MISSING_EXTLIB(GslCholeskyDecomposition, "GSL")
  DECLARE_MISSING_EXTLIB(GslSVD, "GSL")
  DECLARE_MISSING_EXTLIB(GslSingularValues, "GSL")
  DECLARE_MISSING_EXTLIB(GslQR, "GSL")
  DECLARE_MISSING_EXTLIB(GslQRPT, "GSL")
  DECLARE_MISSING_EXTLIB(GslSolverLU, "GSL")
  DECLARE_MISSING_EXTLIB(GslSolverSVD, "GSL")
  DECLARE_MISSING_EXTLIB(GslSolverQR, "GSL")
  DECLARE_MISSING_EXTLIB(GslBidiagDecomposition, "GSL")
  DECLARE_MISSING_EXTLIB(GslBidiagDecompUnpack, "GSL")

#else

  //---- one-liners
  DECLARE_COCOALIB_FUNCTION1(GslLU, MAT)
  DECLARE_COCOALIB_FUNCTION1(GslCholeskyDecomposition, MAT)
  DECLARE_COCOALIB_FUNCTION1(GslSVD, MAT)
  DECLARE_COCOALIB_FUNCTION1(GslSingularValues, MAT)
  DECLARE_COCOALIB_FUNCTION1(GslQR, MAT)
  DECLARE_COCOALIB_FUNCTION1(GslQRPT, MAT) 

  
  //---- other functions
  DECLARE_STD_BUILTIN_FUNCTION(GslSolverLU, 2) {
    intrusive_ptr<MAT> A = runtimeEnv->evalArgAs<MAT>(ARG(0));
    vector<RingElem> b = runtimeEnv->evalArgAsListOfRingElem(ARG(1));
    return Value::from(GslSolverLU(A->theMatrix, b));
  }
  END_STD_BUILTIN_FUNCTION


  DECLARE_STD_BUILTIN_FUNCTION(GslSolverSVD, 2) {
    intrusive_ptr<MAT> A = runtimeEnv->evalArgAs<MAT>(ARG(0));
    vector<RingElem> b = runtimeEnv->evalArgAsListOfRingElem(ARG(1));//, R->theRing);
    return Value::from(GslSolveSVD(A->theMatrix, b));
  }
  END_STD_BUILTIN_FUNCTION


  DECLARE_STD_BUILTIN_FUNCTION(GslSolverQR, 2) {
    intrusive_ptr<MAT> A = runtimeEnv->evalArgAs<MAT>(ARG(0));
    vector<RingElem> b = runtimeEnv->evalArgAsListOfRingElem(ARG(1));
    return Value::from(GslSolverQR(A->theMatrix, b));
  }
  END_STD_BUILTIN_FUNCTION


  DECLARE_STD_BUILTIN_FUNCTION(GslBidiagDecomposition, 4) {
    intrusive_ptr<MAT> A = runtimeEnv->evalArgAs<MAT>(ARG(0));
    intrusive_ptr<MAT> D = runtimeEnv->evalArgAs<MAT>(ARG(1));
    vector<RingElem> tau_U = runtimeEnv->evalArgAsListOf<RINGELEM>(ARG(2));
    vector<RingElem> tau_V = runtimeEnv->evalArgAsListOf<RINGELEM>(ARG(3)); 
    GslBidiagDecomposition(A->theMatrix, D->theMatrix, tau_U, tau_V);
    return VoidValue::theInstance;
  }
  END_STD_BUILTIN_FUNCTION


  DECLARE_STD_BUILTIN_FUNCTION(GslBidiagDecompUnpack, 5) {
    intrusive_ptr<MAT> A = runtimeEnv->evalArgAs<MAT>(ARG(0));
    intrusive_ptr<MAT> U = runtimeEnv->evalArgAs<MAT>(ARG(1));
    intrusive_ptr<MAT> V = runtimeEnv->evalArgAs<MAT>(ARG(2));
    vector<RingElem> tau_U = runtimeEnv->evalArgAsListOf<RINGELEM>(ARG(3));
    vector<RingElem> tau_V = runtimeEnv->evalArgAsListOf<RINGELEM>(ARG(4)); 
    GslBidiagDecompUnpack(A->theMatrix, U->theMatrix, V->theMatrix, tau_U, tau_V);
    return VoidValue::theInstance;
  }
  END_STD_BUILTIN_FUNCTION


#endif //  CoCoA_WITH_GSL

} // namespace AST
} // namespace CoCoA
