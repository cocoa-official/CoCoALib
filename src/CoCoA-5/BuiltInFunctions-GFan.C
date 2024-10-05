//   Copyright (c) 2015  Anna M. Bigatti,  John Abbott
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

#ifndef CoCoA_WITH_GFAN

  DECLARE_MISSING_EXTLIB(GFanRelativeInteriorPoint, "GFAN")
  DECLARE_MISSING_EXTLIB(GFanGeneratorsOfSpan, "GFAN")
  DECLARE_MISSING_EXTLIB(GFanGeneratorsOfLinealitySpace, "GFAN")
  DECLARE_MISSING_EXTLIB(GFanGetFacets, "GFAN")
  DECLARE_MISSING_EXTLIB(GFanGetImpliedEquations, "GFAN")
  DECLARE_MISSING_EXTLIB(GFanGetUniquePoint, "GFAN")
  DECLARE_MISSING_EXTLIB(GFanGetAmbientDimension, "GFAN")
  DECLARE_MISSING_EXTLIB(GFanGetDimension, "GFAN")
  DECLARE_MISSING_EXTLIB(GFanGetCodimension, "GFAN")
  DECLARE_MISSING_EXTLIB(GFanGetDimensionOfLinealitySpace, "GFAN")
  DECLARE_MISSING_EXTLIB(GFanContainsPositiveVector, "GFAN")

#else

//----- CoCoALibSupplement

  matrix GFanRelativeInteriorPoint_forC5(ConstMatrixView EqMat, ConstMatrixView IneqMat)
  { return GFan::RelativeInteriorPoint(GFan::cone(EqMat, IneqMat)); }

  matrix GFanGeneratorsOfSpan_forC5(ConstMatrixView EqMat, ConstMatrixView IneqMat)
  { return GFan::GeneratorsOfSpan(GFan::cone(EqMat, IneqMat)); }

  matrix GFanGeneratorsOfLinealitySpace_forC5(ConstMatrixView EqMat, ConstMatrixView IneqMat)
  { return GFan::GeneratorsOfLinealitySpace(GFan::cone(EqMat, IneqMat)); }

  matrix GFanGetFacets_forC5(ConstMatrixView EqMat, ConstMatrixView IneqMat)
  { return GFan::GetFacets(GFan::cone(EqMat, IneqMat)); }

  matrix GFanGetImpliedEquations_forC5(ConstMatrixView EqMat, ConstMatrixView IneqMat)
  { return GFan::GetFacets(GFan::cone(EqMat, IneqMat)); }

  matrix GFanGetUniquePoint_forC5(ConstMatrixView EqMat, ConstMatrixView IneqMat)
  { return GFan::GetUniquePoint(GFan::cone(EqMat, IneqMat)); }

  long GFanGetAmbientDimension_forC5(ConstMatrixView EqMat, ConstMatrixView IneqMat)
  { return GFan::GetAmbientDimension(GFan::cone(EqMat, IneqMat)); }

  long GFanGetDimension_forC5(ConstMatrixView EqMat, ConstMatrixView IneqMat)
  { return GFan::GetDimension(GFan::cone(EqMat, IneqMat)); }

  long GFanGetCodimension_forC5(ConstMatrixView EqMat, ConstMatrixView IneqMat)
  { return GFan::GetCodimension(GFan::cone(EqMat, IneqMat)); }

  long GFanGetDimensionOfLinealitySpace_forC5(ConstMatrixView EqMat, ConstMatrixView IneqMat)
  { return GFan::GetDimensionOfLinealitySpace(GFan::cone(EqMat, IneqMat)); }

  long GFanContainsPositiveVector_forC5(ConstMatrixView EqMat, ConstMatrixView IneqMat)
  { return GFan::ContainsPositiveVector(GFan::cone(EqMat, IneqMat)); }

//----- one-liners
  DECLARE_COCOALIBFORC5_FUNCTION2(GFanRelativeInteriorPoint, MAT, MAT)
  DECLARE_COCOALIBFORC5_FUNCTION2(GFanGeneratorsOfSpan, MAT, MAT)
  DECLARE_COCOALIBFORC5_FUNCTION2(GFanGeneratorsOfLinealitySpace, MAT, MAT)
  DECLARE_COCOALIBFORC5_FUNCTION2(GFanGetFacets, MAT, MAT)
  DECLARE_COCOALIBFORC5_FUNCTION2(GFanGetImpliedEquations, MAT, MAT)
  DECLARE_COCOALIBFORC5_FUNCTION2(GFanGetUniquePoint, MAT, MAT)
  DECLARE_COCOALIBFORC5_FUNCTION2(GFanGetAmbientDimension, MAT, MAT)
  DECLARE_COCOALIBFORC5_FUNCTION2(GFanGetDimension, MAT, MAT)
  DECLARE_COCOALIBFORC5_FUNCTION2(GFanGetCodimension, MAT, MAT)
  DECLARE_COCOALIBFORC5_FUNCTION2(GFanGetDimensionOfLinealitySpace, MAT, MAT)
  DECLARE_COCOALIBFORC5_FUNCTION2(GFanContainsPositiveVector, MAT, MAT)

#endif // CoCoA_WITH_GFAN

} // namespace AST
} // namespace CoCoA
