//   Copyright (c) 2013  Anna M. Bigatti,  John Abbott
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

#ifndef CoCoA_WITH_FROBBY

  DECLARE_MISSING_EXTLIB(FrbDimension, "FROBBY")
  DECLARE_MISSING_EXTLIB(FrbMaximalStandardMonomials, "FROBBY")
  DECLARE_MISSING_EXTLIB(FrbMultigradedHilbertPoincareNumerator, "FROBBY")
  DECLARE_MISSING_EXTLIB(FrbTotalDegreeHilbertPoincareNumerator, "FROBBY")
  DECLARE_MISSING_EXTLIB(FrbAssociatedPrimes, "FROBBY")
  DECLARE_MISSING_EXTLIB(FrbIrreducibleDecomposition, "FROBBY")
  DECLARE_MISSING_EXTLIB(FrbPrimaryDecomposition, "FROBBY")
  DECLARE_MISSING_EXTLIB(FrbAlexanderDual, "FROBBY")

#else

  //----- supplement
  std::vector<ideal> FrbPrimaryDecomposition_forC5(const ideal& I)
  {
    std::vector<ideal> v;
    FrbPrimaryDecomposition(v,I);
    return v;
  }

  std::vector<ideal> FrbIrreducibleDecomposition_forC5(const ideal& I)
  {
    std::vector<ideal> v;
    FrbIrreducibleDecomposition(v,I);
    return v;
  }

  std::vector<ideal> FrbAssociatedPrimes_forC5(const ideal& I)
  {
    std::vector<ideal> v;
    FrbAssociatedPrimes(v,I);
    return v;
  }

  ideal FrbAlexanderDual_forC5(const ideal& I, ConstRefRingElem t)
  {
    const char* const FnName = "FrbAlexanderDual";
    const ring& R = owner(t);
    if (!IsPolyRing(R))  CoCoA_THROW_ERROR2(ERR::ReqElemPolyRing, FnName);
    if (RingOf(I) != R)  CoCoA_THROW_ERROR2(ERR::MixedRings, FnName);
    if (!IsMonomial(t) || !IsOne(LC(t)))
      CoCoA_THROW_ERROR1(string("Reqire power-product (monic monomial) in ") + FnName);
    return FrbAlexanderDual(I, LPP(t));
  }


//----- one-liners
DECLARE_COCOALIB_FUNCTION1(FrbDimension, IDEAL)
DECLARE_COCOALIB_FUNCTION1(FrbMaximalStandardMonomials, IDEAL)
DECLARE_COCOALIB_FUNCTION1(FrbMultigradedHilbertPoincareNumerator, IDEAL)
DECLARE_COCOALIB_FUNCTION1(FrbTotalDegreeHilbertPoincareNumerator, IDEAL)

DECLARE_COCOALIBFORC5_FUNCTION1(FrbAssociatedPrimes, IDEAL)
DECLARE_COCOALIBFORC5_FUNCTION1(FrbIrreducibleDecomposition, IDEAL)
DECLARE_COCOALIBFORC5_FUNCTION1(FrbPrimaryDecomposition, IDEAL)

//----- other
// variable number of args
DECLARE_ARITYCHECK_FUNCTION(FrbAlexanderDual) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(FrbAlexanderDual) {
  invocationExpression->checkNumberOfArgs(1,2);
	intrusive_ptr<IDEAL> I = runtimeEnv->evalArgAs<IDEAL>(ARG(0));
  if (invocationExpression->args.size()==1)
    return Value::from(FrbAlexanderDual(I->theIdeal));
	intrusive_ptr<RINGELEM> t = runtimeEnv->evalArgAs<RINGELEM>(ARG(1));
  return Value::from(FrbAlexanderDual_forC5(I->theIdeal, t->theRingElem));
}


#endif // CoCoA_WITH_FROBBY

} // namespace AST
} // namespace CoCoA
