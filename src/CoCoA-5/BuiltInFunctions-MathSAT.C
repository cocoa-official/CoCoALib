//   Copyright (c) 2017  John Abbott and Anna M. Bigatti
//   Orig authors: 2017 Anna M. Bigatti
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

#ifndef CoCoA_WITH_MATHSAT

  DECLARE_MISSING_EXTLIB(MSatIsSatisfiable, "MathSAT")
  DECLARE_MISSING_EXTLIB(MSatLinSolve, "MathSAT")
  DECLARE_MISSING_EXTLIB(MSatSolution, "MathSAT")

#else


// ---- MSat STD_BUILTIN_FUNCTION ----

DECLARE_STD_BUILTIN_FUNCTION(MSatLinSolve, 1) {
  intrusive_ptr<RECORD> r = runtimeEnv->evalArgAs<RECORD>(ARG(0));
  vector<string> FieldNames = r->myFieldNamesStrings();

  // create a MathSAT environment
  MathSAT::env E;
  
  for (long i=0; i<len(FieldNames); ++i)
  {
    intrusive_ptr<RightValue> x = boost::dynamic_pointer_cast<RightValue>(r->getField(FieldNames[i], ARG(0).exp));
    if (intrusive_ptr<MAT> M = dynamic_pointer_cast<MAT>(x))
    {
      if (!IsZero(characteristic(RingOf(M->theMatrix))))
        throw RuntimeException("Matrix must have rational entries", ARG(0).exp);
      MathSAT::AddConstraint(E, MathSAT::ToRelOp(FieldNames[i]), M->theMatrix);
    }
    else
      throw WrongTypeException(MAT::type->name
                               // + " or " + LIST::type->name
                               + " as field \"" + FieldNames[i]
                               + "\" of the record",
                               x->getType()->name, ARG(0).exp);
  }
  return Value::from(MathSAT::LinSolve(E));
}
END_STD_BUILTIN_FUNCTION


#endif // CoCoA_WITH_MATHSAT

} // namespace AST
} // namespace CoCoA
