//   Copyright (c)  2012-2021  John Abbott and Anna M. Bigatti
//   Orig authors: 2012-2017 Anna M. Bigatti, Christof Soeger
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

#ifndef CoCoA_WITH_NORMALIZ

  DECLARE_MISSING_EXTLIB(NmzSetVerbosityLevel, "NORMALIZ")
  DECLARE_MISSING_EXTLIB(NmzVerbosityLevel, "NORMALIZ")
  DECLARE_MISSING_EXTLIB(NmzNormalToricRing, "NORMALIZ")
  DECLARE_MISSING_EXTLIB(NmzHilbertBasis,    "NORMALIZ")
  DECLARE_MISSING_EXTLIB(NmzHilbertBasisKer, "NORMALIZ")
  DECLARE_MISSING_EXTLIB(NmzIntClosureToricRing, "NORMALIZ")
  DECLARE_MISSING_EXTLIB(NmzIntClosureMonIdeal, "NORMALIZ")
  DECLARE_MISSING_EXTLIB(NmzEhrhartRing, "NORMALIZ")
  DECLARE_MISSING_EXTLIB(NmzTorusInvariants, "NORMALIZ")
  DECLARE_MISSING_EXTLIB(NmzFiniteDiagInvariants, "NORMALIZ")
  DECLARE_MISSING_EXTLIB(NmzDiagInvariants, "NORMALIZ")
  DECLARE_MISSING_EXTLIB(NmzIntersectionValRings, "NORMALIZ")
  DECLARE_MISSING_EXTLIB(NmzComputation, "NORMALIZ")

#else

// ---- Nmz Supplement ----

  long NmzVerbosityLevel_forC5()
  { return Normaliz::VerbosityLevel(); }

  matrix NmzHilbertBasis_forC5(ConstMatrixView M)
  {
    return NewDenseMat(RingOf(M), Normaliz::HilbertBasis(Normaliz::cone(libnormaliz::Type::normalization, M)));
  }

  matrix NmzHilbertBasisKer_forC5(ConstMatrixView M)
  {
    return NewDenseMat(RingOf(M), Normaliz::HilbertBasis(Normaliz::cone(libnormaliz::Type::equations, M)));
  }

  std::vector<RingElem> NmzNormalToricRing_forC5(const std::vector<RingElem>& v)
  {
    //convert it to a PPVector, run Normaliz and convert back to vector<RingElem>
    std::vector<RingElem> res;
    convert(res,owner(v[0]),Normaliz::NormalToricRing(Normaliz::MonomialsToPPV(v)));
    return res;
  }

  std::vector<RingElem> NmzIntClosureToricRing_forC5(const std::vector<RingElem>& v)
  {
    //convert it to a PPVector
    PPVector ppv = Normaliz::MonomialsToPPV(v);
    //run Normaliz and convert back to vector<RingElem>
    std::vector<RingElem> res;
    convert(res,owner(v[0]),Normaliz::IntClosureToricRing(ppv));
    return res;
  }


  std::vector<RingElem> NmzIntClosureMonIdeal_forC5(const std::vector<RingElem>& v)
  {
    //convert it to a PPVector
    PPVector ppv = Normaliz::MonomialsToPPV(v);
    //run Normaliz and convert back to vector<RingElem>
    std::vector<RingElem> res;
    convert(res,owner(v[0]),Normaliz::IntClosureMonIdeal(ppv));
    return res;
  }


// ---- Nmz OneLiners ----

DECLARE_COCOALIBFORC5_FUNCTION1(NmzHilbertBasis, MAT)
DECLARE_COCOALIBFORC5_FUNCTION1(NmzHilbertBasisKer, MAT)
DECLARE_COCOALIBFORC5_FUNCTION0(NmzVerbosityLevel)

// ---- Nmz STD_BUILTIN_FUNCTION ----

DECLARE_STD_BUILTIN_FUNCTION(NmzSetVerbosityLevel, 1)
{
  const long n = runtimeEnv->evalArgAsLong(ARG(0));
  Normaliz::SetVerbosityLevel(n);
  return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(NmzNormalToricRing, 1) {
  vector<RingElem> v = runtimeEnv->evalArgAsListOfRingElem(ARG(0));
  return Value::from(NmzNormalToricRing_forC5(v));
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(NmzIntClosureToricRing, 1) {
  vector<RingElem> v = runtimeEnv->evalArgAsListOfRingElem(ARG(0));
  return Value::from(NmzIntClosureToricRing_forC5(v));
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(NmzEhrhartRing, 2) {
  vector<RingElem> v = runtimeEnv->evalArgAsListOfRingElem(ARG(0));
  RingElem sym = runtimeEnv->evalArgAs<RINGELEM>(ARG(1))->theRingElem;
  long sym_pos;
  if (!IsIndet(sym_pos, sym))
     throw RuntimeException("not an indet", ARG(1).exp);
  std::vector<RingElem> v_res;
  convert(v_res, owner(sym), Normaliz::EhrhartRing(Normaliz::MonomialsToPPV(v),sym_pos));
  return Value::from(v_res);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(NmzTorusInvariants, 2) {
  vector<vector<BigInt> > Mat = Normaliz::MatrixToVecVecBigInt((runtimeEnv->evalArgAs<MAT>(ARG(0)))->theMatrix);
  SparsePolyRing R = (runtimeEnv->evalArgAs<RING>(ARG(1)))->theRing;
  std::vector<RingElem> v;
  convert(v, R, Normaliz::TorusInvariants(Mat,PPM(R)));
  return Value::from(v);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(NmzFiniteDiagInvariants, 2) {
  vector<vector<BigInt> > Mat = Normaliz::MatrixToVecVecBigInt((runtimeEnv->evalArgAs<MAT>(ARG(0)))->theMatrix);
  SparsePolyRing R = (runtimeEnv->evalArgAs<RING>(ARG(1)))->theRing;
  std::vector<RingElem> v;
  convert(v, R, Normaliz::FiniteDiagInvariants(Mat,PPM(R)));
  return Value::from(v);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(NmzDiagInvariants, 3) {
  vector<vector<BigInt> > MatT = Normaliz::MatrixToVecVecBigInt((runtimeEnv->evalArgAs<MAT>(ARG(0)))->theMatrix);
  vector<vector<BigInt> > MatD = Normaliz::MatrixToVecVecBigInt((runtimeEnv->evalArgAs<MAT>(ARG(1)))->theMatrix);
  SparsePolyRing R = (runtimeEnv->evalArgAs<RING>(ARG(2)))->theRing;
  std::vector<RingElem> v;
  convert(v, R, Normaliz::DiagInvariants(MatT,MatD,PPM(R)));
  return Value::from(v);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(NmzIntersectionValRings, 2) {
  vector<vector<BigInt> > Mat = Normaliz::MatrixToVecVecBigInt((runtimeEnv->evalArgAs<MAT>(ARG(0)))->theMatrix);
  SparsePolyRing R = (runtimeEnv->evalArgAs<RING>(ARG(1)))->theRing;
  std::vector<RingElem> v;
  convert(v, R, Normaliz::IntersectionValRings(Mat,PPM(R)));
  return Value::from(v);
}
END_STD_BUILTIN_FUNCTION

// ---- Nmz BUILTIN_FUNCTION (+ variable number of args) ----

DECLARE_ARITYCHECK_FUNCTION(NmzIntClosureMonIdeal) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(NmzIntClosureMonIdeal) {
  invocationExpression->checkNumberOfArgs(1,2);// variable number of args
  vector<RingElem> v = runtimeEnv->evalArgAsListOfRingElem(ARG(0));
  std::vector<RingElem> v_res;
  if (invocationExpression->args.size()==1)
  {
    v_res = NmzIntClosureMonIdeal_forC5(v);
  }
  else
  {
    RingElem sym = runtimeEnv->evalArgAs<RINGELEM>(ARG(1))->theRingElem;
    long sym_pos;
    if (!IsIndet(sym_pos, sym))
       throw RuntimeException("not an indet", ARG(1).exp);
    convert(v_res, owner(sym), Normaliz::IntClosureMonIdeal(Normaliz::MonomialsToPPV(v),sym_pos));
  }
  return Value::from(v_res);
}

DECLARE_ARITYCHECK_FUNCTION(NmzComputation) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(NmzComputation) {
  invocationExpression->checkNumberOfArgs(1,2);// variable number of args
  intrusive_ptr<RECORD> cone = runtimeEnv->evalArgAs<RECORD>(ARG(0));
  vector<string> fieldnames = cone->myFieldNamesStrings();

  intrusive_ptr<RECORD> out(new RECORD);
  map< libnormaliz::InputType, vector<vector<BigInt> > > InputMap;
  try {
  for (long i=0; i<len(fieldnames); ++i)
  {
    intrusive_ptr<RightValue> x = boost::dynamic_pointer_cast<RightValue>(cone->getField(fieldnames[i], ARG(0).exp));
    if (intrusive_ptr<MAT> M = dynamic_pointer_cast<MAT>(x)) {
      // transform M to a vector< vector<BigInt> >
      // the following line may throw a libnormaliz::BadInputException
      InputMap[Normaliz::ToInputType(fieldnames[i])] = Normaliz::MatrixToVecVecBigInt(M->theMatrix);
    } //TODO allow also [[..],..,[..]] input?
    else
      throw WrongTypeException(MAT::type->name
//                             RINGELEM::type->name + " or " +
                               + " as field \"" + fieldnames[i] + "\" of the record",
                               x->getType()->name, ARG(0).exp);
  }

  Normaliz::cone C(InputMap);

  // collect what to compute
  libnormaliz::ConeProperties CPs;
  if (invocationExpression->args.size()==1)
  {
    CPs.set(libnormaliz::ConeProperty::DefaultMode);
  }
  else
  {
    vector<string> vs = runtimeEnv->evalArgAsListOf<STRING>(ARG(1));
    for (long i=0; i<len(vs); ++i)
    {
      // invalid strings will cause a libnormaliz::BadInputException
      CPs.set(libnormaliz::toConeProperty(vs[i]));
    }
  }

  C.myComputation(CPs);

  // fill the return record with everthing we can
  for (size_t i=0; i<libnormaliz::ConeProperty::EnumSize; i++) {
    libnormaliz::ConeProperty::Enum prop_i = static_cast<libnormaliz::ConeProperty::Enum>(i);
    if (C.isComputed(prop_i)) switch (prop_i) {
      case libnormaliz::ConeProperty::HilbertBasis:
        out->setField(libnormaliz::toString(prop_i), Value::from(HilbertBasis(C)), ARG(0).exp);
        break;
      case libnormaliz::ConeProperty::ModuleGenerators:
        out->setField(libnormaliz::toString(prop_i), Value::from(ModuleGenerators(C)), ARG(0).exp);
        break;
        // case libnormaliz::ConeProperty::TriangulationGenerators: // v 3.8.9
      case libnormaliz::ConeProperty::Triangulation: // v 3.8.10
        out->setField(libnormaliz::toString(prop_i), Value::from(TriangulationGenerators(C)), ARG(0).exp);
        break;
      case libnormaliz::ConeProperty::ExtremeRays:
        out->setField(libnormaliz::toString(prop_i), Value::from(ExtremeRays(C)), ARG(0).exp);
        out->setField("Rank", Value::from(Normaliz::rank(C)), ARG(0).exp);
        break;
      case libnormaliz::ConeProperty::VerticesOfPolyhedron:
        out->setField(libnormaliz::toString(prop_i), Value::from(VerticesOfPolyhedron(C)), ARG(0).exp);
        break;
      case libnormaliz::ConeProperty::Deg1Elements:
        out->setField(libnormaliz::toString(prop_i), Value::from(Deg1Elements(C)), ARG(0).exp);
        break;
      case libnormaliz::ConeProperty::OriginalMonoidGenerators:
        out->setField(libnormaliz::toString(prop_i), Value::from(GeneratorsOfToricRing(C)), ARG(0).exp);
        break;
      case libnormaliz::ConeProperty::SupportHyperplanes:
        out->setField(libnormaliz::toString(prop_i), Value::from(SupportHyperplanes(C)), ARG(0).exp);
        out->setField("Equations", Value::from(Equations(C)), ARG(0).exp);
        out->setField("Congruences", Value::from(Congruences(C)), ARG(0).exp);
        break;
      case libnormaliz::ConeProperty::ExcludedFaces:
        out->setField(libnormaliz::toString(prop_i), Value::from(ExcludedFaces(C)), ARG(0).exp);
        break;
      case libnormaliz::ConeProperty::HilbertSeries:
        out->setField(libnormaliz::toString(prop_i), Value::from(HilbertSeries(C)), ARG(0).exp);
        out->setField("HilbertQuasiPolynomial", Value::from(HilbertQuasiPoly(C)), ARG(0).exp);
        break;
      case libnormaliz::ConeProperty::Multiplicity:
        out->setField(libnormaliz::toString(prop_i), Value::from(multiplicity(C)), ARG(0).exp);
        break;
      case libnormaliz::ConeProperty::Grading:
        out->setField(libnormaliz::toString(prop_i), Value::from(grading(C)), ARG(0).exp);
        break;
      case libnormaliz::ConeProperty::IsDeg1HilbertBasis:
        out->setField(libnormaliz::toString(prop_i), Value::from(IsDeg1HilbertBasis(C)), ARG(0).exp);
        break;
      case libnormaliz::ConeProperty::IsPointed:
        out->setField(libnormaliz::toString(prop_i), Value::from(IsPointed(C)), ARG(0).exp);
        break;
      case libnormaliz::ConeProperty::IsIntegrallyClosed:
        out->setField(libnormaliz::toString(prop_i), Value::from(IsIntegrallyClosed(C)), ARG(0).exp);
        break;
      case libnormaliz::ConeProperty::RecessionRank:
        out->setField(libnormaliz::toString(prop_i), Value::from(RecessionRank(C)), ARG(0).exp);
        break;
      case libnormaliz::ConeProperty::AffineDim:
        out->setField(libnormaliz::toString(prop_i), Value::from(AffineDim(C)), ARG(0).exp);
        break;
      case libnormaliz::ConeProperty::ModuleRank:
        out->setField(libnormaliz::toString(prop_i), Value::from(ModuleRank(C)), ARG(0).exp);
        break;
      case libnormaliz::ConeProperty::Dehomogenization:
        out->setField(libnormaliz::toString(prop_i), Value::from(dehomogenization(C)), ARG(0).exp);
        break;
      default:
        break;
    }
  }

  out->setField("EmbeddingDim", Value::from(EmbeddingDim(C)), ARG(0).exp);
  out->setField("IsInhomogeneous", Value::from(IsInhomogeneous(C)), ARG(0).exp);

  } catch (const libnormaliz::NormalizException& e) {
    CoCoA_THROW_ERROR(e.what(),"NmzComputation");
  }

  return out;
}

#endif // CoCoA_WITH_NORMALIZ

} // namespace AST
} // namespace CoCoA
