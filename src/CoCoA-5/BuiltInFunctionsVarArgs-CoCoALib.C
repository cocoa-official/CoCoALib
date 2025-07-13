//   Copyright (c) 2010-2017,2022   John Abbott, Anna M. Bigatti
//   Authors: 2010-2011 Giovanni Lagorio
//   Authors: 2011--  John Abbott, Anna M. Bigatti
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
#include "CoCoALibSupplement.H"

using namespace std;
using namespace boost;
using namespace boost::iostreams;
using namespace CoCoA::AST;
using namespace CoCoA::LexerNS;
using namespace CoCoA::ParserNS;

namespace CoCoA {
namespace InterpreterNS {

//  extern std::vector<NameFunPair> builtIns; // declared in BuiltInFunctions.C
//----------------------------------------------------------------------

namespace  // anonymous  =============== IsVectorBigInt =================
{
  bool IsVectorBigInt(std::vector<BigInt>& BigIntVec, const intrusive_ptr<LIST> l)
  {
    vector<BigInt> v;
    LIST::ContainerType::size_type size = l->size();
    for (unsigned long i=0; i<size; ++i)
      if (const boost::intrusive_ptr<INT> n = boost::dynamic_pointer_cast<INT>(l->getValue(i)))
        v.push_back(n->theBigInt);
      else
        return false;
    swap(v, BigIntVec);
    return true;
  }  
} // anonymous namespace



//------------------------------------------------------------------

DECLARE_ARITYCHECK_FUNCTION(RandomSubsetIndices) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(RandomSubsetIndices) {  // JAA
  invocationExpression->checkNumberOfArgs(1,2);// variable number of args
  const long n = runtimeEnv->evalArgAsLong(ARG(0));
  vector<long> ans;
  if (invocationExpression->args.size()==1)
  {
    ans = RandomSubsetIndices(n);
  }
  else
  {
    const long r = runtimeEnv->evalArgAsLong(ARG(1));
    ans = RandomSubsetIndices(n,r);
  }
  const long card = len(ans);
  for (long i=0; i < card; ++i)
    ++ans[i];
  return Value::from(ans);
}


DECLARE_ARITYCHECK_FUNCTION(deg) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(deg) {  // AMB
  invocationExpression->checkNumberOfArgs(1,2);// variable number of args
  intrusive_ptr<RINGELEM> a = runtimeEnv->evalArgAs<RINGELEM>(ARG(0));
  if (invocationExpression->args.size()==1)
    return Value::from(deg(a->theRingElem));
  intrusive_ptr<RINGELEM> b = runtimeEnv->evalArgAs<RINGELEM>(ARG(1));
  long i;
  if (!IsIndet(i,b->theRingElem))
    throw RuntimeException("must be an indet", ARG(1).exp);
  return Value::from(deg(a->theRingElem, i));
}


DECLARE_ARITYCHECK_FUNCTION(coefficients) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(coefficients) {  // AMB
  invocationExpression->checkNumberOfArgs(1,2);// variable number of args
  intrusive_ptr<RINGELEM> a = runtimeEnv->evalArgAs<RINGELEM>(ARG(0));
  if (invocationExpression->args.size()==1)
    return Value::from(coefficients_forC5(a->theRingElem));
  vector<RingElem> v = runtimeEnv->evalArgAsListOfRingElem(ARG(1));
  vector<RingElem> res;
  for (long i=0; i<len(v); ++i)
    res.push_back(CoeffOfTerm_forC5(a->theRingElem, v[i]));
  return Value::from(res);
}


// DECLARE_ARITYCHECK_FUNCTION(NewMat) { return (2<=nArg) && (nArg<=3); }
// DECLARE_BUILTIN_FUNCTION(NewMat) {
//  invocationExpression->checkNumberOfArgs(2,3);
//   if (invocationExpression->args.size()==2)
//     throw RuntimeException("NewMat(NR,NC) not allowed, use NewMat(R:RING,NR:INT,NC:INT) instead", ARG(0).exp);
//   intrusive_ptr<RING> R = runtimeEnv->evalArgAs<RING>(ARG(0));
//   intrusive_ptr<INT> nr = runtimeEnv->evalArgAs<INT>(ARG(1));
//   intrusive_ptr<INT> nc = runtimeEnv->evalArgAs<INT>(ARG(2));
//   return Value::from(ZeroMat(R->theRing, ConvertTo<long>(nr->theBigInt), ConvertTo<long>(nc->theBigInt)));
// }

// DECLARE_ARITYCHECK_FUNCTION(GrammSchmidtRows) { return (1<=nArg) && (nArg<=2); }
// DECLARE_BUILTIN_FUNCTION(GrammSchmidtRows) {  // AMB
//  invocationExpression->checkNumberOfArgs(1,2);
//  intrusive_ptr<MAT> M = runtimeEnv->obtainUnshared<MAT>(ARG(0));
//   if (invocationExpression->args.size()==1)
//   {
//     GrammSchmidtRows(M->theMatrix);
//     return VoidValue::theInstance;
//   }
//  intrusive_ptr<INT> N = runtimeEnv->evalArgAs<INT>(ARG(1));
//   long n;
//   if (!IsConvertible(n, N->theBigInt))
//     throw RuntimeException("invalid row index", ARG(1).exp);
//   GrammSchmidtRows(M->theMatrix, n);
//   return VoidValue::theInstance;
// }


DECLARE_ARITYCHECK_FUNCTION(RatReconstructByContFrac) { return (2<=nArg) && (nArg<=3); }
DECLARE_BUILTIN_FUNCTION(RatReconstructByContFrac) {
  invocationExpression->checkNumberOfArgs(2,3);
  intrusive_ptr<INT> X = runtimeEnv->evalArgAs<INT>(ARG(0));
  intrusive_ptr<INT> M = runtimeEnv->evalArgAs<INT>(ARG(1));  
  //        BigInt threshold; // Determine threshold: 0 means use default value
  long LogEps = 20; // default value
  if (invocationExpression->args.size()==3)
  {
    LogEps = runtimeEnv->evalArgAsLong(ARG(2));
    //          intrusive_ptr<INT> thresh = runtimeEnv->evalArgAs<INT>(ARG(2));
    //          if (thresh->theBigInt < 0) throw RuntimeException("Threshold must be >= 0", ARG(2).exp);
    //          threshold = thresh->theBigInt;
  }
  RatReconstructByContFrac reconstructor(LogEps);
  reconstructor.myAddInfo(X->theBigInt, M->theBigInt);
  // Create return value: ... record
  intrusive_ptr<RECORD> ans(new RECORD);
  ans->setFieldNoCheck("failed", Value::from(!IsConvincing(reconstructor)));
  if (IsConvincing(reconstructor))
    ans->setFieldNoCheck("ReconstructedRat", Value::from(ReconstructedRat(reconstructor)));
  return ans;
}

 
DECLARE_ARITYCHECK_FUNCTION(RatReconstructByLattice) { return (2<=nArg) && (nArg<=3); }
DECLARE_BUILTIN_FUNCTION(RatReconstructByLattice) {
  invocationExpression->checkNumberOfArgs(2,3);
  intrusive_ptr<INT> X = runtimeEnv->evalArgAs<INT>(ARG(0));
  intrusive_ptr<INT> M = runtimeEnv->evalArgAs<INT>(ARG(1));
  BigInt SafetyFactor; // Determine threshold: 0 means use default value
  if (invocationExpression->args.size()==3)
  {
    intrusive_ptr<INT> safety = runtimeEnv->evalArgAs<INT>(ARG(2));
    //if (safety->theBigInt < 0) throw RuntimeException("SafetyFactor must be >= 0", ARG(2).exp);
    SafetyFactor = safety->theBigInt;
  }
  RatReconstructByLattice reconstructor(SafetyFactor);
  reconstructor.myAddInfo(X->theBigInt, M->theBigInt);
  // Create return value: ... record
  intrusive_ptr<RECORD> ans(new RECORD);
  ans->setFieldNoCheck("failed", Value::from(!IsConvincing(reconstructor)));
  if (IsConvincing(reconstructor))
    ans->setFieldNoCheck("ReconstructedRat", Value::from(ReconstructedRat(reconstructor)));
  return ans;
}


DECLARE_ARITYCHECK_FUNCTION(MantissaAndExponent2) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(MantissaAndExponent2)
{
  invocationExpression->checkNumberOfArgs(1,2);// variable number of args
  MantExp2 ME;
  if (invocationExpression->args.size()==1)
  {
    intrusive_ptr<RightValue> x = runtimeEnv->evalArgAs<RINGELEM>(ARG(0));
    ME = MantissaAndExponent2(RefTo<RingElem>(x));
  }
  else
  {
  const long n = runtimeEnv->evalArgAsLong(ARG(1));
  if (n < 1) throw RuntimeException("Precision must be positive (and not too large)", ARG(1).exp);

  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<INT, RAT>(ARG(0), which);
  switch (which)
  {
  case 1: ME = MantissaAndExponent2(RefTo<BigInt>(x), n); break;
  case 2: ME = MantissaAndExponent2(RefTo<BigRat>(x), n); break;
  default: throw RuntimeException(ERRORMissingCode(x),invocationExpression);
  }
  }
  // Create return value: ... record
  intrusive_ptr<RECORD> ans(new RECORD);
  ans->setFieldNoCheck("mantissa", Value::from(ME.mySign * ME.myMantissa));
  ans->setFieldNoCheck("exponent", Value::from(ME.myExponent));
  ans->setFieldNoCheck("NumDigits", Value::from(ME.myNumDigits));
  return ans;
}


DECLARE_ARITYCHECK_FUNCTION(FloatStr) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(FloatStr)
{
  invocationExpression->checkNumberOfArgs(1,2);
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2orT3<INT, RAT, RINGELEM>(ARG(0), which);

  BigRat q;
  switch (which)
  {
  case 1: q = RefTo<BigInt>(x); break;
  case 2: q = RefTo<BigRat>(x); break;
  case 3: q = ConvertTo<BigRat>(RefTo<RingElem>(x)); break;
  default: return Value::from(false); // just to keep the compiler quiet
  }
  if (invocationExpression->args.size()==1)
    return Value::from(FloatStr(q));
  long n = runtimeEnv->evalArgAsLong(ARG(1));
  //  intrusive_ptr<INT> Ndigits = runtimeEnv->evalArgAs<INT>(ARG(1));
  if (n < 1) throw RuntimeException("nonsensical value", ARG(1).exp);
  return Value::from(FloatStr(q, n));
//???  return Value::from(ScientificStr(q, ConvertTo<long>(Ndigits, "nonsensical value")));
}


DECLARE_ARITYCHECK_FUNCTION(ScientificStr) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(ScientificStr)
{
  invocationExpression->checkNumberOfArgs(1,2);
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2orT3<INT, RAT, RINGELEM>(ARG(0), which);

  BigRat q;
  switch (which)
  {
  case 1: q = RefTo<BigInt>(x); break;
  case 2: q = RefTo<BigRat>(x); break;
  case 3: q = ConvertTo<BigRat>(RefTo<RingElem>(x)); break;
  default: return Value::from(false); // just to keep the compiler quiet
  }

  if (invocationExpression->args.size()==1)
    return Value::from(ScientificStr(q));
  long n = runtimeEnv->evalArgAsLong(ARG(1));
  //  intrusive_ptr<INT> Ndigits = runtimeEnv->evalArgAs<INT>(ARG(1));
  if (n < 1) throw RuntimeException("nonsensical value", ARG(1).exp);
  return Value::from(ScientificStr(q, n));
//???  return Value::from(ScientificStr(q, ConvertTo<long>(Ndigits, "nonsensical value")));
}


DECLARE_ARITYCHECK_FUNCTION(DecimalStr) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(DecimalStr)
{
  invocationExpression->checkNumberOfArgs(1,2);
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2orT3<INT, RAT, RINGELEM>(ARG(0), which);

  BigRat q;
  switch (which)
  {
  case 1: q = RefTo<BigInt>(x); break;
  case 2: q = RefTo<BigRat>(x); break;
  case 3: q = ConvertTo<BigRat>(RefTo<RingElem>(x)); break;
  }

  if (invocationExpression->args.size()==1)
    return Value::from(DecimalStr(q));
  long n = runtimeEnv->evalArgAsLong(ARG(1));
  //  intrusive_ptr<INT> Ndigits = runtimeEnv->evalArgAs<INT>(ARG(1));
  if (n < 1) throw RuntimeException("nonsensical value", ARG(1).exp);
  return Value::from(DecimalStr(q, n));
//???  return Value::from(DecimalStr(q, ConvertTo<long>(Ndigits, "nonsensical value")));
}


// still valid?  I do not think it does anything useful AMB 2015-09
// // variable number of args
// DECLARE_ARITYCHECK_FUNCTION(janet) { return 1<=nArg; }
// DECLARE_BUILTIN_FUNCTION(janet) {
//  const int nArgs = invocationExpression->args.size();
//  if (nArgs==0)
//    throw RuntimeException("Wrong number of arguments; found: 0, expecting: at least 1", invocationExpression);
//  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAs<RightValue>(ARG(0));
//   if (invocationExpression->args.size()==2)
//     if (intrusive_ptr<RING> R = dynamic_pointer_cast<RING>(x))
//       return Value::from(ExtendedJanetBasis(runtimeEnv->evalArgAsListOfRingElem(ARG(1), R->theRing)));
//  return Value::from(ExtendedJanetBasis(runtimeEnv->evalAllArgsAsListOf<RINGELEM>(invocationExpression)));
// }


DECLARE_ARITYCHECK_FUNCTION(SyzOfGens) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(SyzOfGens) {
  invocationExpression->checkNumberOfArgs(1,2);
  long n = invocationExpression->args.size();
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<IDEAL, MODULE>(ARG(n-1), which);
  if (n==1)
    switch (which) {
    case 1: return Value::from(SyzOfGens(RefTo<ideal>(x)));
    case 2: return Value::from(SyzOfGens(RefTo<module>(x)));
    }
  // n=2
  intrusive_ptr<MODULE> M = runtimeEnv->evalArgAs<MODULE>(ARG(0));
  switch (which) {
  case 1: return Value::from(SyzOfGens(M->theModule, RefTo<ideal>(x)));
  case 2: return Value::from(SyzOfGens(M->theModule, RefTo<module>(x)));
  default: return Value::from(false); // just to keep the compiler quiet
  }
}


DECLARE_ARITYCHECK_FUNCTION(syz) { return (1<=nArg) && (nArg<=2); } // AMB 2024-03-26
DECLARE_BUILTIN_FUNCTION(syz) {
  invocationExpression->checkNumberOfArgs(1,2);
  long n = invocationExpression->args.size();
  int which;
  intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3<IDEAL, MODULE, LIST>(ARG(n-1), which);
  if (n==1)
    switch (which) {
    case 1: CoCoA_THROW_ERROR2(ERR::NYI, "minimal syzygies");
      return Value::from(NewFreeModule(RingZZ(),1));
    case 2: CoCoA_THROW_ERROR2(ERR::NYI, "minimal syzygies");
      return Value::from(NewFreeModule(RingZZ(),1));
    case 3: {
      //     if ()
      //     {
      //       vector<ModuleElem> x = runtimeEnv->evalArgAsListOf<MODULEELEM>(ARG(0));
      //       return Value::from(SyzOfGens(submodule(x)));
      //     }
      //     else
      vector<RingElem> L = runtimeEnv->evalRVAsListOfRingElem(v, ARG(0));
      return Value::from(syz(L));
    }
    default: return Value::from(false); // just to keep the compiler quiet
    }
  // n=2
  intrusive_ptr<MODULE> M = runtimeEnv->evalArgAs<MODULE>(ARG(0));
  switch (which) {
  case 1: CoCoA_THROW_ERROR2(ERR::NYI, "minimal syzygies");
    return Value::from(NewFreeModule(RingZZ(),1));
  case 2: CoCoA_THROW_ERROR2(ERR::NYI, "minimal syzygies");
    return Value::from(NewFreeModule(RingZZ(),1));
  case 3: {
    //     if ()
    //     {
    //       vector<ModuleElem> x = runtimeEnv->evalArgAsListOf<MODULEELEM>(ARG(0));
    //       return Value::from(SyzOfGens(submodule(x)));
    //     }
    //     else
    vector<RingElem> L = runtimeEnv->evalRVAsListOfRingElem(v, ARG(0));
    return Value::from(syz(M->theModule, L));
  }
  default: return Value::from(false); // just to keep the compiler quiet
  }
}


DECLARE_ARITYCHECK_FUNCTION(gcd) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(gcd) {
  invocationExpression->checkNumberOfArgs(1,2);
  if (invocationExpression->args.size()==1)
  {
    intrusive_ptr<LIST> L = runtimeEnv->evalArgAs<LIST>(ARG(0));
      if (L->size() == 0) return INT::zero;
      vector<BigInt> v1;
      if (IsVectorBigInt(v1, L))  return Value::from(gcd_forC5(v1));
      return Value::from(gcd_forC5(runtimeEnv->SpecializeToListOfRingElem(L, ARG(0))));
  }
  int which0;
  intrusive_ptr<RightValue> v0 = runtimeEnv->evalArgAsT1orT2<INT, RINGELEM>(ARG(0), which0);
  int which1;
  intrusive_ptr<RightValue> v1 = runtimeEnv->evalArgAsT1orT2<INT, RINGELEM>(ARG(1), which1);
  if (which0==1 && which1==1)
    return Value::from(gcd(RefTo<BigInt>(v0), RefTo<BigInt>(v1)));
  if (which0==1)
  {
    RingElem x1 = RefTo<RingElem>(v1);
    return Value::from(gcd(RingElem(owner(x1), RefTo<BigInt>(v0)), x1));
  }
  if (which1==1)
  {
    RingElem x0 = RefTo<RingElem>(v0);
    return Value::from(gcd(x0, RingElem(owner(x0), RefTo<BigInt>(v1))));
  }
  return Value::from(gcd(RefTo<RingElem>(v0), RefTo<RingElem>(v1)));
}


DECLARE_ARITYCHECK_FUNCTION(lcm) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(lcm)
{
  invocationExpression->checkNumberOfArgs(1,2);
  if (invocationExpression->args.size()==1)
  {
    intrusive_ptr<LIST> L = runtimeEnv->evalArgAs<LIST>(ARG(0));
      if (L->size() == 0) return INT::one;
      vector<BigInt> v1;
      if (IsVectorBigInt(v1, L))  return Value::from(lcm_forC5(v1));
      return Value::from(lcm_forC5(runtimeEnv->SpecializeToListOfRingElem(L, ARG(0))));
  }
  int which0, which1;
  intrusive_ptr<RightValue> v0 = runtimeEnv->evalArgAsT1orT2<INT, RINGELEM>(ARG(0), which0);
  intrusive_ptr<RightValue> v1 = runtimeEnv->evalArgAsT1orT2<INT, RINGELEM>(ARG(1), which1);
  if (which0==1 && which1==1)
    return Value::from(lcm(RefTo<BigInt>(v0), RefTo<BigInt>(v1)));
  if (which0==1)
  {
    RingElem x1 = RefTo<RingElem>(v1);
    return Value::from(lcm(RingElem(owner(x1), RefTo<BigInt>(v0)), x1));
  }
  if (which1==1)
  {
    RingElem x0 = RefTo<RingElem>(v0);
    return Value::from(lcm(x0, RingElem(owner(x0), RefTo<BigInt>(v1))));
  }
  return Value::from(lcm(RefTo<RingElem>(v0), RefTo<RingElem>(v1)));
}


DECLARE_ARITYCHECK_FUNCTION(RootBound) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(RootBound) {  // JAA
  invocationExpression->checkNumberOfArgs(1,2);// variable number of args
  intrusive_ptr<RINGELEM> poly = runtimeEnv->evalArgAs<RINGELEM>(ARG(0));
  const bool OneArg = (invocationExpression->args.size()==1);
  const long NumIters = OneArg?-1:runtimeEnv->evalArgAsLong(ARG(1));
  return Value::from(RootBound(poly->theRingElem, NumIters));
}


DECLARE_ARITYCHECK_FUNCTION(RootBound2) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(RootBound2) {  // JAA
  invocationExpression->checkNumberOfArgs(1,2);// variable number of args
  intrusive_ptr<RINGELEM> poly = runtimeEnv->evalArgAs<RINGELEM>(ARG(0));
  const bool OneArg = (invocationExpression->args.size()==1);
  const long NumIters = OneArg?-1:runtimeEnv->evalArgAsLong(ARG(1));
  return Value::from(RootBound2(poly->theRingElem, NumIters));
}

//---- RINGHOM -------------------------------------------------------

DECLARE_ARITYCHECK_FUNCTION(indets) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(indets) { // AMB+JAA
  invocationExpression->checkNumberOfArgs(1,2);
  intrusive_ptr<RING> R = runtimeEnv->evalArgAs<RING>(ARG(0));
  if (invocationExpression->args.size()==1)
    return Value::from(indets((runtimeEnv->evalArgAs<RING>(ARG(0)))->theRing));
  return Value::from(indets((runtimeEnv->evalArgAs<RING>(ARG(0)))->theRing,
                            runtimeEnv->evalArgAs<STRING>(ARG(1))->theString));
}


//---- MAT -------------------------------------------------------

DECLARE_ARITYCHECK_FUNCTION(JacobianMat) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(JacobianMat) {  // AMB
  invocationExpression->checkNumberOfArgs(1,2);// variable number of args
  const vector<RingElem> polys = runtimeEnv->evalArgAsListOfRingElem(ARG(0));
  if (polys.size()==0) throw RuntimeException("Empty list", ARG(0).exp);
  if (invocationExpression->args.size()==1)
    return new MAT(JacobianMat(polys, indets(owner(polys[0]))));
  // There were 2 args; 2nd arg is list of indets
  const vector<RingElem> indets = runtimeEnv->evalArgAsListOfRingElem(ARG(1));
  // Let CoCoALib fn JacobianMat check that the lists are good
  return new MAT(JacobianMat(polys, indets));
}


// variable number of args
DECLARE_ARITYCHECK_FUNCTION(MakeTermOrdMat) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(MakeTermOrdMat) {  // AMB
  invocationExpression->checkNumberOfArgs(1,2);// variable number of args
  intrusive_ptr<const MAT> M0 = runtimeEnv->evalArgAs<const MAT>(ARG(0));
  if (invocationExpression->args.size()==1)
    return Value::from(MakeTermOrdMat(M0->theMatrix));
  intrusive_ptr<INT> GrDim = runtimeEnv->evalArgAs<INT>(ARG(1));
  return Value::from(MakeTermOrdMat(M0->theMatrix, ConvertTo<long>(GrDim->theBigInt)));
}

// variable number of args
DECLARE_ARITYCHECK_FUNCTION(RandomUnimodularMat) { return (2<=nArg) && (nArg<=3); }
DECLARE_BUILTIN_FUNCTION(RandomUnimodularMat) {  // JAA
  invocationExpression->checkNumberOfArgs(2,3);// variable number of args
  intrusive_ptr<RING> R = runtimeEnv->evalArgAs<RING>(ARG(0));
  const long n = runtimeEnv->evalArgAsLong(ARG(1));
  long niters = 0;
  if (invocationExpression->args.size()==3)
  {
    niters = runtimeEnv->evalArgAsLong(ARG(2));
  }
  return Value::from(RandomUnimodularMat(R->theRing,n,niters));
}


// variable number of args
DECLARE_ARITYCHECK_FUNCTION(FrobeniusMat) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(FrobeniusMat) {  // AMB
  invocationExpression->checkNumberOfArgs(1,2);// variable number of args
  intrusive_ptr<const IDEAL> I = runtimeEnv->evalArgAs<const IDEAL>(ARG(0));
  if (invocationExpression->args.size()==1)
    return Value::from(FrobeniusMat(I->theIdeal));
  vector<RingElem> QB2 = runtimeEnv->evalArgAsListOfRingElem(RingOf(I->theIdeal), ARG(1));
  vector<PPMonoidElem> QB2pp;
  for (long i = 0; i<len(QB2); ++i)
  {
    if (IsZero(QB2[i]))  CoCoA_THROW_ERROR2(ERR::ReqNonZero, "FrobeniusMat");  
    if (!IsMonomial(QB2[i]))  CoCoA_THROW_ERROR2("expected list of PP", "FrobeniusMat");  
    QB2pp.push_back(LPP(QB2[i]));
  }
  return Value::from(FrobeniusMat(I->theIdeal, QB2pp));
}


// variable number of args
DECLARE_ARITYCHECK_FUNCTION(ideal) { return 1<=nArg; }
DECLARE_BUILTIN_FUNCTION(ideal) {
  const int nArgs = invocationExpression->args.size();
  if (nArgs==0)
    throw RuntimeException("Wrong number of arguments; found: 0, expecting: at least 1", invocationExpression);
  int which;
  if (nArgs == 1)
  {
    intrusive_ptr<const RightValue> x = runtimeEnv->evalArgAsT1orT2<LIST, RINGELEM>(ARG(0), which);
    switch (which)
    {
    case 1:
      return Value::from(ideal(runtimeEnv->evalRVAsListOfRingElem(x, ARG(0))));

    case 2:
      return Value::from(ideal(RefTo<RingElem>(x)));
    }
  }
  if (nArgs > 2)
    return Value::from(ideal(runtimeEnv->evalAllArgsAsListOf<RINGELEM>(invocationExpression)));
  // Now nArgs == 2: could be ideal(RING, LIST) or ideal(RINGELEM, RINGELEM)
  intrusive_ptr<const RightValue> x = runtimeEnv->evalArgAsT1orT2<RING, RINGELEM>(ARG(0), which);
  switch (which)
  {
  default: // never executed -- just to keep compiler quiet???
  case 1:
    return Value::from(ideal(RefTo<ring>(x), runtimeEnv->evalArgAsListOfRingElem(RefTo<ring>(x), ARG(1))));

  case 2:
    intrusive_ptr<RINGELEM> x2 = runtimeEnv->evalArgAs<RINGELEM>(ARG(1));
    return Value::from(ideal(RefTo<RingElem>(x), x2->theRingElem));
  }
  // NEVER GET HERE
  return Value::from(false); // JUST TO KEEP COMPILER QUIET!!!
  // intrusive_ptr<const RightValue> x = runtimeEnv->evalArgAsT1orT2orT3<RING, LIST, RINGELEM>(ARG(0), which);
  // if (nArgs==2 && which==1)
  //   return Value::from(ideal(RefTo<ring>(x), runtimeEnv->evalArgAsListOfRingElem(RefTo<ring>(x), ARG(1))));
  // if (nArgs==1 && which==2)
  //   return Value::from(ideal(runtimeEnv->evalRVAsListOfRingElem(x, ARG(0))));
  // return Value::from(ideal(runtimeEnv->evalAllArgsAsListOf<RINGELEM>(invocationExpression)));
}

//---- MODULE -------------------------------------------------------

DECLARE_ARITYCHECK_FUNCTION(submodule) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(submodule) {  // AMB
  invocationExpression->checkNumberOfArgs(1,2);
  int which;
  intrusive_ptr<RightValue> v0 = runtimeEnv->evalArgAsT1orT2<MODULE, LIST>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(submodule(RefTo<module>(v0), runtimeEnv->evalArgAsListOf<MODULEELEM>(ARG(1))));
  case 2: return Value::from(submodule(runtimeEnv->evalArgAsListOf<MODULEELEM>(ARG(0))));
  default: throw RuntimeException(ERRORMissingCode(v0),invocationExpression);
  }
}


//---- RING -------------------------------------------------------

DECLARE_ARITYCHECK_FUNCTION(NewPolyRing) { return (nArg>=2) && (nArg<=4); }
DECLARE_BUILTIN_FUNCTION(NewPolyRing) {
  invocationExpression->checkNumberOfArgs(2,4);  // min,max
  const int nArgs = invocationExpression->args.size();
  if (nArgs!=2 && nArgs!=4)
    throw RuntimeException("Wrong number of arguments; found: "+boost::lexical_cast<std::string>(nArgs)+", expecting: 2 or 4", invocationExpression);
  intrusive_ptr<RING> R = runtimeEnv->evalArgAs<RING>(ARG(0));
  vector<symbol> syms = runtimeEnv->evalArgAsListOfSymbols(ARG(1));
  if (nArgs==2) return Value::from(NewPolyRing(R->theRing, syms));
  intrusive_ptr<MAT> M = runtimeEnv->evalArgAs<MAT>(ARG(2));
  //  if (nArgs==3) return Value::from(NewPolyRing(R->theRing, syms, M->theMatrix)); now NewPolyRingWeights
    // nArgs == 4
  long d = runtimeEnv->evalArgAsLong(ARG(3));
  const PPOrdering PPO = NewMatrixOrdering(M->theMatrix, d);
  //return new RING(NewPolyRing(R->theRing, NewPPMonoid(syms, PPO))); //SLOW!
  return new RING(NewPolyRing(R->theRing, syms, PPO));
}


DECLARE_ARITYCHECK_FUNCTION(NewWeylAlgebra) { return (nArg==2) || (nArg==3); }
DECLARE_BUILTIN_FUNCTION(NewWeylAlgebra) {
  //  invocationExpression->checkNumberOfArgs(2,4);
  const int nArgs = invocationExpression->args.size();
  if (nArgs!=2 && nArgs!=3)
    throw RuntimeException("Wrong number of arguments; found: "+boost::lexical_cast<std::string>(nArgs)+", expecting: 2 or 3", invocationExpression);
  intrusive_ptr<RING> R = runtimeEnv->evalArgAs<RING>(ARG(0));
  vector<symbol> syms = runtimeEnv->evalArgAsListOfSymbols(ARG(1));
  if (nArgs==2)
    return Value::from(NewWeylAlgebra(R->theRing, syms, vector<long>(0)));
  if (nArgs==3)
    CoCoA_THROW_ERROR2(ERR::NYI, "builtin function NewWeylAlgebra with 3 args");
  CoCoA_THROW_ERROR2(ERR::NYI, "builtin function NewWeylAlgebra with 4+ args");
  return Value::from(long(0)); // BUG??? JUST TO KEEP COMPILER QUIET
}


DECLARE_ARITYCHECK_FUNCTION(RandomLinearForm) {return (1<=nArg) && (nArg<=2);}
DECLARE_BUILTIN_FUNCTION(RandomLinearForm) { // AMB 2018-03
  invocationExpression->checkNumberOfArgs(1,2);
  intrusive_ptr<RING> v = runtimeEnv->evalArgAs<RING>(ARG(0));
  ring P = (runtimeEnv->evalArgAs<RING>(ARG(0)))->theRing;
  if (invocationExpression->args.size()==1)
    return Value::from(RandomLinearForm(P));
  //  if (invocationExpression->args.size()==2)
  return Value::from(RandomLinearForm(P, runtimeEnv->evalArgAsLong(ARG(1))));
//   return Value::from(RandomLinearForm(P,
//                                       runtimeEnv->evalArgAsLong(ARG(1)),
//                                       runtimeEnv->evalArgAsLong(ARG(2))));
}


//---- IDEAL -------------------------------------------------------

// variable number of args // AMB 2018-03
DECLARE_ARITYCHECK_FUNCTION(MinPolyQuot) { return (nArg==3) || (nArg==4); }
DECLARE_BUILTIN_FUNCTION(MinPolyQuot) {
  invocationExpression->checkNumberOfArgs(3,4);
  const int nArgs = invocationExpression->args.size();
  intrusive_ptr<RINGELEM> f = runtimeEnv->evalArgAs<RINGELEM>(ARG(0));
  intrusive_ptr<IDEAL> I = runtimeEnv->evalArgAs<IDEAL>(ARG(1));
  intrusive_ptr<RINGELEM> x = runtimeEnv->evalArgAs<RINGELEM>(ARG(2));
  if (nArgs==3)
    return Value::from(MinPolyQuot(f->theRingElem,I->theIdeal, x->theRingElem));
  return Value::from(MinPolyQuot(f->theRingElem, I->theIdeal, x->theRingElem,
                                 VerificationLevel(runtimeEnv->evalArgAsLong(ARG(3)))));
}


DECLARE_ARITYCHECK_FUNCTION(StarRoot) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(StarRoot) {  // JAA
  invocationExpression->checkNumberOfArgs(1,2);// variable number of args
  long UPBexp = 0;
  if (invocationExpression->args.size()==2)
  {
    UPBexp = runtimeEnv->evalArgAsLong(ARG(1));
  }
  intrusive_ptr<INT> N = runtimeEnv->evalArgAs<INT>(ARG(0));
  BigInt ans = StarRoot(N->theBigInt,UPBexp);
  return Value::from(ans);
}



//---- IDEAL (points) -----------------------------------------------



//---- IDEAL (implicit -- temporary)



//------ POLY ---------------------------------------------------------

DECLARE_ARITYCHECK_FUNCTION(resultant) { return (2<=nArg) && (nArg<=3); }
DECLARE_BUILTIN_FUNCTION(resultant) {  // JAA
  invocationExpression->checkNumberOfArgs(2,3);// variable number of args
  intrusive_ptr<RINGELEM> f = runtimeEnv->evalArgAs<RINGELEM>(ARG(0));
  intrusive_ptr<RINGELEM> g = runtimeEnv->evalArgAs<RINGELEM>(ARG(1));
  if (invocationExpression->args.size()==2)
    return Value::from(resultant(f->theRingElem, g->theRingElem));
  intrusive_ptr<RINGELEM> x = runtimeEnv->evalArgAs<RINGELEM>(ARG(2));
  return Value::from(resultant_forC5(f->theRingElem, g->theRingElem, x->theRingElem));
}


DECLARE_ARITYCHECK_FUNCTION(discriminant) { return (1<=nArg) && (nArg<=2); }
DECLARE_BUILTIN_FUNCTION(discriminant) {  // JAA
  invocationExpression->checkNumberOfArgs(1,2);// variable number of args
  intrusive_ptr<RINGELEM> f = runtimeEnv->evalArgAs<RINGELEM>(ARG(0));
  if (invocationExpression->args.size()==1)
    return Value::from(discriminant(f->theRingElem));
  intrusive_ptr<RINGELEM> x = runtimeEnv->evalArgAs<RINGELEM>(ARG(1));
  return Value::from(discriminant_forC5(f->theRingElem, x->theRingElem));
}


// -- QuasiPoly -------------------------------------------------------


//-- UTILITIES ---------------------------------------------------------



} // namespace InterpreterNS
} // namespace CoCoA
