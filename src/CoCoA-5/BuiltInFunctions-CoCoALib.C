//   Copyright (c) 2010-2020 Giovanni Lagorio, John Abbott, Anna M. Bigatti
//   Authors: 2010-2011 Giovanni Lagorio
//   Authors: 2011-2020 John Abbott, Anna M. Bigatti
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
#include "Interpreter.H"  // for SpecializeToListOfRingElem
#include "CoCoALibSupplement.H"

#include <boost/foreach.hpp>

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
    {
      if (const auto n = boost::dynamic_pointer_cast<INT>(l->getValue(i)))
        v.push_back(n->theBigInt);
      else
        return false;
    }
    swap(v, BigIntVec);
    return true;
  }  

  bool IsVectorBigRat(std::vector<BigRat>& BigRatVec, const intrusive_ptr<LIST> l)
  {
    vector<BigRat> v;
    LIST::ContainerType::size_type size = l->size();
    for (unsigned long i=0; i<size; ++i)
    {
      if (const boost::intrusive_ptr<INT> n = boost::dynamic_pointer_cast<INT>(l->getValue(i)))
        v.push_back(BigRat(n->theBigInt));
      else if (const boost::intrusive_ptr<RAT> q = boost::dynamic_pointer_cast<RAT>(l->getValue(i)))
        v.push_back(q->theBigRat);
      else
        return false;
    }
    swap(v, BigRatVec);
    return true;
  }  
} // anonymous namespace


DECLARE_STD_BUILTIN_FUNCTION(RingElemList, 2)
{ // JAA 2020-11
  intrusive_ptr<RING> R = runtimeEnv->evalArgAs<RING>(ARG(0));
  intrusive_ptr<STRING> str = runtimeEnv->evalArgAs<STRING>(ARG(1));
  vector<RingElem> v = RingElemVec(R->theRing, str->theString);
  return Value::from(v);
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(reseed, 1)
{ // JAA 2017
  intrusive_ptr<INT> seed = runtimeEnv->evalArgAs<INT>(ARG(0));
  reseed_forC5(seed->theBigInt);
  return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(fold, 2)
{ // JAA 2020
  intrusive_ptr<STRING> str = runtimeEnv->evalArgAs<STRING>(ARG(0));
  intrusive_ptr<INT> width = runtimeEnv->evalArgAs<INT>(ARG(1));
  return Value::from(fold_forC5(str->theString, width->theBigInt));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(NumEvalUPoly, 2)
{ // JAA
  const BigRat q = runtimeEnv->evalArgAs<RAT>(ARG(1))->theBigRat;
  intrusive_ptr<LIST> CoeffList = runtimeEnv->evalArgAs<LIST>(ARG(0));
  LIST::ContainerType::size_type size = CoeffList->size();
  vector<BigInt> C; C.reserve(size);
  if (!IsVectorBigInt(C, CoeffList)) CoCoA_THROW_ERROR("Bad CoeffList","NumEvalUPoly");
  if (IsZero(q)) return new RAT(BigRat(C.front()));
  BigInt numer;
  HornerRecursiveIterQQ2(numer, C, num(q), den(q));
  return new RAT(BigRat(numer, power(den(q),size-1)));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(IsInteger, 2) {
  intrusive_ptr<LeftValue> resultRef = intrusive_ptr_cast<LeftValue>(runtimeEnv->evalArg(ARG(0), RuntimeEnvironment::EVAL_BY_REF));
  intrusive_ptr<RINGELEM> poly = runtimeEnv->evalArgAs<RINGELEM>(ARG(1));
  BigInt N;
  if (!IsInteger(N, poly->theRingElem)) return Value::from(false);
  if (resultRef->assignmentNeedsOwnership()) resultRef->obtainOwnership();
  const intrusive_ptr<const Expression> resultRefExp = ARG(0).exp;
  resultRef->assign(new INT(N), resultRefExp->getBegin(), resultRefExp->getEnd(), runtimeEnv);
  return Value::from(true);
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(IsRational, 2) {
  intrusive_ptr<LeftValue> resultRef = intrusive_ptr_cast<LeftValue>(runtimeEnv->evalArg(ARG(0), RuntimeEnvironment::EVAL_BY_REF));
  intrusive_ptr<RINGELEM> poly = runtimeEnv->evalArgAs<RINGELEM>(ARG(1));
  BigRat qq;
  bool isRational(IsRational(qq, poly->theRingElem));
  if (isRational) {
    if (resultRef->assignmentNeedsOwnership()) resultRef->obtainOwnership();
    const intrusive_ptr<const Expression> resultRefExp = ARG(0).exp;
    resultRef->assign(new RAT(qq), resultRefExp->getBegin(), resultRefExp->getEnd(), runtimeEnv);
  }
  return Value::from(isRational);
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(IsDivisible, 2) {
  int which0;
  int which1;
  intrusive_ptr<const RightValue> v0 = runtimeEnv->evalArgAsT1orT2<INT, RINGELEM>(ARG(0), which0);
  intrusive_ptr<const RightValue> v1 = runtimeEnv->evalArgAsT1orT2<INT, RINGELEM>(ARG(1), which1);
  switch (which0) {
  case 1: // INT
    if (which1==1)
      return Value::from(IsDivisible(RefTo<BigInt>(v0), RefTo<BigInt>(v1)));
    return Value::from(IsDivisible(RefTo<BigInt>(v0), RefTo<RingElem>(v1)));
  case 2: // RINGELEM
    if (which1==1)
      return Value::from(IsDivisible(RefTo<RingElem>(v0), RefTo<BigInt>(v1)));
    return Value::from(IsDivisible(RefTo<RingElem>(v0), RefTo<RingElem>(v1)));
  default: return Value::from(false); // just to keep the compiler quiet
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(IsDivisible_AllowFields, 2) {
  int which0;
  int which1;
  intrusive_ptr<const RightValue> v0 = runtimeEnv->evalArgAsT1orT2<INT, RINGELEM>(ARG(0), which0);
  intrusive_ptr<const RightValue> v1 = runtimeEnv->evalArgAsT1orT2<INT, RINGELEM>(ARG(1), which1);
  switch (which0) {
  case 1: // INT
    if (which1==1)
      return Value::from(IsDivisible(RefTo<BigInt>(v0), RefTo<BigInt>(v1)));
    return Value::from(IsDivisible_AllowFields(RefTo<BigInt>(v0), RefTo<RingElem>(v1)));
  case 2: // RINGELEM
    if (which1==1)
      return Value::from(IsDivisible_AllowFields(RefTo<RingElem>(v0), RefTo<BigInt>(v1)));
    return Value::from(IsDivisible_AllowFields(RefTo<RingElem>(v0), RefTo<RingElem>(v1)));
  default: return Value::from(false); // just to keep the compiler quiet
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(len, 1) {
  int which;
  intrusive_ptr<const RightValue> v = runtimeEnv->evalArgAsT1orT2orT3orT4orT5<STRING, LIST, RINGELEM, MAT, MatrixRowValue>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(BigInt(RefTo<string>(v).length()));
  case 2: return Value::from(BigInt(intrusive_ptr_cast<const LIST>(v)->size()));
  // meaningful error messages for obsolete uses of "len" ---->
  case 3: throw RuntimeException("len(RINGELEM) not allowed, use NumTerms(RINGELEM) instead", ARG(0).exp);
  case 4: throw RuntimeException("len(MAT) not allowed, use NumRows(MAT) instead", ARG(0).exp);
  case 5: throw RuntimeException("len(MATRIXROW) not allowed, use NumCols(MAT) instead", ARG(0).exp);
  default: return Value::from(false); // just to keep the compiler quiet
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(ScalarProduct, 2) {
  const Argument &arg1 = ARG(0);
  intrusive_ptr<LIST> l1 = runtimeEnv->evalArgAs<LIST>(arg1);
  const Argument &arg2 = ARG(1);
  intrusive_ptr<LIST> l2 = runtimeEnv->evalArgAs<LIST>(arg2);
  const LIST::ContainerType::size_type size = l1->size();
  const CharPointer &begin = arg1.exp->getBegin();
  const CharPointer &end = arg2.exp->getEnd();
  if (l2->size()!=size)
    throw RuntimeException("The arguments must be lists of the same size", begin, end);
  intrusive_ptr<RightValue> result = INT::zero;
  try {
    for (LIST::ContainerType::size_type a=0; a<size; ++a)
      result = runtimeEnv->binaryOperatorDispatch(
            result,
            runtimeEnv->binaryOperatorDispatch(l1->getValue(a), l2->getValue(a), RuntimeEnvironment::opStarMap, begin, end),
            RuntimeEnvironment::opPlusMap,
            begin,
            end);
  } catch (const InterruptException &) {
    throw;
  } catch (const RuntimeException &) {
    throw RuntimeException("Some elements have incompatible types for the scalar product", begin, end);
  }
  return result;
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(RandomTupleIndices, 2) {  // JAA
  const long n = runtimeEnv->evalArgAsLong(ARG(0));
  const long r = runtimeEnv->evalArgAsLong(ARG(1));
  vector<long> ans = RandomTupleIndices(n,r);
  const long card = len(ans);
  for (long i=0; i < card; ++i)
    ++ans[i];
  return Value::from(ans);
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(RandomPermutation, 1) {  // JAA
  const long n = runtimeEnv->evalArgAsLong(ARG(0));
  vector<long> ans = RandomPermutation(n);
  for (long i=0; i < n; ++i)
    ++ans[i];
  return Value::from(ans);
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(div, 2) {
  intrusive_ptr<INT> a = runtimeEnv->evalArgAs<INT>(ARG(0));
  intrusive_ptr<INT> b = runtimeEnv->evalArgAs<INT>(ARG(1));
  BigInt result;
  div(result, a->theBigInt, b->theBigInt);
  return new INT(result);
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(mod, 2) {
  intrusive_ptr<INT> a = runtimeEnv->evalArgAs<INT>(ARG(0));
  intrusive_ptr<INT> b = runtimeEnv->evalArgAs<INT>(ARG(1));
  BigInt result;
  mod(result, a->theBigInt, b->theBigInt);
  return new INT(result);
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(IsCoprime, 2) { // JAA 2023-03-06
  int which0;
  int which1;
  intrusive_ptr<const RightValue> v0 = runtimeEnv->evalArgAsT1orT2<INT, RINGELEM>(ARG(0), which0);
  intrusive_ptr<const RightValue> v1 = runtimeEnv->evalArgAsT1orT2<INT, RINGELEM>(ARG(1), which1);
  switch (which0)
  {
  case 1: // ARG(0) is INT
    if (which1==1)
      return Value::from(IsCoprime(RefTo<BigInt>(v0), RefTo<BigInt>(v1)));
    return Value::from(IsCoprime(RefTo<BigInt>(v0), RefTo<RingElem>(v1)));
  case 2: // ARG(0) is RINGELEM
  {
    if (which1==1)
      return Value::from(IsCoprime(RefTo<RingElem>(v0), RefTo<BigInt>(v1)));
    // Special handling for PPs (which CoCoA-5 reprs as polys)
    const RingElem& f = RefTo<RingElem>(v0);
    const RingElem& g = RefTo<RingElem>(v1);
    if (IsPolyRing(owner(f)) && IsMonomial(f) && IsOne(LC(f)) &&
        IsPolyRing(owner(g)) && IsMonomial(g) && IsOne(LC(g)))
      return Value::from(IsCoprime(LPP(f), LPP(g))); // power products;
    return Value::from(IsCoprime(f,g));
//    return Value::from(IsCoprime(RefTo<RingElem>(v0), RefTo<RingElem>(v1)));
  }
  default: return Value::from(false); // just to keep the compiler quiet
  }
}
END_STD_BUILTIN_FUNCTION




// DECLARE_STD_BUILTIN_FUNCTION(NumDigits, 2) {
//  intrusive_ptr<INT> a = runtimeEnv->evalArgAs<INT>(ARG(0));
//   // intrusive_ptr<INT> b = runtimeEnv->evalArgAs<INT>(ARG(1));
//  long base = runtimeEnv->evalArgAsLong(ARG(1));
//  if (base<2 || base>36)
//    throw RuntimeException("Base must be in the range 2..36", ARG(1).exp);
//  return new INT(NumDigits(a->theBigInt, base));
// }
// END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(ContFracToRat, 1) { // JAA
  const vector<BigInt> CFQuots = runtimeEnv->evalArgAsListOf<INT>(ARG(0));
    ContFracApproximant ans;
    const long n = len(CFQuots);
    for (long i=0; i < n; ++i)
    {
      if (i > 0 && sign(CFQuots[i])<=0)
        throw RuntimeException("Cont frac quotients must  be positive", ARG(0).exp);
      ans.myAppendQuot(CFQuots[i]);
    }
    return Value::from(ans.myRational());
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(CRT, 4) {
  intrusive_ptr<INT> R1 = runtimeEnv->evalArgAs<INT>(ARG(0));
  intrusive_ptr<INT> M1 = runtimeEnv->evalArgAs<INT>(ARG(1));
  intrusive_ptr<INT> R2 = runtimeEnv->evalArgAs<INT>(ARG(2));
  intrusive_ptr<INT> M2 = runtimeEnv->evalArgAs<INT>(ARG(3));
        CRTMill CRT;
        CRT.myAddInfo(R1->theBigInt, M1->theBigInt);
        CRT.myAddInfo(R2->theBigInt, M2->theBigInt);

        // Create return value: ... record
        intrusive_ptr<RECORD> ans(new RECORD);
        ans->setFieldNoCheck("residue", Value::from(CombinedResidue(CRT)));
        ans->setFieldNoCheck("modulus", Value::from(CombinedModulus(CRT)));
        return ans;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(RatReconstructWithBounds, 5) {
  intrusive_ptr<INT> E = runtimeEnv->evalArgAs<INT>(ARG(0));
  intrusive_ptr<INT> P = runtimeEnv->evalArgAs<INT>(ARG(1));
  intrusive_ptr<INT> Q = runtimeEnv->evalArgAs<INT>(ARG(2));
  const vector<BigInt> ResList = runtimeEnv->evalArgAsListOf<INT>(ARG(3));
  const vector<BigInt> ModList = runtimeEnv->evalArgAsListOf<INT>(ARG(4));

        const long NumMod = len(ModList);
        const long NumRes = len(ResList);
        if (NumMod != NumRes)
          CoCoA_THROW_ERROR2(ERR:IncompatDims,  "residue list and modulus list");
        const long e = ConvertTo<long>(E->theBigInt);
        vector<long> mod(NumMod);
        for (long i=0; i < NumMod; ++i)
          mod[i] = ConvertTo<long>(ModList[i]);
        vector<long> res(NumRes);
        for (long i=0; i < NumRes; ++i)
          res[i] = ConvertTo<long>(ResList[i]%mod[i]);
        const BigRat result = RatReconstructWithBounds(e, P->theBigInt, Q->theBigInt, res, mod);
        // Create return value: ... record
        intrusive_ptr<RECORD> ans(new RECORD);
        ans->setFieldNoCheck("failed", Value::from(den(result) == 0));
        if (den(result) != 0)
          ans->setFieldNoCheck("ReconstructedRat", Value::from(result));
        return ans;
}
END_STD_BUILTIN_FUNCTION



DECLARE_STD_BUILTIN_FUNCTION(MantissaAndExponent10, 2)
{
  const long n = runtimeEnv->evalArgAsLong(ARG(1));
  if (n < 1) throw RuntimeException("Precision must be positive (and not too large)", ARG(1).exp);

  MantExp10 ME;
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<INT, RAT>(ARG(0), which);
  switch (which)
  {
  case 1: ME = MantissaAndExponent10(RefTo<BigInt>(x), n); break;
  case 2: ME = MantissaAndExponent10(RefTo<BigRat>(x), n); break;
  default: return Value::from(false); // just to keep the compiler quiet
  }

  // Create return value: ... record
  intrusive_ptr<RECORD> ans(new RECORD);
  ans->setFieldNoCheck("mantissa", Value::from(ME.mySign * ME.myMantissa));
  ans->setFieldNoCheck("exponent", Value::from(ME.myExponent));
  ans->setFieldNoCheck("NumDigits", Value::from(ME.myNumDigits));
  return ans;
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(FloatApprox, 2)
{
  int which;
  const auto x = runtimeEnv->evalArgAsT1orT2orT3<INT, RAT, RINGELEM>(ARG(0), which);
  BigRat q;
  switch (which)
  {
  case 1: q = RefTo<BigInt>(x); break;
  case 2: q = RefTo<BigRat>(x); break;
  case 3: q = ConvertTo<BigRat>(RefTo<RingElem>(x)); break;
  default: return Value::from(false); // just to keep the compiler quiet
  }
  // body
  long PrecBits = runtimeEnv->evalArgAsLong(ARG(1));
  if (PrecBits < 2) throw RuntimeException("nonsensical value", ARG(1).exp);
  return Value::from(FloatApprox(q, PrecBits));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(sign, 1) {
  int which;
  const auto v = runtimeEnv->evalArgAsT1orT2orT3<INT, RAT, RINGELEM>(ARG(0), which);
  switch (which) {
  case 1: return INT::fromInt(sign(RefTo<BigInt>(v)));
  case 2: return INT::fromInt(sign(RefTo<BigRat>(v)));
  case 3: return INT::fromInt(sign(RefTo<RingElem>(v)));
  default: return Value::from(false); // just to keep the compiler quiet
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(IsZero, 1) {
  int which;
  const auto v = runtimeEnv->evalArgAsT1orT2orT3orT4orT5orT6orT7<INT, RAT, RINGELEM, MODULEELEM, IDEAL, MAT, MODULE>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(IsZero(RefTo<BigInt>(v)));
  case 2: return Value::from(IsZero(RefTo<BigRat>(v)));
  case 3: return Value::from(IsZero(RefTo<RingElem>(v)));
  case 4: return Value::from(IsZero(RefTo<ModuleElem>(v)));
  case 5: return Value::from(IsZero(RefTo<ideal>(v)));
  case 6: return Value::from(IsZero(RefTo<matrix>(v)));
  case 7: return Value::from(IsZero(RefTo<module>(v)));
  default: return Value::from(false); // just to keep the compiler quiet
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(IsOne, 1) {
  int which;
  const auto v = runtimeEnv->evalArgAsT1orT2orT3orT4<INT, RAT, RINGELEM, IDEAL>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(IsOne(RefTo<BigInt>(v)));
  case 2: return Value::from(IsOne(RefTo<BigRat>(v)));
  case 3: return Value::from(IsOne(RefTo<RingElem>(v)));
  case 4: return Value::from(IsOne(RefTo<ideal>(v)));
  default: return Value::from(false); // just to keep the compiler quiet
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(IsMinusOne, 1) {
  int which;
  const auto v = runtimeEnv->evalArgAsT1orT2orT3<INT, RAT, RINGELEM>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(IsMinusOne(RefTo<BigInt>(v)));
  case 2: return Value::from(IsMinusOne(RefTo<BigRat>(v)));
  case 3: return Value::from(IsMinusOne(RefTo<RingElem>(v)));
  default: return Value::from(false); // just to keep the compiler quiet
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(binomial, 2) {
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<RINGELEM, INT>(ARG(0), which);
  intrusive_ptr<INT> n = runtimeEnv->evalArgAs<INT>(ARG(1));
  switch (which) {
  case 1: return Value::from(binomial(RefTo<RingElem>(x), n->theBigInt));
  case 2: return Value::from(binomial(RefTo<BigInt>(x), n->theBigInt));
  default: return Value::from(false); // just to keep the compiler quiet
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(radical, 1) {
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2orT3<RINGELEM, INT, IDEAL>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(radical(RefTo<RingElem>(x)));
  case 2: return Value::from(radical(RefTo<BigInt>(x)));
  case 3: return Value::from(radical(RefTo<ideal>(x)));
  default: return Value::from(false); // just to keep the compiler quiet
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(IsSqFree, 1) {
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<RINGELEM, INT>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(IsSqFree(RefTo<RingElem>(x)));
  case 2: return Value::from(!IsFalse3(IsSqFree(RefTo<BigInt>(x))));  // !!! BUG BUG BUG should handle uncertain3 properly!  BUG BUG BUG
  default: return Value::from(false); // just to keep the compiler quiet
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(zero, 1) {
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<RING, MODULE>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(zero(RefTo<ring>(x)));
  case 2: return Value::from(zero(RefTo<module>(x)));
  default: return Value::from(false); // just to keep the compiler quiet
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(gens, 1) {
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<IDEAL, MODULE>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(gens(RefTo<ideal>(x)));
  case 2: return Value::from(gens(RefTo<module>(x)));
  default: return Value::from(false); // just to keep the compiler quiet
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(NumGens, 1) {
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<IDEAL, MODULE>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(NumGens(RefTo<ideal>(x)));
  case 2: return Value::from(NumGens(RefTo<module>(x)));
  default: return Value::from(false); // just to keep the compiler quiet
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(IsContained, 2) {
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<IDEAL, MODULE>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(IsContained(RefTo<ideal>(x), runtimeEnv->evalArgAs<IDEAL>(ARG(1))->theIdeal));
  case 2: return Value::from(IsContained(RefTo<module>(x), runtimeEnv->evalArgAs<MODULE>(ARG(1))->theModule));
  default: return Value::from(false); // just to keep the compiler quiet
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(IsElem, 2) {
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<RINGELEM, MODULEELEM>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(IsElem(RefTo<RingElem>(x), runtimeEnv->evalArgAs<IDEAL>(ARG(1))->theIdeal));
  case 2: return Value::from(IsElem(RefTo<ModuleElem>(x), runtimeEnv->evalArgAs<MODULE>(ARG(1))->theModule));
  default: return Value::from(false); // just to keep the compiler quiet
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(IsInRadical, 2)
{
  int which;
  intrusive_ptr<RightValue> arg0 = runtimeEnv->evalArgAsT1orT2<RINGELEM, IDEAL>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(IsInRadical(RefTo<RingElem>(arg0), runtimeEnv->evalArgAs<IDEAL>(ARG(1))->theIdeal));
  case 2: return Value::from(IsInRadical(RefTo<ideal>(arg0), runtimeEnv->evalArgAs<IDEAL>(ARG(1))->theIdeal));
  default: return Value::from(false); // just to keep the compiler quiet
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(GBasis, 1) {
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<IDEAL, MODULE>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(GBasis(RefTo<ideal>(x)));
  case 2: return Value::from(TidyGens(RefTo<module>(x)));
  default: return Value::from(false); // just to keep the compiler quiet
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(GBasisTrunc, 2) {
  intrusive_ptr<IDEAL> a = runtimeEnv->evalArgAs<IDEAL>(ARG(0));
  long n = runtimeEnv->evalArgAsLong(ARG(1));
  return Value::from(GBasisTrunc(a->theIdeal, n));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(GBasis_timeout, 2) {
  int which;
  intrusive_ptr<RightValue> I = runtimeEnv->evalArgAsT1orT2<IDEAL, MODULE>(ARG(0), which);
  int INTorRAT;
  intrusive_ptr<RightValue> TimeLimit = runtimeEnv->evalArgAsT1orT2<INT, RAT>(ARG(1), INTorRAT);
  double Tmax;
  switch (INTorRAT)
  {
  case 1: Tmax = ConvertTo<double>(RefTo<BigInt>(TimeLimit)); break;
  case 2: Tmax = ConvertTo<double>(RefTo<BigRat>(TimeLimit)); break;
  default: return Value::from(false); // just to keep the compiler quiet
  }
  switch (which) {
  case 1: return Value::from(GBasis(RefTo<ideal>(I), CpuTimeLimit(Tmax,IterationVariability::high)));
  case 2: return Value::from(TidyGens(RefTo<module>(I)));
  default: return Value::from(false); // just to keep the compiler quiet
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(MinGens, 1) {
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<IDEAL, MODULE>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(MinGens(RefTo<ideal>(x)));
  case 2: return Value::from(MinGens(RefTo<module>(x)));
  default: return Value::from(false); // just to keep the compiler quiet
  }
}
END_STD_BUILTIN_FUNCTION


// DECLARE_STD_BUILTIN_FUNCTION(syz, 1) --> moved to BuiltInFunctionsVarArgs-CoCoALib.C AMB 2024-03-26


DECLARE_STD_BUILTIN_FUNCTION(interreduced, 1)
{ // JAA
  const vector<RingElem> v = runtimeEnv->evalArgAsListOfRingElem(ARG(0));
  return Value::from(interreduced(v));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(RingElem, 2) {
  intrusive_ptr<const RING> R(runtimeEnv->evalArgAs<RING>(ARG(0)));
  int which;
  intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3orT4orT5<INT, RAT, STRING, RINGELEM, LIST>(ARG(1), which);
  switch (which) {
  case 1: return Value::from(RingElem(R->theRing, RefTo<BigInt>(v)));
  case 2: return Value::from(RingElem(R->theRing, RefTo<BigRat>(v)));
  case 3: return Value::from(RingElem(R->theRing, RefTo<string>(v)));
  case 4: return Value::from(RingElem(R->theRing, RefTo<RingElem>(v)));
  case 5: {// Expecting list of [STRING, INT, INT, ..., INT] specifying a symbol
    const intrusive_ptr<const LIST> l = dynamic_pointer_cast<const LIST>(v);
    LIST::ContainerType::size_type ListLen = l->size();
    string SymbolHead;
    if (ListLen == 0) throw RuntimeException("Non-empty list expected", ARG(1).exp);
    if ( const boost::intrusive_ptr<STRING> s = boost::dynamic_pointer_cast<STRING>(l->getValue(0)))
      SymbolHead = s->theString;
    else
      throw RuntimeException("String expected as first entry in list", ARG(1).exp);
    std::vector<BigInt> indices;
    for (LIST::ContainerType::size_type a=1; a<ListLen; ++a)
    {
      if (const boost::intrusive_ptr<INT> elem = boost::dynamic_pointer_cast<INT>(l->getValue(a)))
        indices.push_back(theValue(elem));
      else
        throw RuntimeException("Integer expected for symbol index", ARG(1).exp);
    }
    return Value::from(RingElem(R->theRing, symbol(SymbolHead, VectorLong(indices, "RingElem"))));
  }
  default: return Value::from(false); // just to keep the compiler quiet
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(homog, 2) {
  int which;
  intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3orT4orT5<RINGELEM, MODULEELEM, IDEAL, MODULE, LIST>(ARG(0), which);
  intrusive_ptr<const RINGELEM> h(runtimeEnv->evalArgAs<RINGELEM>(ARG(1)));
  switch (which)
  {
  case 1: return Value::from(homog(RefTo<RingElem>(v), h->theRingElem));
  case 2: return Value::from(homog(RefTo<ModuleElem>(v), h->theRingElem));
  case 3: return Value::from(homog(RefTo<ideal>(v), h->theRingElem));
  case 4: CoCoA_THROW_ERROR(ERR::NYI, "homog(module)"); break;
//return Value::from(homog(intrusive_ptr_cast<MODULE>(v)->theModule, h->theRingElem));
  case 5:
  {
    const vector<RingElem> v1 = runtimeEnv->evalRVAsListOfRingElem(v, ARG(0));
    intrusive_ptr<LIST> returnValue = new LIST();
    for (long a=0; a<len(v1); ++a)
      returnValue->addValue(Value::from(homog(v1[a], h->theRingElem)));
    return returnValue;
  }
  }
  throw RuntimeException(ERRORMissingCode(v),invocationExpression);
//  CoCoA_THROW_ERROR(ERR::ShouldNeverGetHere, "BuiltinFunction-homog");
//  return Value::from(false); // just to keep the compiler quiet
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(RandomSmallPrime, 1) {
  const BigInt N = runtimeEnv->evalArgAs<INT>(ARG(0))->theBigInt;
  const char* const FnName = "RandomSmallPrime";
  long n;
  if (!IsConvertible(n, N))
    CoCoA_THROW_ERROR(ERR::ArgTooBig, FnName);
  if (n < 0)
    CoCoA_THROW_ERROR(ERR::ReqNonNegative, FnName);
  return Value::from(RandomSmallPrime(n)); // might throw -- checks arg
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(RandomNBitPrime, 1) {
  const BigInt N = runtimeEnv->evalArgAs<INT>(ARG(0))->theBigInt;
  const char* const FnName = "RandomNBitPrime";
  long n;
  if (!IsConvertible(n, N))
    CoCoA_THROW_ERROR(ERR::ArgTooBig, FnName);
  if (n < 0)
    CoCoA_THROW_ERROR(ERR::ReqNonNegative, FnName);
  return Value::from(RandomNBitPrime(n)); // might throw -- checks arg
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(abs, 1) {
  int which;
  intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3<INT, RAT, RINGELEM>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(abs(RefTo<BigInt>(v)));
  case 2: return Value::from(abs(RefTo<BigRat>(v)));
  case 3: return Value::from(abs(RefTo<RingElem>(v)));
  default: return Value::from(false); // just to keep the compiler quiet
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(floor, 1) {
  int which;
  intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2<INT, RAT>(ARG(0), which);
  switch (which) {
  case 1: return v;
  case 2: return Value::from(floor(RefTo<BigRat>(v)));
  default: return Value::from(false); // just to keep the compiler quiet
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(ceil, 1) {
  int which;
  intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2<INT, RAT>(ARG(0), which);
  switch (which) {
  case 1: return v;
  case 2: return Value::from(ceil(RefTo<BigRat>(v)));
  default: return Value::from(false); // just to keep the compiler quiet
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(round, 1) {
  int which;
  intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2<INT, RAT>(ARG(0), which);
  switch (which) {
  case 1: return v;
  case 2: return Value::from(round(RefTo<BigRat>(v)));
  default: return Value::from(false); // just to keep the compiler quiet
  }
}
END_STD_BUILTIN_FUNCTION



DECLARE_STD_BUILTIN_FUNCTION(ExtGcd, 1) // TEMPORARY JAA
{
  intrusive_ptr<LIST> l = runtimeEnv->evalArgAs<LIST>(ARG(0));
//  if (l->size() == 0) return l;
  vector<BigInt> VecBigInt;
  if (IsVectorBigInt(VecBigInt, l))  return Value::from(ExtGcd(VecBigInt).myCofacs);
  CoCoA_THROW_ERROR("Expected list of integers","ExtGcd");
  return Value::from(false); // just to keep compiler quiet
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(CoprimeFactorBasis, 1)
{
  intrusive_ptr<LIST> L = runtimeEnv->evalArgAs<LIST>(ARG(0));
  if (L->size() == 0) return L;
  vector<BigInt> VecBigInt;
  if (IsVectorBigInt(VecBigInt, L))  return Value::from(CoprimeFactorBasis_forC5(VecBigInt));
  return Value::from(CoprimeFactorBasis_forC5(runtimeEnv->SpecializeToListOfRingElem(L, ARG(0))));
//runtimeEnv->evalArgAsListOfRingElem(ARG(0))));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(num, 1) { // AMB (changed)
  const Argument &arg = ARG(0);
  int which;
  intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3<INT, RAT, RINGELEM>(arg, which);
  switch (which) {
  case 1: return v;
  case 2: return Value::from(num(RefTo<BigRat>(v)));
  case 3: return Value::from(num(RefTo<RingElem>(v)));
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(den, 1) {
  const Argument &arg = ARG(0);
  int which;
  intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3<INT, RAT, RINGELEM>(arg, which);
  switch (which) {
  case 1: return Value::from(BigInt(1));
  case 2: return Value::from(den(RefTo<BigRat>(v)));
  case 3: return Value::from(den(RefTo<RingElem>(v)));
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(SimplestRatBetween, 2) { //JAA 2012-12-11
  const Argument &arg0 = ARG(0);
  const Argument &arg1 = ARG(1);
        BigRat A, B;
        int which;
  intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2<INT, RAT>(arg0, which);
  switch (which) {
  case 1: A = RefTo<BigInt>(v); break;
  case 2: A = RefTo<BigRat>(v); break;
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
  v = runtimeEnv->evalArgAsT1orT2<INT, RAT>(arg1, which);
  switch (which) {
  case 1: B = RefTo<BigInt>(v); break;
  case 2: B = RefTo<BigRat>(v); break;
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
  return new RAT(SimplestBigRatBetween(A, B));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(SimplestBinaryRatBetween, 2) { //JAA 2012-12-11
  const Argument &arg0 = ARG(0);
  const Argument &arg1 = ARG(1);
        BigRat A, B;
        int which;
  intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2<INT, RAT>(arg0, which);
  switch (which) {
  case 1: A = RefTo<BigInt>(v); break;
  case 2: A = RefTo<BigRat>(v); break;
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
  v = runtimeEnv->evalArgAsT1orT2<INT, RAT>(arg1, which);
  switch (which) {
  case 1: B = RefTo<BigInt>(v); break;
  case 2: B = RefTo<BigRat>(v); break;
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
  return new RAT(SimplestBinaryRatBetween(A, B));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(FloorLog2, 1) { // JAA
  int which;
  intrusive_ptr<RightValue> v0 = runtimeEnv->evalArgAsT1orT2<INT, RAT>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(FloorLog2(RefTo<BigInt>(v0)));
  case 2: return Value::from(FloorLog2(RefTo<BigRat>(v0)));
  default: throw RuntimeException(ERRORMissingCode(v0),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(FloorLog10, 1) { // JAA
  int which;
  intrusive_ptr<RightValue> v0 = runtimeEnv->evalArgAsT1orT2<INT, RAT>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(FloorLog10(RefTo<BigInt>(v0)));
  case 2: return Value::from(FloorLog10(RefTo<BigRat>(v0)));
  default: throw RuntimeException(ERRORMissingCode(v0),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(FloorLogBase, 2) { // AMB
  intrusive_ptr<INT> base = runtimeEnv->evalArgAs<INT>(ARG(1));
  int which;
  intrusive_ptr<RightValue> v0 = runtimeEnv->evalArgAsT1orT2<INT, RAT>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(FloorLogBase(RefTo<BigInt>(v0), base->theBigInt));
  case 2: return Value::from(FloorLogBase(RefTo<BigRat>(v0), base->theBigInt));
  default: throw RuntimeException(ERRORMissingCode(v0),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION





DECLARE_STD_BUILTIN_FUNCTION(exponents, 1) { // AMB
  intrusive_ptr<RINGELEM> f = runtimeEnv->evalArgAs<RINGELEM>(ARG(0));
  vector<BigInt> expv;
  if (NumIndets(owner(f->theRingElem))==1)
    expv.push_back(BigInt(deg(f->theRingElem)));
  else 
    BigExponents(expv, LPP(f->theRingElem));
  return new LIST(expv);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(wdeg, 1) { // AMB
  int which;
  intrusive_ptr<const RightValue> x = runtimeEnv->evalArgAsT1orT2<RINGELEM, MODULEELEM>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(wdeg(RefTo<RingElem>(x)));
  case 2: return Value::from(wdeg(RefTo<ModuleElem>(x)));
  default: throw RuntimeException(ERRORMissingCode(x),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(IsHomog, 1) { // AMB
  int which;
  intrusive_ptr<const RightValue> v = runtimeEnv->evalArgAsT1orT2orT3orT4orT5<RINGELEM, MODULEELEM, IDEAL, MODULE, LIST>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(IsHomog(RefTo<RingElem>(v)));
  case 2: return Value::from(IsHomog(RefTo<ModuleElem>(v)));
  case 3: return Value::from(IsHomog(RefTo<ideal>(v)));
  case 4: return Value::from(IsHomog(RefTo<module>(v)));
  case 5: {
    const intrusive_ptr<const LIST> l = intrusive_ptr_cast<const LIST>(v);
    LIST::ContainerType::size_type size = l->size();
    for (LIST::ContainerType::size_type a=0; a<size; ++a)
    {
      if (const boost::intrusive_ptr<const RINGELEM> f = boost::dynamic_pointer_cast<const RINGELEM>(l->getValue(a)))
      { if (!IsHomog(f->theRingElem)) return Value::from(false); }
      else throw RuntimeException("RingElem expected", ARG(1).exp);
    }
    return Value::from(true);
  }
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(IsFactorClosed, 1)
{ // JAA
  const vector<RingElem> v = runtimeEnv->evalArgAsListOfRingElem(ARG(0));
  return Value::from(IsFactorClosed_forC5(v));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(CommonDenom, 1)
{ // AMB 2018-03
  int which;
  intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2<RINGELEM, LIST>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(CommonDenom(RefTo<RingElem>(v)));
  case 2: return Value::from(CommonDenom(runtimeEnv->evalRVAsListOfRingElem(v, ARG(0))));
  }
  return Value::from(false); // just to keep the compiler quiet
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(LPP, 1) { // AMB
  int which;
  intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2<RINGELEM, MODULEELEM>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(LPP_forC5(RefTo<RingElem>(v)));
  case 2: return Value::from(LPP_forC5(RefTo<ModuleElem>(v)));
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(LC, 1) { // AMB
  int which;
  intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2<RINGELEM, MODULEELEM>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(LC(RefTo<RingElem>(v)));
  case 2: return Value::from(LC(RefTo<ModuleElem>(v)));
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(LM, 1) { // AMB
  int which;
  intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2<RINGELEM, MODULEELEM>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(LM_forC5(RefTo<RingElem>(v)));
  case 2: return Value::from(LM_forC5(RefTo<ModuleElem>(v)));
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(LT, 1) { // AMB
  int which;
  intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3orT4<RINGELEM, MODULEELEM, IDEAL, MODULE>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(LT_forC5(RefTo<RingElem>(v)));
  case 2: return Value::from(LT_forC5(RefTo<ModuleElem>(v)));
  case 3: return Value::from(LT(RefTo<ideal>(v)));
  case 4: return Value::from(LT(RefTo<module>(v)));
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(LF, 1) { // AMB
  int which;
  intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3orT4<RINGELEM, MODULEELEM, IDEAL, MODULE>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(LF(RefTo<RingElem>(v)));
  case 2: CoCoA_THROW_ERROR(ERR::NYI,"LF for MODULEELEM");
    return Value::from(false); // just to keep the compiler quiet
    //return Value::from(LF(RefTo<ModuleElem>(v)));
  case 3: return Value::from(LF(RefTo<ideal>(v)));
  case 4: CoCoA_THROW_ERROR(ERR::NYI,"LF for MODULE");
    //return Value::from(LF(RefTo<module>(v)));
    return Value::from(false); // just to keep the compiler quiet
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(RingOf, 1) { // AMB
  int which;
  intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3orT4orT5<RINGELEM,IDEAL,MAT,MODULE,MODULEELEM>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(owner(RefTo<RingElem>(v)));
  case 2: return Value::from(RingOf(RefTo<ideal>(v)));
  case 3: return Value::from(RingOf(RefTo<matrix>(v)));
  case 4: return Value::from(RingOf(RefTo<module>(v)));
  case 5: return Value::from(RingOf(owner(RefTo<ModuleElem>(v))));
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(ModuleOf, 1) { // AMB
  int which;
  intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2<MODULEELEM, MODULE>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(owner(RefTo<ModuleElem>(v)));
  case 2: return Value::from(AmbientFreeModule(RefTo<module>(v)));
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(CRTPoly, 4) { // AMB
  intrusive_ptr<RINGELEM> f1 = runtimeEnv->evalArgAs<RINGELEM>(ARG(0));
  intrusive_ptr<INT> M1 = runtimeEnv->evalArgAs<INT>(ARG(1));
  intrusive_ptr<RINGELEM> f2 = runtimeEnv->evalArgAs<RINGELEM>(ARG(2));
  intrusive_ptr<INT> M2 = runtimeEnv->evalArgAs<INT>(ARG(3));
  RingElem f;
  BigInt M;
  CRTPoly(f,M,  f1->theRingElem,M1->theBigInt,  f2->theRingElem,M2->theBigInt);
  // Create return value: ... record
  intrusive_ptr<RECORD> ans(new RECORD);
  ans->setFieldNoCheck("residue", Value::from(f));
  ans->setFieldNoCheck("modulus", Value::from(M));
  return ans;
}
END_STD_BUILTIN_FUNCTION


//---- RINGHOM -------------------------------------------------------

DECLARE_STD_BUILTIN_FUNCTION(InducedHom, 2) { // AMB
  intrusive_ptr<RING> R = runtimeEnv->evalArgAs<RING>(ARG(0));
  intrusive_ptr<RINGHOM> phi = runtimeEnv->evalArgAs<RINGHOM>(ARG(1));
  if (IsFractionField(R->theRing))
    return Value::from(InducedHom(FractionField(R->theRing),phi->theRingHom));
  if (IsQuotientRing(R->theRing))
    return Value::from(InducedHom(QuotientRing(R->theRing), phi->theRingHom));
  throw RuntimeException("FractionField or QuotientRing expected", ARG(0).exp);
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(PolyAlgebraHom, 3) { // AMB
  intrusive_ptr<RING> P = runtimeEnv->evalArgAs<RING>(ARG(0));
  intrusive_ptr<RING> R = runtimeEnv->evalArgAs<RING>(ARG(1));
  int which;
  intrusive_ptr<const RightValue> v = runtimeEnv->evalArgAsT1orT2<LIST,STRING>(ARG(2), which);
  switch (which) {
  case 1: return Value::from(PolyAlgebraHom(P->theRing, R->theRing,
                                            runtimeEnv->evalRVAsListOfRingElem(R->theRing, v, ARG(2))));
  case 2: return Value::from(PolyAlgebraHom(P->theRing, R->theRing, RefTo<string>(v)));
  default: throw RuntimeException("second arg must be LIST/STRING",ARG(2).exp);
  }
  return Value::from(RingZZ()); // just to keep the compiler quiet
}
END_STD_BUILTIN_FUNCTION 


DECLARE_STD_BUILTIN_FUNCTION(PolyRingHom, 4) { // AMB
  intrusive_ptr<RING> P = runtimeEnv->evalArgAs<RING>(ARG(0));
  intrusive_ptr<RING> R = runtimeEnv->evalArgAs<RING>(ARG(1));
  intrusive_ptr<RINGHOM> CoeffHom = runtimeEnv->evalArgAs<RINGHOM>(ARG(2));
  int which;
  intrusive_ptr<const RightValue> v = runtimeEnv->evalArgAsT1orT2<LIST,STRING>(ARG(3), which);
  switch (which) {
  case 1: return Value::from(PolyRingHom(P->theRing, R->theRing, CoeffHom->theRingHom,
                                         runtimeEnv->evalRVAsListOfRingElem(R->theRing, v, ARG(3))));
  case 2: return Value::from(PolyRingHom(P->theRing, R->theRing, CoeffHom->theRingHom, RefTo<string>(v)));
  default: throw RuntimeException("second arg must be LIST/STRING",ARG(3).exp);
  }
  return Value::from(RingZZ()); // just to keep the compiler quiet

}
END_STD_BUILTIN_FUNCTION


// DECLARE_STD_BUILTIN_FUNCTION(apply, 2) {
//   intrusive_ptr<RINGHOM> phi = runtimeEnv->evalArgAs<RINGHOM>(ARG(0));
//   int which;
//   intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3<RINGELEM, MAT, LIST>(ARG(1), which);
//   switch (which) {
// //  case 1: return Value::from(phi->theRingHom(RefTo<RingElem>(v)));
//   case 1: return Value::from(phi->theRingHom(RefTo<RingElem>(v)));
//   case 2: return Value::from(phi->theRingHom(RefTo<matrix>(v)));
//   case 3: {
//     vector<RingElem> w = runtimeEnv->evalRVAsListOfRingElem(v, ARG(1));
//     return Value::from(phi->theRingHom(w));
//   }
//   default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
//   }
// }
// END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(IsInImage, 2) { // JAA
  intrusive_ptr<RINGHOM> phi = runtimeEnv->evalArgAs<RINGHOM>(ARG(0));
  int which;
  intrusive_ptr<RightValue> arg = runtimeEnv->evalArgAsT1orT2orT3<INT, RAT, RINGELEM>(ARG(1), which);
  RingElem y(codomain(phi->theRingHom));
  switch (which) {
  case 1: y = RefTo<BigInt>(arg); break;
  case 2: y = RefTo<BigRat>(arg); break;
  case 3: y = RefTo<RingElem>(arg); break;
  default: throw RuntimeException(ERRORMissingCode(arg), invocationExpression);
  }
  return Value::from(IsInImage(phi->theRingHom, y));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(preimage0, 2) { // JAA
  intrusive_ptr<RINGHOM> phi = runtimeEnv->evalArgAs<RINGHOM>(ARG(0));
  int which;
  intrusive_ptr<RightValue> arg = runtimeEnv->evalArgAsT1orT2orT3<INT, RAT, RINGELEM>(ARG(1), which);
  RingElem y(codomain(phi->theRingHom));
  switch (which) {
  case 1: y = RefTo<BigInt>(arg); break;
  case 2: y = RefTo<BigRat>(arg); break;
  case 3: y = RefTo<RingElem>(arg); break;
  default: throw RuntimeException(ERRORMissingCode(arg), invocationExpression);
  }
  return Value::from(preimage0(phi->theRingHom, y));
}
END_STD_BUILTIN_FUNCTION



//---- INDET -------------------------------------------------------

DECLARE_STD_BUILTIN_FUNCTION(indet, 2) { // AMB
  intrusive_ptr<RING> a = runtimeEnv->evalArgAs<RING>(ARG(0));
  //  intrusive_ptr<INT> b = runtimeEnv->evalArgAs<INT>(ARG(1));
  long n = runtimeEnv->evalArgAsLong(ARG(1));
  //  if (!IsConvertible(n, b->theBigInt))
  //    throw RuntimeException("invalid indet index", ARG(1).exp);
  return Value::from(indet(a->theRing, n-1));
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(IndetIndex, 1) { // AMB
  long i;
  if (!IsIndet(i, runtimeEnv->evalArgAs<RINGELEM>(ARG(0))->theRingElem))
    throw RuntimeException("not an indet", ARG(0).exp);
  return Value::from(i+1);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(IndetName, 1) { // AMB
  intrusive_ptr<RINGELEM> a = runtimeEnv->evalArgAs<RINGELEM>(ARG(0));
  long i;
  if (!IsIndet(i,a->theRingElem))
    throw RuntimeException("not an indet", ARG(0).exp);
  const symbol s=IndetSymbol(owner(a->theRingElem), i);
  return Value::from(head(s));
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(IndetSubscripts, 1) { // AMB
  intrusive_ptr<RINGELEM> a = runtimeEnv->evalArgAs<RINGELEM>(ARG(0));
  long i;
  if (!IsIndet(i,a->theRingElem))
    throw RuntimeException("not an indet", ARG(0).exp);
  const symbol s=IndetSymbol(owner(a->theRingElem), i);
  vector<long> iss;
  for (long n=0; n<NumSubscripts(s); ++n) iss.push_back(subscript(s,n));
  return Value::from(iss);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(IndetsProd, 1) {
  int which;
  intrusive_ptr<const RightValue> v = runtimeEnv->evalArgAsT1orT2<RINGELEM, LIST>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(IndetsProd(RefTo<RingElem>(v)));
  case 2: {
    vector<RingElem> x = runtimeEnv->evalRVAsListOfRingElem(v, ARG(0));
    return Value::from(IndetsProd(x));
  }
  default: return Value::from(false); // just to keep the compiler quiet
  }
}
END_STD_BUILTIN_FUNCTION



DECLARE_STD_BUILTIN_FUNCTION(ContentWRT, 2) { // AMB
  intrusive_ptr<RINGELEM> f = runtimeEnv->evalArgAs<RINGELEM>(ARG(0));
  int which;
  intrusive_ptr<RightValue> idt = runtimeEnv->evalArgAsT1orT2<RINGELEM, LIST>(ARG(1), which);
  switch (which) {
  case 1: return Value::from(ContentWRT(f->theRingElem, RefTo<RingElem>(idt)));
  case 2: {
//     vector<long> v1 = VectorLongDecr1(runtimeEnv->evalArgAsListOf<INT>(ARG(1)), ERR::BadIndetIndex, "ContentWRT");
//     return Value::from(ContentWRT(f->theRingElem, v1));
    vector<RingElem> v1 = runtimeEnv->evalRVAsListOfRingElem(idt, ARG(1));
    return Value::from(ContentWRT_forC5(f->theRingElem, v1));
  }
  default: throw RuntimeException(ERRORMissingCode(idt),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(CoefficientsWRT, 2) { // AMB
  intrusive_ptr<RINGELEM> f = runtimeEnv->evalArgAs<RINGELEM>(ARG(0));
  int which;
  intrusive_ptr<RightValue> idt = runtimeEnv->evalArgAsT1orT2<RINGELEM, LIST>(ARG(1), which);
  std::vector<CoeffPP> CoeffWRT;
  switch (which) {
  case 1: CoeffWRT=CoefficientsWRT(f->theRingElem, RefTo<RingElem>(idt)); break;
  case 2: {
//     vector<long> v1 = VectorLongDecr1(runtimeEnv->evalArgAsListOf<INT>(ARG(1)), ERR::BadIndetIndex, "ContentWRT");
//     CoeffWRT = CoefficientsWRT(f->theRingElem, v1);
    vector<RingElem> v1 = runtimeEnv->evalRVAsListOfRingElem(idt, ARG(1));
    CoeffWRT = CoefficientsWRT_forC5(f->theRingElem, v1);
    break;
  }
  default: throw RuntimeException(ERRORMissingCode(idt),invocationExpression);
  }
  // create return value: list of ...
  const SparsePolyRing P = owner(f->theRingElem);
  intrusive_ptr<LIST> returnValue(new LIST);
  for (const CoeffPP& cpp: CoeffWRT)
  {
    // create return value: ... record
    intrusive_ptr<RECORD> cpp5(new RECORD);
    cpp5->setFieldNoCheck("coeff", Value::from(cpp.myCoeff));
    cpp5->setFieldNoCheck("PP",    Value::from(monomial(P,cpp.myPP)));
    returnValue->addValue(cpp5);
  }
  return returnValue;
}
END_STD_BUILTIN_FUNCTION



//---- MAT -------------------------------------------------------

DECLARE_STD_BUILTIN_FUNCTION(SetEntry, 4) { // AMB
  intrusive_ptr<MAT> M = runtimeEnv->obtainUnshared<MAT>(ARG(0));
  intrusive_ptr<INT> I = runtimeEnv->evalArgAs<INT>(ARG(1));
  intrusive_ptr<INT> J = runtimeEnv->evalArgAs<INT>(ARG(2));
  int which;
  intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3<INT, RAT, RINGELEM>(ARG(3), which);
  RingElem x(RingOf(M->theMatrix));
  switch (which) {
  case 1: x = RefTo<BigInt>(v); break;
  case 2: x = RefTo<BigRat>(v); break;
  case 3: x = RefTo<RingElem>(v); break;
  default: throw RuntimeException(ERRORMissingCode(v), invocationExpression);
  }
  SetEntry_forC5(M->theMatrix, I->theBigInt, J->theBigInt, x);
  return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(SetRow, 3) { // AMB
  intrusive_ptr<MAT> M = runtimeEnv->obtainUnshared<MAT>(ARG(0));
  intrusive_ptr<INT> I = runtimeEnv->evalArgAs<INT>(ARG(1));
  vector<RingElem> v = runtimeEnv->evalArgAsListOfRingElem(RingOf(M->theMatrix), ARG(2));
  SetRow_forC5(M->theMatrix, I->theBigInt, v);
  return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(SetCol, 3) { // AMB
  intrusive_ptr<MAT> M = runtimeEnv->obtainUnshared<MAT>(ARG(0));
  intrusive_ptr<INT> I = runtimeEnv->evalArgAs<INT>(ARG(1));
  vector<RingElem> v = runtimeEnv->evalArgAsListOfRingElem(RingOf(M->theMatrix), ARG(2));
  SetCol_forC5(M->theMatrix, I->theBigInt, v);
  return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(SwapRows, 3) { // AMB
  intrusive_ptr<MAT> M = runtimeEnv->obtainUnshared<MAT>(ARG(0));
  intrusive_ptr<INT> I1 = runtimeEnv->evalArgAs<INT>(ARG(1));
  intrusive_ptr<INT> I2 = runtimeEnv->evalArgAs<INT>(ARG(2));
  SwapRows_forC5(M->theMatrix, I1->theBigInt, I2->theBigInt);
  return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(SwapCols, 3) { // AMB
  intrusive_ptr<MAT> M = runtimeEnv->obtainUnshared<MAT>(ARG(0));
  intrusive_ptr<INT> I1 = runtimeEnv->evalArgAs<INT>(ARG(1));
  intrusive_ptr<INT> I2 = runtimeEnv->evalArgAs<INT>(ARG(2));
  SwapCols_forC5(M->theMatrix, I1->theBigInt, I2->theBigInt);
  return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(AddRowMul, 4) { // AMB
  intrusive_ptr<MAT> M = runtimeEnv->obtainUnshared<MAT>(ARG(0));
  long i1 = runtimeEnv->evalArgAsLong(ARG(1));
  long i2 = runtimeEnv->evalArgAsLong(ARG(2));
  intrusive_ptr<RINGELEM> c = runtimeEnv->evalArgAs<RINGELEM>(ARG(3));
  AddRowMul(M->theMatrix, i1-1, i2-1, c->theRingElem);
  return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(AddColMul, 4) { // AMB
  intrusive_ptr<MAT> M = runtimeEnv->obtainUnshared<MAT>(ARG(0));
  long j1 = runtimeEnv->evalArgAsLong(ARG(1));
  long j2 = runtimeEnv->evalArgAsLong(ARG(2));
  intrusive_ptr<RINGELEM> c = runtimeEnv->evalArgAs<RINGELEM>(ARG(3));
  AddColMul(M->theMatrix, j1-1, j2-1, c->theRingElem);
  return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(IsZeroRow, 2) { // AMB
  intrusive_ptr<MAT> M = runtimeEnv->evalArgAs<MAT>(ARG(0));
  long N = runtimeEnv->evalArgAsLong(ARG(1));
  return Value::from(IsZeroRow(M->theMatrix, N-1)); ///check for underflow???
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(IsZeroCol, 2) { // AMB
  intrusive_ptr<MAT> M = runtimeEnv->evalArgAs<MAT>(ARG(0));
  long N = runtimeEnv->evalArgAsLong(ARG(1));
  return Value::from(IsZeroCol(M->theMatrix, N-1)); ///check for underflow???
}
END_STD_BUILTIN_FUNCTION



DECLARE_STD_BUILTIN_FUNCTION(HadamardBoundSq, 1) { // JAA
  intrusive_ptr<MAT> M = runtimeEnv->evalArgAs<MAT>(ARG(0));

  const HadamardRowCol bounds = HadamardBoundSq(M->theMatrix);
  // Create return value: ... record
  intrusive_ptr<RECORD> ans(new RECORD);
  ans->setFieldNoCheck("RowBoundSq", Value::from(bounds.myRowBoundSq));
  ans->setFieldNoCheck("ColBoundSq", Value::from(bounds.myColBoundSq));
  return ans;
}
END_STD_BUILTIN_FUNCTION




DECLARE_STD_BUILTIN_FUNCTION(ElimHomogMat, 2) { // AMB
  vector<long> v = VectorLongDecr1(runtimeEnv->evalArgAsListOf<INT>(ARG(0)),
                                   ERR::BadIndex,   "HomogElimMat");
  intrusive_ptr<MAT> M = runtimeEnv->evalArgAs<MAT>(ARG(1));
  return Value::from(ElimHomogMat(v, M->theMatrix));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(ElimMat, 2) {
  int which;
  intrusive_ptr<RightValue> arg0 = runtimeEnv->evalArgAsT1orT2orT3<INT, MAT, LIST>(ARG(0), which);
  switch (which) {
  case 1: 
  case 2: 
  throw RuntimeException("ElimMat(INT/MAT, ElimIndets) obsolete, use ElimMat(ElimIndets, INT/MAT) instead", invocationExpression);
  }
  vector<BigInt> ListOfBigInt;
  if (!IsVectorBigInt(ListOfBigInt, intrusive_ptr_cast<LIST>(arg0)))
  { RuntimeException("ElimMat(ElimIndets, INT/MAT) first arg must be list of integer indices", invocationExpression); }
  vector<long> elim = VectorLongDecr1(ListOfBigInt, ERR::BadIndex, "ElimMat");
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<INT, MAT>(ARG(1), which);
  switch (which) {
  case 1:
  { // long-winded... isn't there a better way?
    const long num_indets = ConvertTo<long>(RefTo<BigInt>(x));
    return Value::from(ElimMat(elim, num_indets));
  }
  case 2: return Value::from(ElimMat(elim, RefTo<matrix>(x)));
  default: throw RuntimeException(ERRORMissingCode(x),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION


// variable number of args (???)
DECLARE_STD_BUILTIN_FUNCTION(IsPositiveGrading,1) {
  intrusive_ptr<MAT> M = runtimeEnv->evalArgAs<MAT>(ARG(0));
  return Value::from(IsPositiveGrading(M->theMatrix));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(submat, 3) { // AMB
  intrusive_ptr<MAT> M = runtimeEnv->evalArgAs<MAT>(ARG(0));
  vector<long> rows = VectorLongDecr1(runtimeEnv->evalArgAsListOf<INT>(ARG(1)), ERR::BadIndex, "submat rows");
  vector<long> cols = VectorLongDecr1(runtimeEnv->evalArgAsListOf<INT>(ARG(2)), ERR::BadIndex, "submat cols");
  return new MAT(NewDenseMat(submat(M->theMatrix, rows, cols)));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(FirstCols, 2) { // AMB-2024-03
  intrusive_ptr<MAT> M = runtimeEnv->evalArgAs<MAT>(ARG(0));
  const long n = runtimeEnv->evalArgAsLong(ARG(1));
  return new MAT(NewDenseMat(FirstCols(M->theMatrix, n)));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(FirstRows, 2) { // AMB-2024-03
  intrusive_ptr<MAT> M = runtimeEnv->evalArgAs<MAT>(ARG(0));
  const long n = runtimeEnv->evalArgAsLong(ARG(1));
  return new MAT(NewDenseMat(FirstRows(M->theMatrix, n)));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(RandomSparseNonSing01Mat, 2) { // JAA
  intrusive_ptr<RING> R = runtimeEnv->evalArgAs<RING>(ARG(0));
  const long n = runtimeEnv->evalArgAsLong(ARG(1));
  return Value::from(RandomSparseNonSing01Mat(R->theRing, n));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(SylvesterMat, 3) {  // JAA
  intrusive_ptr<RINGELEM> f = runtimeEnv->evalArgAs<RINGELEM>(ARG(0));
  intrusive_ptr<RINGELEM> g = runtimeEnv->evalArgAs<RINGELEM>(ARG(1));
  intrusive_ptr<RINGELEM> x = runtimeEnv->evalArgAs<RINGELEM>(ARG(2));
  return Value::from(SylvesterMat(f->theRingElem, g->theRingElem, x->theRingElem));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(HilbertMat, 1) {  // JAA
  const long n = runtimeEnv->evalArgAsLong(ARG(0));
  return Value::from(HilbertMat(n));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(NR, 2) { // AMB
  intrusive_ptr<RINGELEM> f = runtimeEnv->evalArgAs<RINGELEM>(ARG(0));
  vector<RingElem> v = runtimeEnv->evalArgAsListOf<RINGELEM>(ARG(1));
  return Value::from(NR(f->theRingElem,v));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(TmpDivAlg, 2) { // AMB
  intrusive_ptr<RINGELEM> f = runtimeEnv->evalArgAs<RINGELEM>(ARG(0));
  vector<RingElem> v = runtimeEnv->evalArgAsListOf<RINGELEM>(ARG(1));
  return Value::from(TmpDivAlg(f->theRingElem,v));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(BinSequentialToric, 2) { // AMB
  intrusive_ptr<IDEAL> I = runtimeEnv->evalArgAs<IDEAL>(ARG(0));
  vector<BigInt> indices = runtimeEnv->evalArgAsListOf<INT>(ARG(1));
  return Value::from(SequentialToric_C(I->theIdeal, VectorLong(indices, "BinSequentialToric")));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(MatSequentialToric, 2) { // AMB
  intrusive_ptr<RING> R = runtimeEnv->evalArgAs<RING>(ARG(0));
  intrusive_ptr<MAT> M = runtimeEnv->evalArgAs<MAT>(ARG(1));
  return Value::from(SequentialToric_C(R->theRing, M->theMatrix));
}
END_STD_BUILTIN_FUNCTION



//---- MODULE -------------------------------------------------------

DECLARE_STD_BUILTIN_FUNCTION(ModuleElem, 2) { // AMB
  intrusive_ptr<MODULE> M = runtimeEnv->evalArgAs<MODULE>(ARG(0));
  vector<RingElem> v = runtimeEnv->evalArgAsListOfRingElem(RingOf(M->theModule), ARG(1));
  return Value::from(NewFreeModuleElem(M->theModule, v));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(NewFreeModule, 2) { // AMB
  intrusive_ptr<RING> R = runtimeEnv->evalArgAs<RING>(ARG(0));
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<INT, MAT>(ARG(1), which);
  switch (which) {
  case 1: return Value::from(NewFreeModule(R->theRing, ConvertTo<long>(RefTo<BigInt>(x))));
  case 2: return Value::from(NewFreeModule_forC5(R->theRing, RefTo<matrix>(x)));
  default: throw RuntimeException(ERRORMissingCode(x), invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(shifts, 1) { // AMB
  intrusive_ptr<MODULE> v = runtimeEnv->evalArgAs<MODULE>(ARG(0));
  return Value::from(shifts(v->theModule));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(NewFreeModuleForSyz, 1) { // AMB
  int which;
  intrusive_ptr<RightValue> v = runtimeEnv->evalArgAsT1orT2orT3<IDEAL, MODULE, LIST>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(NewFreeModuleForSyz(gens(RefTo<ideal>(v))));
  case 2: return Value::from(NewFreeModuleForSyz(gens(RefTo<module>(v))));
  case 3: {
    vector<RingElem> x = runtimeEnv->evalRVAsListOfRingElem(v, ARG(0));
    return Value::from(NewFreeModuleForSyz(x));
  }
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression); 
  }
}
END_STD_BUILTIN_FUNCTION



DECLARE_STD_BUILTIN_FUNCTION(GensAsRows, 1) { // AMB
  intrusive_ptr<MODULE> F = runtimeEnv->evalArgAs<MODULE>(ARG(0));
  return Value::from(GensAsRows(F->theModule));
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(GensAsCols, 1) { // AMB
  intrusive_ptr<MODULE> F = runtimeEnv->evalArgAs<MODULE>(ARG(0));
  return Value::from(GensAsCols(F->theModule));
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(NumCompts, 1) { // AMB
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<MODULE, MODULEELEM>(ARG(0), which);
  switch (which) {
  case 1: return Value::from(NumCompts_forC5(RefTo<module>(x)));
  case 2: return Value::from(NumCompts(RefTo<ModuleElem>(x)));
  default: throw RuntimeException(ERRORMissingCode(x),invocationExpression);
  }
}
END_STD_BUILTIN_FUNCTION


//---- RING -------------------------------------------------------

DECLARE_STD_BUILTIN_FUNCTION(NewRingTwinFloat, 1) { // AMB
  //  intrusive_ptr<INT> prec = runtimeEnv->evalArgAs<INT>(ARG(0));
  long d = runtimeEnv->evalArgAsLong(ARG(0));
  //  if (!IsConvertible(d, prec->theBigInt))
  //    throw RuntimeException("invalid precision", ARG(0).exp);
  return Value::from(NewRingTwinFloat(d));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(NewExtAlgebra, 2) { // AMB 2018-02
  intrusive_ptr<RING> R = runtimeEnv->evalArgAs<RING>(ARG(0));
  string IndetNames = runtimeEnv->evalArgAs<STRING>(ARG(1))->theString;
  return Value::from(NewExtAlgebra(R->theRing, symbols(IndetNames)));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(NewPolyRingWeights, 3) { // AMB 2024-03
  intrusive_ptr<RING> R = runtimeEnv->evalArgAs<RING>(ARG(0));
  vector<symbol> syms = runtimeEnv->evalArgAsListOfSymbols(ARG(1));
  intrusive_ptr<MAT> W = runtimeEnv->evalArgAs<MAT>(ARG(2));
  return new RING(NewPolyRingWeights(R->theRing, syms, W->theMatrix));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(NewQuotientRing, 2) { // AMB 2019-10
  intrusive_ptr<RING> R = runtimeEnv->evalArgAs<RING>(ARG(0));
  int which;
  intrusive_ptr<const RightValue> v = runtimeEnv->evalArgAsT1orT2<IDEAL,STRING>(ARG(1), which);
  switch (which) {
  case 1: return Value::from(NewQuotientRing(R->theRing, RefTo<ideal>(v)));
  case 2: return Value::from(NewQuotientRing(R->theRing, RefTo<string>(v)));
  default: throw RuntimeException("second arg must be IDEAL/STRING",ARG(1).exp);
  }
  return Value::from(RingZZ()); // just to keep the compiler quiet
}
END_STD_BUILTIN_FUNCTION


//---- IDEAL -------------------------------------------------------

DECLARE_STD_BUILTIN_FUNCTION(elim, 2) { // AMB
  intrusive_ptr<IDEAL> I = runtimeEnv->evalArgAs<IDEAL>(ARG(1));
  ideal J = I->theIdeal;
  int which;
  intrusive_ptr<RightValue> x = runtimeEnv->evalArgAsT1orT2<LIST, RINGELEM>(ARG(0), which);
  vector<RingElem> ElimIndets;
  switch (which) {
  case 1: ElimIndets = runtimeEnv->evalRVAsListOfRingElem(x, ARG(0));break;
  case 2: ElimIndets.push_back(RefTo<RingElem>(x));break;
  default: throw RuntimeException(ERRORMissingCode(x),invocationExpression);
  }
  //  MakeUnique(J)->myElim(ElimIndets); // I->myElim(v1) ~~~> const
  //  return Value::from(J);
  return Value::from(elim(J, ElimIndets));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(ComputeElimFirst, 2) { // AMB
  intrusive_ptr<IDEAL> I = runtimeEnv->evalArgAs<IDEAL>(ARG(1));
  intrusive_ptr<RINGELEM> ElimIndets = runtimeEnv->evalArgAs<RINGELEM>(ARG(0));
  return Value::from(ComputeElimFirst(gens(I->theIdeal), LPP(ElimIndets->theRingElem)));
}
END_STD_BUILTIN_FUNCTION


// now oneliner
// DECLARE_STD_BUILTIN_FUNCTION(saturate, 2) { // AMB
//   intrusive_ptr<IDEAL> I0 = runtimeEnv->evalArgAs<IDEAL>(ARG(0));
//   intrusive_ptr<IDEAL> I1 = runtimeEnv->evalArgAs<IDEAL>(ARG(1));
//   return Value::from(saturate(I0->theIdeal, I1->theIdeal));
// }
// END_STD_BUILTIN_FUNCTION




//---- IDEAL (points) -----------------------------------------------

// DECLARE_STD_BUILTIN_FUNCTION(IdealOfPoints, 2) { // JAA 2013-01-19
//  intrusive_ptr<RING> P = runtimeEnv->evalArgAs<RING>(ARG(0));
//  intrusive_ptr<MAT> pts = runtimeEnv->evalArgAs<MAT>(ARG(1));
//   return Value::from(IdealOfPoints(P->theRing, pts->theMatrix));
// }
// END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(ApproxPointsNBM, 3) { // AMB
  intrusive_ptr<RING> P = runtimeEnv->evalArgAs<RING>(ARG(0));
  intrusive_ptr<MAT> Pts = runtimeEnv->evalArgAs<MAT>(ARG(1));
  intrusive_ptr<MAT> Tols = runtimeEnv->evalArgAs<MAT>(ARG(2));
  vector<RingElem> QB;
  vector<RingElem> BB;
  vector<RingElem> AV;
  ApproxPointsNBM_forC5(QB, BB, AV, P->theRing, Pts->theMatrix, Tols->theMatrix);
  intrusive_ptr<RECORD> rec(new RECORD);
  rec->setFieldNoCheck("QuotientBasis", Value::from(QB));
  rec->setFieldNoCheck("AlmostVanishing", Value::from(AV));
  if (BB.empty())
    rec->setFieldNoCheck("StableBBasisFound", Value::from(false));
  else
  {
    rec->setFieldNoCheck("StableBBasisFound", Value::from(true));
    rec->setFieldNoCheck("BBasis", Value::from(BB));
  }
  return rec;
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(ApproxPointsSOI, 3) { // AMB
  intrusive_ptr<RING> P = runtimeEnv->evalArgAs<RING>(ARG(0));
  intrusive_ptr<MAT> Pts = runtimeEnv->evalArgAs<MAT>(ARG(1));
  intrusive_ptr<MAT> Tols = runtimeEnv->evalArgAs<MAT>(ARG(2));
  vector<RingElem> QB;
  vector<RingElem> BB;
  vector<RingElem> AV;
  ApproxPointsSOI_forC5(QB, BB, AV, P->theRing, Pts->theMatrix, Tols->theMatrix);
  intrusive_ptr<RECORD> rec(new RECORD);
  rec->setFieldNoCheck("QuotientBasis", Value::from(QB));
  rec->setFieldNoCheck("AlmostVanishing", Value::from(AV));
  if (BB.empty())
    rec->setFieldNoCheck("StableBBasisFound", Value::from(false));
  else
  {
    rec->setFieldNoCheck("StableBBasisFound", Value::from(true));
    rec->setFieldNoCheck("BBasis", Value::from(BB));
  }
  return rec;
}
END_STD_BUILTIN_FUNCTION


// DECLARE_STD_BUILTIN_FUNCTION(ClosePassingPoly, 3) { // AMB
//  intrusive_ptr<RING> P = runtimeEnv->evalArgAs<RING>(ARG(0));
//  intrusive_ptr<MAT> Pts = runtimeEnv->evalArgAs<MAT>(ARG(1));
//  intrusive_ptr<MAT> Tols = runtimeEnv->evalArgAs<MAT>(ARG(2));
//  intrusive_ptr<RAT> MaxTol = runtimeEnv->evalArgAs<RAT>(ARG(3));
//   ClosePassingPoly_forC5(P->theRing, Pts->theMatrix, Tols->theMatrix);
//   return Value::from(ClosePassingPoly_forC5(P->theRing, Pts->theMatrix, Tols->theMatrix);
// }
// END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(PreprocessPts, 2) { // JAA
  intrusive_ptr<MAT> OrigPts = runtimeEnv->evalArgAs<MAT>(ARG(0));
  intrusive_ptr<MAT> epsilon = runtimeEnv->evalArgAs<MAT>(ARG(1));
        const ring R = RingOf(OrigPts->theMatrix);
        vector< vector<RingElem> > NewPts;
        vector<long> weights;
        PreprocessPts_forC5("auto", NewPts, weights, OrigPts->theMatrix, epsilon->theMatrix);
//  NBM_forC5(QB, BB, AV, P->theRing, Pts->theMatrix, Tols->theMatrix);
  intrusive_ptr<RECORD> rec(new RECORD);
  rec->setFieldNoCheck("NewPoints", Value::from(NewDenseMat(R,NewPts)));
  rec->setFieldNoCheck("weights", Value::from(weights));
  return rec;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(PreprocessPtsGrid, 2) { // JAA
  intrusive_ptr<MAT> OrigPts = runtimeEnv->evalArgAs<MAT>(ARG(0));
  intrusive_ptr<MAT> epsilon = runtimeEnv->evalArgAs<MAT>(ARG(1));
        const ring R = RingOf(OrigPts->theMatrix);
        vector< vector<RingElem> > NewPts;
        vector<long> weights;
        PreprocessPts_forC5("grid", NewPts, weights, OrigPts->theMatrix, epsilon->theMatrix);
//  NBM_forC5(QB, BB, AV, P->theRing, Pts->theMatrix, Tols->theMatrix);
  intrusive_ptr<RECORD> rec(new RECORD);
  rec->setFieldNoCheck("NewPoints", Value::from(NewDenseMat(R,NewPts)));
  rec->setFieldNoCheck("weights", Value::from(weights));
  return rec;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(PreprocessPtsAggr, 2) { // JAA
  intrusive_ptr<MAT> OrigPts = runtimeEnv->evalArgAs<MAT>(ARG(0));
  intrusive_ptr<MAT> epsilon = runtimeEnv->evalArgAs<MAT>(ARG(1));
        const ring R = RingOf(OrigPts->theMatrix);
        vector< vector<RingElem> > NewPts;
        vector<long> weights;
        PreprocessPts_forC5("aggr", NewPts, weights, OrigPts->theMatrix, epsilon->theMatrix);
//  NBM_forC5(QB, BB, AV, P->theRing, Pts->theMatrix, Tols->theMatrix);
  intrusive_ptr<RECORD> rec(new RECORD);
  rec->setFieldNoCheck("NewPoints", Value::from(NewDenseMat(R,NewPts)));
  rec->setFieldNoCheck("weights", Value::from(weights));
  return rec;
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(PreprocessPtsSubdiv, 2) { // JAA
  intrusive_ptr<MAT> OrigPts = runtimeEnv->evalArgAs<MAT>(ARG(0));
  intrusive_ptr<MAT> epsilon = runtimeEnv->evalArgAs<MAT>(ARG(1));
        const ring R = RingOf(OrigPts->theMatrix);
        vector< vector<RingElem> > NewPts;
        vector<long> weights;
        PreprocessPts_forC5("subdiv", NewPts, weights, OrigPts->theMatrix, epsilon->theMatrix);
//  NBM_forC5(QB, BB, AV, P->theRing, Pts->theMatrix, Tols->theMatrix);
  intrusive_ptr<RECORD> rec(new RECORD);
  rec->setFieldNoCheck("NewPoints", Value::from(NewDenseMat(R,NewPts)));
  rec->setFieldNoCheck("weights", Value::from(weights));
  return rec;
}
END_STD_BUILTIN_FUNCTION

//---- IDEAL (implicit -- temporary)

DECLARE_STD_BUILTIN_FUNCTION(ImplicitDirect, 1) {  // JAA
  vector<RingElem> ParamDescr = runtimeEnv->evalArgAsListOfRingElem(ARG(0));
        RingElem ans = ImplicitDirect(ParamDescr);
        return Value::from(ans);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(ImplicitDirectLPP, 1) {  // JAA
  vector<RingElem> ParamDescr = runtimeEnv->evalArgAsListOfRingElem(ARG(0));
        RingElem ans = ImplicitDirectLPP(ParamDescr);
        return Value::from(ans);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(ImplicitDirectLPP2, 1) {  // JAA
  vector<RingElem> ParamDescr = runtimeEnv->evalArgAsListOfRingElem(ARG(0));
        RingElem ans = ImplicitDirectLPP2(ParamDescr);
        return Value::from(ans);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(ImplicitDirectWithCond, 2) {  // JAA
  vector<RingElem> ParamDescr = runtimeEnv->evalArgAsListOfRingElem(ARG(0));
  vector<RingElem> relations = runtimeEnv->evalArgAsListOfRingElem(ARG(1));
        RingElem ans = ImplicitDirectWithCond(ParamDescr, relations);
        return Value::from(ans);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(ImplicitDirectWithCondLPP, 2) {  // JAA
  vector<RingElem> ParamDescr = runtimeEnv->evalArgAsListOfRingElem(ARG(0));
  vector<RingElem> relations = runtimeEnv->evalArgAsListOfRingElem(ARG(1));
        RingElem ans = ImplicitDirectWithCondLPP(ParamDescr, relations);
        return Value::from(ans);
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(ImplicitByPoints, 1) {  // JAA
  vector<RingElem> ParamDescr = runtimeEnv->evalArgAsListOfRingElem(ARG(0));
        RingElem ans = ImplicitByPoints(ParamDescr);
        return Value::from(ans);
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(ImplicitByPoints2, 1) {  // JAA
  vector<RingElem> ParamDescr = runtimeEnv->evalArgAsListOfRingElem(ARG(0));
        RingElem ans = ImplicitByPoints2(ParamDescr);
        return Value::from(ans);
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(ImplicitByPoints3, 1) {  // JAA
  vector<RingElem> ParamDescr = runtimeEnv->evalArgAsListOfRingElem(ARG(0));
        RingElem ans = ImplicitByPoints3(ParamDescr);
        return Value::from(ans);
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(SliceCore, 3) {  // AMB
  vector<RingElem> ParamDescr = runtimeEnv->evalArgAsListOfRingElem(ARG(0));
  long RecDepth =runtimeEnv->evalArgAsLong(ARG(1));
  string FinalCalls = runtimeEnv->evalArgAs<STRING>(ARG(2))->theString;
  return Value::from(SliceCore(ParamDescr, RecDepth, FinalCalls));
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(SliceCoreQQ, 3) {  // AMB
  vector<RingElem> ParamDescr = runtimeEnv->evalArgAsListOfRingElem(ARG(0));
  long RecDepth =runtimeEnv->evalArgAsLong(ARG(1));
  string FinalCalls = runtimeEnv->evalArgAs<STRING>(ARG(2))->theString;
  return Value::from(SliceCoreQQ(ParamDescr, RecDepth, FinalCalls));
}
END_STD_BUILTIN_FUNCTION


//------ POLY ---------------------------------------------------------

DECLARE_STD_BUILTIN_FUNCTION(ChebyshevPoly, 2) {
  const BigInt N = runtimeEnv->evalArgAs<INT>(ARG(0))->theBigInt;
  const RingElem x = runtimeEnv->evalArgAs<RINGELEM>(ARG(1))->theRingElem;
  const char* const FnName = "ChebyshevPoly";
  long n;
  if (!IsConvertible(n, N))
    CoCoA_THROW_ERROR(ERR::ArgTooBig, FnName);
  if (n < 0)
    CoCoA_THROW_ERROR(ERR::ReqNonNegative, FnName);
  return Value::from(ChebyshevPoly(n,x));
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(ChebyshevPoly2, 2) {
  const BigInt N = runtimeEnv->evalArgAs<INT>(ARG(0))->theBigInt;
  const RingElem x = runtimeEnv->evalArgAs<RINGELEM>(ARG(1))->theRingElem;
  const char* const FnName = "ChebyshevPoly2";
  long n;
  if (!IsConvertible(n, N))
    CoCoA_THROW_ERROR(ERR::ArgTooBig, FnName);
  if (n < 0)
    CoCoA_THROW_ERROR(ERR::ReqNonNegative, FnName);
  return Value::from(ChebyshevPoly2(n,x));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(HermitePoly, 2) {
  const BigInt N = runtimeEnv->evalArgAs<INT>(ARG(0))->theBigInt;
  const RingElem x = runtimeEnv->evalArgAs<RINGELEM>(ARG(1))->theRingElem;
  const char* const FnName = "HermitePoly";
  long n;
  if (!IsConvertible(n, N))
    CoCoA_THROW_ERROR(ERR::ArgTooBig, FnName);
  if (n < 0)
    CoCoA_THROW_ERROR(ERR::ReqNonNegative, FnName);
  return Value::from(HermitePoly(n,x));
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(HermitePoly2, 2) {
  const BigInt N = runtimeEnv->evalArgAs<INT>(ARG(0))->theBigInt;
  const RingElem x = runtimeEnv->evalArgAs<RINGELEM>(ARG(1))->theRingElem;
  const char* const FnName = "HermitePoly2";
  long n;
  if (!IsConvertible(n, N))
    CoCoA_THROW_ERROR(ERR::ArgTooBig, FnName);
  if (n < 0)
    CoCoA_THROW_ERROR(ERR::ReqNonNegative, FnName);
  return Value::from(HermitePoly2(n,x));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(LaguerrePoly, 2) {
  const BigInt N = runtimeEnv->evalArgAs<INT>(ARG(0))->theBigInt;
  const RingElem x = runtimeEnv->evalArgAs<RINGELEM>(ARG(1))->theRingElem;
  const char* const FnName = "LaguerrePoly";
  long n;
  if (!IsConvertible(n, N))
    CoCoA_THROW_ERROR(ERR::ArgTooBig, FnName);
  if (n < 0)
    CoCoA_THROW_ERROR(ERR::ReqNonNegative, FnName);
  return Value::from(LaguerrePoly(n,x));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(BesselPoly, 2) {
  const BigInt N = runtimeEnv->evalArgAs<INT>(ARG(0))->theBigInt;
  const RingElem x = runtimeEnv->evalArgAs<RINGELEM>(ARG(1))->theRingElem;
  const char* const FnName = "BesselPoly";
  long n;
  if (!IsConvertible(n, N))
    CoCoA_THROW_ERROR(ERR::ArgTooBig, FnName);
  if (n < 0)
    CoCoA_THROW_ERROR(ERR::ReqNonNegative, FnName);
  return Value::from(BesselPoly(n,x));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(DicksonPoly, 3) {
  const RingElem x = runtimeEnv->evalArgAs<RINGELEM>(ARG(0))->theRingElem;
  const BigInt N = runtimeEnv->evalArgAs<INT>(ARG(1))->theBigInt;
  int WhichType;
  intrusive_ptr<const RightValue> v = runtimeEnv->evalArgAsT1orT2orT3<RINGELEM,INT,RAT>(ARG(2), WhichType);
  RingElem alpha;
  switch (WhichType)
  {
  case 1: alpha = RefTo<RingElem>(v); break;
  case 2: alpha = RingElem(owner(x), RefTo<BigInt>(v)); break;
  case 3: alpha = RingElem(owner(x), RefTo<BigRat>(v)); break;
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression);
  }
  const char* const FnName = "DicksonPoly";
  long n;
  if (!IsConvertible(n, N) || n < 0)
    CoCoA_THROW_ERROR(ERR::ArgTooBig, FnName);
  if (n < 0)
    CoCoA_THROW_ERROR(ERR::ReqNonNegative, FnName);
  return Value::from(DicksonPoly(x,n,alpha));
}
END_STD_BUILTIN_FUNCTION

DECLARE_STD_BUILTIN_FUNCTION(DicksonPoly2, 3) {
  const RingElem x = runtimeEnv->evalArgAs<RINGELEM>(ARG(0))->theRingElem;
  const BigInt N = runtimeEnv->evalArgAs<INT>(ARG(1))->theBigInt;
  int WhichType;
  intrusive_ptr<const RightValue> v = runtimeEnv->evalArgAsT1orT2orT3<RINGELEM,INT,RAT>(ARG(2), WhichType);
  RingElem alpha;
  switch (WhichType)
  {
  case 1: alpha = RefTo<RingElem>(v); break;
  case 2: alpha = RingElem(owner(x), RefTo<BigInt>(v)); break;
  case 3: alpha = RingElem(owner(x), RefTo<BigRat>(v)); break;
  default: throw RuntimeException(ERRORMissingCode(v),invocationExpression);
  }
  const char* const FnName = "DicksonPoly2";
  long n;
  if (!IsConvertible(n, N) || n < 0)
    CoCoA_THROW_ERROR(ERR::ArgTooBig, FnName);
  if (n < 0)
    CoCoA_THROW_ERROR(ERR::ReqNonNegative, FnName);
  return Value::from(DicksonPoly2(x,n,alpha));
}
END_STD_BUILTIN_FUNCTION


DECLARE_STD_BUILTIN_FUNCTION(SwinnertonDyerPoly, 2) { // JAA
  const RingElem x = runtimeEnv->evalArgAs<RINGELEM>(ARG(0))->theRingElem;
  intrusive_ptr<LIST> L = runtimeEnv->evalArgAs<LIST>(ARG(1));
  vector<BigRat> roots;
  if (!IsVectorBigRat(roots, L)) CoCoA_THROW_ERROR("Expected list of INT/RAT", "SwinnertonDyerPoly");
  return Value::from(SwinnertonDyerPoly(x,roots));
}
END_STD_BUILTIN_FUNCTION


// -- QuasiPoly -------------------------------------------------------
DECLARE_STD_BUILTIN_FUNCTION(EvalQuasiPoly, 2) {
  const vector<RingElem> QPoly = runtimeEnv->evalArgAsListOfRingElem(ARG(0));
  const BigInt N = runtimeEnv->evalArgAs<INT>(ARG(1))->theBigInt;
  const char* const FnName = "EvalQuasiPoly";
  if (QPoly.empty())
    CoCoA_THROW_ERROR(ERR::ReqNonEmpty, FnName);
  if (!HasUniqueOwner(QPoly))
    CoCoA_THROW_ERROR(ERR::MixedRings, FnName);
  const long i = ConvertTo<long>(LeastNNegRemainder(N, len(QPoly)));
  const RingHom EvalAtN = EvalHom(owner(QPoly[i]), N);
  return Value::from(EvalAtN(QPoly[i]));
}
END_STD_BUILTIN_FUNCTION


//-- UTILITIES ---------------------------------------------------------

DECLARE_STD_BUILTIN_FUNCTION(SetVerbosityLevel, 1) { // AMB
  long n =runtimeEnv->evalArgAsLong(ARG(0));
  SetVerbosityLevel(n);
  return VoidValue::theInstance;
}
END_STD_BUILTIN_FUNCTION


} // namespace InterpreterNS
} // namespace CoCoA

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/CoCoA-5/BuiltInFunctions-CoCoALib.C,v 1.111 2024/07/15 17:18:49 bigatti Exp $
// $Log: BuiltInFunctions-CoCoALib.C,v $
// Revision 1.111  2024/07/15 17:18:49  bigatti
// Summary: Error codes:  added CoCoA_THROW_ERROR1(ErrID);
// Removed BadRowIndex, BadColIndex;  Renamed  ERR::Empty --> ERR::ReqNonEmpty
//
// Revision 1.110  2024/07/03 09:44:11  bigatti
// Summary: error codes (redmine 92):
// removed BadIntedIndex, BadDegIndex, BadComptIndex
//
// Revision 1.109  2024/07/02 20:42:08  abbott
// Summary: Renaming of errors (redmine 308)
//
// Revision 1.108  2024/05/28 13:08:45  bigatti
// Summary: added GBasisTrunc
//
// Revision 1.107  2024/03/26 12:35:28  bigatti
// Summary: moved syz from BuiltInFunctions-CoCoALib to BuiltInFunctionsVarArgs-CoCoALib
//
// Revision 1.106  2024/03/21 12:29:48  bigatti
// Summary: saturate moved to OneLiners
//
// Revision 1.105  2024/03/20 15:37:13  bigatti
// Summary: added FirstCols/FirstRows
//
// Revision 1.104  2024/03/15 19:06:00  abbott
// Summary: Corrected body of saturate
//
// Revision 1.103  2024/03/15 16:01:10  bigatti
// Summary: added NewPolyRingWeights;  now using new function elim
//
// Revision 1.102  2024/02/05 13:32:45  bigatti
// Summary: now radical(ideal)  is in C++
//
// Revision 1.101  2023/12/01 20:27:59  abbott
// Summary: Changed names of fields in record produced by HadamardBoundSq (redmine 1771)
//
// Revision 1.100  2023/07/04 09:21:11  abbott
// Summary: Changed TimeOut to Timeout
//
// Revision 1.99  2023/06/06 19:20:59  abbott
// Summary: Now have CyclotomicTest and CyclotomicIndex
//
// Revision 1.98  2023/05/11 19:29:44  abbott
// Summary: Added temporary fn RECOG
//
// Revision 1.97  2023/03/14 22:35:51  abbott
// Summary: Renamed GBasis_TimeOut
//
// Revision 1.96  2023/03/06 22:10:56  abbott
// Summary: Updated IsCoprime (redmine 1710); added new fn LRSDegerancyOrders
//
// Revision 1.95  2022/06/08 15:24:33  abbott
// Summary: Added new fn RandomNBitPrime
//
// Revision 1.94  2022/03/17 14:51:37  abbott
// Summary: Removed several superfluous include directives
//
// Revision 1.93  2022/03/09 07:58:25  bigatti
// Summary: added TmpDivAlg
//
// Revision 1.92  2022/03/07 14:32:21  abbott
// Summary: Removed some cruft from "syz"
//
// Revision 1.91  2022/03/07 09:47:10  bigatti
// Summary: added func syz(const FreeModule& F, const std::vector<RingElem>& g);
//
// Revision 1.90  2022/03/04 11:26:39  bigatti
// Summary: homogenized vertical spacing; auto instead of intrusive_ptr
//
// Revision 1.89  2022/02/22 20:39:26  abbott
// Summary: Updated copyright message (redmine 1555)
//
// Revision 1.88  2022/02/17 20:20:09  abbott
// Summary: Commented out "apply" (redmine 1598)
//
// Revision 1.87  2021/10/14 12:11:49  abbott
// Summary: Added RandomSmallPrime
//
// Revision 1.86  2021/10/04 09:04:10  abbott
// Summary: Replaced calls to apply by direct applic of ringhom (redmine 1467, 1598)
//
// Revision 1.85  2021/07/30 15:19:47  bigatti
// Summary: using evalRVAsListOfRingElem Redmine #946
//
// Revision 1.84  2021/07/20 11:29:23  abbott
// Summary: Replaced apply(phi,x) by phi(x) where x is RingElem.
//
// Revision 1.83  2021/07/19 11:29:27  abbott
// Summary: Avoid double evaluation (redmine 946)
//
// Revision 1.82  2021/04/26 14:12:43  abbott
// Summary: ExtGcd now gives helpful error rather than returning false on bad input.
//
// Revision 1.81  2021/02/18 16:43:27  abbott
// Summary: Revised impl of CpuTimeLimit (redmine 1558)
//
// Revision 1.80  2021/01/15 16:59:34  abbott
// Summary: Redmine 1563: new ctor for BigRat direct from integer
//
// Revision 1.79  2021/01/08 17:42:41  abbott
// Summary: Added "high variability" for CpuTimeLimit in GBasisTimeout (redmine 1376)
//
// Revision 1.78  2020/12/04 10:15:22  abbott
// Summary: Added space after keyword "for"
//
// Revision 1.77  2020/11/23 09:10:28  abbott
// Summary: Added new fn RingElemList; BOOST_FOREACH changed to for (redmine 1520)
//
// Revision 1.76  2020/10/14 20:06:03  abbott
// Summary: Added IsDivisible_AllowFields
//
// Revision 1.75  2020/09/22 18:22:07  abbott
// Summary: Added buultin fn interreduced
//
// Revision 1.74  2020/06/17 15:49:30  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.73  2020/02/28 08:57:43  abbott
// Summary: Added PROTOTYPE ExtGcd
//
// Revision 1.72  2020/02/13 15:50:56  bigatti
// -- added IndetsProd for LIST
//
// Revision 1.71  2020/02/13 14:03:50  bigatti
// -- now IsDivisible works for INT too
//
// Revision 1.70  2020/02/03 17:06:58  abbott
// Summary: Added new fn fold
//
// Revision 1.69  2019/10/17 20:20:42  abbott
// Summary: Disabled len(RINGELEM)
//
// Revision 1.68  2019/10/11 12:55:36  abbott
// Summary: Added SylvesterMat
//
// Revision 1.67  2019/10/10 06:59:25  bigatti
// -- added PolyAlgebraHom PolyRingHom with string
//
// Revision 1.66  2019/10/08 13:01:22  bigatti
// -- added NewQuotientRing with string
//
// Revision 1.65  2019/09/25 16:08:09  bigatti
// -- functions with variable number of args moved to BuiltInFunctionsVarArgs-CoCoALib
//
// Revision 1.64  2019/09/16 17:46:21  abbott
// Summary: Added IsVectorBigInt; better args checks for ChebyshevPoly, HermitePoly, LaguerrePoly
//
// Revision 1.63  2019/04/01 13:33:36  bigatti
// -- some includes moved here from BuiltInFunctions.H
//
// Revision 1.62  2019/03/27 14:48:13  bigatti
// -- added NumGens (also for MODULE)
//
// Revision 1.61  2019/03/27 14:22:08  bigatti
// (abbott) renamed GCDFreeBasis --> CoprimeFactorBasis
//
// Revision 1.60  2019/03/04 13:10:05  abbott
// Summary: Added DicksonPoly
//
// Revision 1.59  2018/08/28 12:38:46  abbott
// Summary: Added RandomPermutation
//
// Revision 1.58  2018/06/27 12:15:42  abbott
// Summary: Removed cruft
//
// Revision 1.57  2018/06/27 08:51:27  abbott
// Summary: Changed to work with new CpuTimeLimit
//
// Revision 1.56  2018/06/25 12:33:14  abbott
// Summary: Added GCDFreeBasis
//
// Revision 1.55  2018/06/13 14:54:04  abbott
// Summary: Cleaned impl of IsHomog
//
// Revision 1.54  2018/03/20 15:54:21  bigatti
// -- fixed arity RandomLinearForm
//
// Revision 1.53  2018/03/20 13:49:47  bigatti
// -- added CommonDenom
// -- added RandomLinearForm
// -- changed unused default behaviour after which
//
// Revision 1.52  2018/03/13 18:03:52  bigatti
// -- MinPolyModular: now using VerificationLevel class
//
// Revision 1.51  2018/03/13 17:43:38  bigatti
// -- now MinPolyQuot takes a verification level
//    instead of being called MinPolyQuotHeuristic
//
// Revision 1.50  2018/02/22 14:59:12  bigatti
// -- untabified
//
// Revision 1.49  2018/02/22 14:34:59  abbott
// Summary: Cleaned up impl for jacobian
//
// Revision 1.48  2018/02/19 10:18:21  abbott
// Summary: Updated RatReconstructByContFrac (after changing ctor)
//
// Revision 1.47  2018/02/02 02:06:49  bigatti
// -- added NewExtAlgebra
//
// Revision 1.46  2017/12/21 10:53:13  bigatti
// -- fixed bug IsInImage
//
// Revision 1.45  2017/11/20 15:51:56  bigatti
// -- removed minimalized (obsolescent)
//
// Revision 1.44  2017/11/08 14:07:32  abbott
// Summary: Added new fns HilbertMat and RandomSparseNonSing01Mat
//
// Revision 1.43  2017/10/17 10:00:29  abbott
// Summary: Added new fns: ChebyshevPoly(good impl), HermitePoly, LaguerrePoly
//
// Revision 1.42  2017/09/14 15:57:49  abbott
// Summary: Added RootBound
//
// Revision 1.41  2017/07/24 12:07:08  abbott
// Summary: Corrected silly bug
//
// Revision 1.40  2017/07/23 15:31:34  abbott
// Summary: Added GBasisTimeout (just for ideals)
//
// Revision 1.39  2017/07/19 16:39:54  abbott
// Summary: IsInRadical & MinPowerInIeal are now built-in fns
//
// Revision 1.38  2017/05/22 16:16:27  abbott
// Summary: Added reseed
//
// Revision 1.37  2017/05/16 16:24:08  bigatti
// -- fixed double evaluation of some arguments (/redmine/issues/946)
//
// Revision 1.36  2017/05/02 12:06:22  bigatti
// -- just a comment
//
// Revision 1.35  2017/04/27 14:55:00  bigatti
// -- ReadExpr --> RingElem
//
// Revision 1.34  2017/04/26 15:57:22  bigatti
// -- now IdealOfGBasis, IdealOfMinGens are in cocoalib (moved from Supplement)
//
// Revision 1.33  2017/04/26 09:05:22  bigatti
// -- updated RingElem (same as ReadExpr for STRING)
// -- updated author
//
// Revision 1.32  2017/04/18 09:22:32  bigatti
// -- fixed slug in NewPolyRing
//
// Revision 1.31  2017/03/20 08:38:45  bigatti
// -- TmpNBM --> ApproxPointsNBM
// -- first prototype for SOI (not working yet)
//
// Revision 1.30  2017/03/13 17:23:41  bigatti
// -- now using IdealOfMinGens_forC5
//
// Revision 1.29  2017/03/02 10:04:22  bigatti
// -- modified interface for VerbosityLevel
//
// Revision 1.28  2016/10/27 14:07:58  abbott
// Summary: Added RandomUnimodularMat
//
// Revision 1.27  2016/10/25 20:54:50  abbott
// Summary: Added new fn IsSqFree
//
// Revision 1.26  2016/10/20 18:05:55  abbott
// Summary: Exposed "radical" under temporary name "rad"
//
// Revision 1.25  2016/09/22 15:33:37  bigatti
// -- renamed HomogElimMat into ElimHomogMat
// -- improved readability for ElimHomogMat/ElimMat (removed auxiliary functions)
//
// Revision 1.24  2016/09/22 14:38:15  bigatti
// -- removed HomogElimMat (now in obsolescent.cpkg5)
//
// Revision 1.23  2016/09/22 14:14:35  bigatti
// -- modified ElimMat and ElimHomogMat
//
// Revision 1.22  2016/09/21 16:34:59  bigatti
// -- changed year
// -- changed implementation for (homog)ElimMat using VectorLongDecr1
//    (and removed from CoCoALibSupplement)
//
// Revision 1.21  2016/08/02 09:54:17  bigatti
// -- commented out  NumDigits  (and changed manual suggesting FloorLog10)
//
// Revision 1.20  2016/06/27 14:49:21  bigatti
// -- now FrobeniusMat may take two args
//
// Revision 1.19  2016/06/24 14:27:41  bigatti
// -- renamed CRT_poly --> CRTPoly
//
// Revision 1.18  2016/06/10 15:55:52  bigatti
// now FrobeniusMat in BuiltInOneLiners-CoCoALib
//
// Revision 1.17  2016/04/14 11:33:13  bigatti
// -- added FrobeniusMat
//
// Revision 1.16  2016/04/14 08:07:25  bigatti
// -- in CoefficientsWRT now using monomial(P, pp)  without coeff
//
// Revision 1.15  2016/03/25 20:15:23  abbott
// Summary: New impls for MantissaAndExponent (2 & 10); removed some cruft.
//
// Revision 1.14  2016/02/17 20:02:55  abbott
// Summary: Corrected fn name in error mesg (in EvalQuasiPoly)
//
// Revision 1.13  2016/02/17 10:16:09  abbott
// Summary: Renamed EvalQuasiPoly and moved here (from BuiltInFns-Normaliz)
//
// Revision 1.12  2016/02/09 15:04:59  bigatti
// -- added AddRowMul, AddColMul
//
// Revision 1.11  2016/02/01 13:17:37  abbott
// Summary: Added NewRingFq fns
//
// Revision 1.10  2016/01/27 13:35:52  bigatti
// -- added NewRingFqVec, NewRingFqLog
//
// Revision 1.9  2016/01/26 13:56:29  bigatti
// -- just some spaces
//
// Revision 1.8  2015/12/09 10:19:52  abbott
// Summary: Changed name CompleteToOrd to MakeTermOrd
//
// Revision 1.7  2015/12/01 16:54:00  abbott
// Summary: Hasty hack just to get cocoa5 to compile: modified CompleteToOrd, IsPositive, removed ExtractOrdMat
//
// Revision 1.6  2015/12/01 13:43:15  abbott
// Summary: Removed AssignZero (for matrix)
//
// Revision 1.5  2015/11/30 21:53:56  abbott
// Summary: Major update to matrices for orderings (not yet complete, some tests fail)
//
// Revision 1.4  2015/11/23 18:23:05  abbott
// Summary: Renamed ILogBase -> FloorLogBase; added FloorLog2, FloorLog10
//
// Revision 1.3  2015/11/21 19:19:32  abbott
// Summary: Added SimplestBinaryRatBetween
//
// Revision 1.2  2015/10/08 13:10:04  bigatti
// -- added revision log at the end
//
