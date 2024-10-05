//   Copyright (c)  2005-2018,2022  John Abbott and Anna M. Bigatti
//   Authors:  2005-2018  John Abbott and Anna M. Bigatti

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


// Source code for abstract class SparsePolyRing and friends

#include "CoCoA/SparsePolyOps-RingElem.H"

#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/DUPFp.H"
#include "CoCoA/MatrixOps.H"  // for LinSolve
#include "CoCoA/MatrixView.H" // for ZeroMat
#include "CoCoA/NumTheory-RatReconstruct.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/QuotientRing.H" // for IsQuotientRing
#include "CoCoA/RingDistrMPolyClean.H" // for NewPolyRing_DMP
#include "CoCoA/RingFp.H" // for IsRingFp
#include "CoCoA/RingQQ.H" // for IsQQ
#include "CoCoA/RingTwinFloat.H" // for IsRingTwinFloat
#include "CoCoA/RingZZ.H" // for IsZZ
#include "CoCoA/SmallFpImpl.H"
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/SparsePolyOps-ideal.H"
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/factor.H"    // for myGcd
#include "CoCoA/geobucket.H" // for myMul
#include "CoCoA/ideal.H"     // for myGcd
#include "CoCoA/interrupt.H"
#include "CoCoA/matrix.H"    // for OrdMat, myDivMod
#include "CoCoA/module.H"    // for myGcd
#include "CoCoA/random.H"    // for RandomLong
#include "CoCoA/submodule.H" // for myGcd
#include "CoCoA/symbol.H"


#include <algorithm>
using std::min;     // for PPContent
using std::max;     // for MaxExponent, StdDeg
#include <iostream>
// using std::ostream in SparsePolyRingBase::myOutput
#include <utility>
using std::pair;
//#include <vector>
using std::vector;

namespace CoCoA
{


  //---- Functions for creating/building polynomials

  RingElem monomial(const SparsePolyRing& P, ConstRefRingElem c, ConstRefPPMonoidElem pp)
  {
    if (owner(c) != CoeffRing(P)) CoCoA_THROW_ERROR(ERR::MixedCoeffRings, "monomial(P,c,pp)");
    if (owner(pp) != PPM(P)) CoCoA_THROW_ERROR(ERR::MixedPPMs, "monomial(P,c,pp)");
    if (IsZero(c)) return zero(P);
    return P->myMonomial(raw(c), raw(pp));
  }

  RingElem monomial(const SparsePolyRing& P, ConstRefPPMonoidElem pp)
  { return monomial(P, one(CoeffRing(P)), pp); }

  RingElem monomial(const SparsePolyRing& P, const MachineInt& n, ConstRefPPMonoidElem pp)
  { return monomial(P, RingElem(CoeffRing(P), n), pp); }

  RingElem monomial(const SparsePolyRing& P, const BigInt& N, ConstRefPPMonoidElem pp)
  { return monomial(P, RingElem(CoeffRing(P), N), pp); }

  RingElem monomial(const SparsePolyRing& P, const BigRat& Q, ConstRefPPMonoidElem pp)
  { return monomial(P, RingElem(CoeffRing(P), Q), pp); }

  RingElem monomial(const SparsePolyRing& P, ConstRefRingElem c, const std::vector<long>& expv)
  { return monomial(P, c, PPMonoidElem(PPM(P), expv)); }

  RingElem monomial(const SparsePolyRing& P, const std::vector<long>& expv)
  { return monomial(P, PPMonoidElem(PPM(P), expv)); }

  RingElem monomial(const SparsePolyRing& P, const MachineInt& n, const std::vector<long>& expv)
  { return monomial(P, RingElem(CoeffRing(P), n), PPMonoidElem(PPM(P), expv)); }

  RingElem monomial(const SparsePolyRing& P, const BigInt& N, const std::vector<long>& expv)
  { return monomial(P, RingElem(CoeffRing(P), N), PPMonoidElem(PPM(P), expv)); }

  RingElem monomial(const SparsePolyRing& P, const BigRat& Q, const std::vector<long>& expv)
  { return monomial(P, RingElem(CoeffRing(P), Q), PPMonoidElem(PPM(P), expv)); }


  RingElem SparsePolyRingBase::mySymbolValue(const symbol& s) const
  {
    std::vector<symbol> syms = symbols(myPPM());    // makes a copy of vec
    syms.push_back(s);
    if (!AreDistinct(syms))   // should be IsInSymbols(...)
      return myMonomial(raw(one(myCoeffRing())), raw(myPPM()->mySymbolValue(s)));
    if (AreArityConsistent(syms))
      return myCoeffEmbeddingHomCtor()(myCoeffRing()->mySymbolValue(s));
    CoCoA_THROW_ERROR(ERR::BadIndetNames, "SparsePolyRingBase::mySymbolValue");
    return myZero(); // just to keep the compiler quiet
  }


  //---- RandomLinearForm ---------------------
  namespace  { // anonymous

    RingElem RandomLinearForm(const ring& P, long lo, long hi)
    {
      if (!IsSparsePolyRing(P))
        CoCoA_THROW_ERROR(ERR::ReqSparsePolyRing,"RandomLinearForm");
      if (lo >= hi)
        CoCoA_THROW_ERROR("Bad range lo-hi","RandomLinearForm");
      const long nvars = NumIndets(P);
      while (true)
      {
        RingElem L(P); // SLUG:  better to use a geobucket???  Or PushFront???
        for (long i=nvars-1; i >=0; --i) L += RandomLong(lo,hi)*indet(P,i);
        if (!IsZero(L)) return L;
      }
      return zero(P); // just to keep the compiler quiet
    }

  } // anonymous namespace -- end
  

  RingElem RandomLinearForm(const ring& P)
  {
    const char* fn("RandomLinearForm(P)");
    if (!IsSparsePolyRing(P))  CoCoA_THROW_ERROR(ERR::ReqSparsePolyRing, fn);
    if (!IsRingFp(CoeffRing(P)))
      CoCoA_THROW_ERROR("only for coefficient ring = Fp", fn);
    return RandomLinearForm(P, 1, ConvertTo<long>(characteristic(P)));
  }


  RingElem RandomLinearForm(const ring& P, long N)
  { return RandomLinearForm(P, -N, N); }


  //-----------------------------------------------------------------
  namespace // for functions local to this file/compilation unit.
  {
    // inline void CheckCompatible(ConstRefRingElem x, ConstRefRingElem y, const char* const FnName)
    // {
    //   if (owner(x) != owner(y))  CoCoA_THROW_ERROR(ERR::MixedRings, FnName);
    // }

    inline void CheckElemSparsePolyRing(ConstRefRingElem f, const ErrorContext& ErrCtx)
    {
      if (!IsSparsePolyRing(owner(f)))
        CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::ReqSparsePolyRing, ErrCtx);
    }

    
    // inline void CheckElemSparsePolyRing(ConstRefRingElem f, const char* const FnName)
    // {
    //   if (!IsSparsePolyRing(owner(f))) CoCoA_THROW_ERROR(ERR::ReqElemSparsePolyRing, FnName);
    // }

    
    void CheckCoeffExpv(const SparsePolyRing& P,
                        ConstRefRingElem c, const std::vector<long>& expv,
                        const char* const FnName)
    {
      if (CoeffRing(P) != owner(c))    CoCoA_THROW_ERROR(ERR::MixedCoeffRings, FnName);
      if (NumIndets(P) != len(expv))   CoCoA_THROW_ERROR(ERR::BadArraySize, FnName);
    }
    
    void CheckCoeffPP(const SparsePolyRing& P,
                      ConstRefRingElem c, ConstRefPPMonoidElem pp,
                      const char* const FnName)
    {
      if (CoeffRing(P) != owner(c))    CoCoA_THROW_ERROR(ERR::MixedCoeffRings, FnName);
      if (PPM(P) != owner(pp))         CoCoA_THROW_ERROR(ERR::MixedPPMs, FnName);
    }
  }
  
  
  RingElem& PushFront(RingElem& f, ConstRefRingElem c, const std::vector<long>& expv) /// SHOULD BE vector<BigInt> ????
  {
    CheckElemSparsePolyRing(f, CoCoA_ERROR_CONTEXT);
    const SparsePolyRing Rx = owner(f);
    CheckCoeffExpv(Rx, c, expv, "PushFront(f, c, expv)");
    PPMonoidElem pp(PPM(Rx), expv);
    if (!IsZero(f) && pp <= LPP(f))
      CoCoA_THROW_ERROR(ERR::PPOrder, "PushFront(f, c, expv)");
    Rx->myPushFront(raw(f), raw(c), raw(pp)); // OK 'cos makes a copy of raw(pp)
    return f;
  }


  RingElem& PushFront(RingElem& f, ConstRefRingElem c, ConstRefPPMonoidElem pp)
  {
    CheckElemSparsePolyRing(f, CoCoA_ERROR_CONTEXT);
    const SparsePolyRing Rx = owner(f);
    CheckCoeffPP(Rx, c, pp, "PushFront(f, c, pp)");
    if (!IsZero(f) && pp <= LPP(f))
      CoCoA_THROW_ERROR(ERR::PPOrder, "PushFront(f, c, pp)");
    Rx->myPushFront(raw(f), raw(c), raw(pp));
    return f;
  }


  RingElem& PushBack(RingElem& f, ConstRefRingElem c, const std::vector<long>& expv) /// SHOULD BE vector<BigInt> ????
  {
    CheckElemSparsePolyRing(f, CoCoA_ERROR_CONTEXT);
    const SparsePolyRing Rx = owner(f);
    CheckCoeffExpv(Rx, c, expv, "PushBack(f, c, expv)");
    PPMonoidElem pp(PPM(Rx), expv);
    Rx->myPushBack(raw(f), raw(c), raw(pp)); // OK 'cos makes a copy of raw(pp)
    return f;
  }


  RingElem& PushBack(RingElem& f, ConstRefRingElem c, ConstRefPPMonoidElem pp)
  {
    CheckElemSparsePolyRing(f, CoCoA_ERROR_CONTEXT);
    const SparsePolyRing Rx = owner(f);
    CheckCoeffPP(Rx, c, pp, "PushBack(f, c, pp)");
    Rx->myPushBack(raw(f), raw(c), raw(pp));
    return f;
  }


  RingElem ClearDenom(const SparsePolyRing& ZZx, const RingElem& f)
  {
    const ring& P = owner(f);
    if (!IsSparsePolyRing(P) ||
        !IsFractionField(CoeffRing(P)) ||
        BaseRing(CoeffRing(P)) != CoeffRing(ZZx) ||
        PPM(P) != PPM(ZZx))
      CoCoA_THROW_ERROR(ERR::BadArg, "ClearDenom(NewRing, f)");
    const RingElem D = CommonDenom(f);
    RingElem ans(ZZx);
    for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it)
      PushBack(ans, num(coeff(it))*(D/den(coeff(it))), PP(it));
    return ans;
  }



  /*----------------------------------------------------------------------
    Member functions every concrete SparsePolyRing implementation
    must have in addition to those of PolyRingBase.
    ----------------------------------------------------------------------*/

  bool SparsePolyRingBase::myIsValid(ConstRawPtr rawf) const
  {
    if (myIsZero(rawf)) return true;
    SparsePolyIter itf = myBeginIter(rawf);
    if (IsZero(coeff(itf))) return false;
    PPMonoidElem PrevPP = PP(itf);
    for (++itf; !IsEnded(itf); ++itf)
    {
      if (IsZero(coeff(itf))) return false;
      if (PrevPP <= PP(itf)) return false;
      PrevPP = PP(itf);
    }
    return true;
  }



  long SparsePolyRingBase::myStdDeg(ConstRawPtr rawf) const
  {
    if (myIsZero(rawf)) CoCoA_THROW_ERROR(ERR::ReqNonZeroRingElem, "myStdDeg(rawf)");
    if (IsStdGraded(ordering(myPPM())))
      return deg(PP(myBeginIter(rawf))); // same as deg(LPP(f))
    long PolyDegree = 0;
    for (SparsePolyIter i=myBeginIter(rawf); !IsEnded(i); ++i)
      PolyDegree = max(PolyDegree, StdDeg(PP(i)));
    return PolyDegree;
  }


  long SparsePolyRingBase::myDeg(ConstRawPtr rawf, long index) const
  {
    if (myIsZero(rawf))  CoCoA_THROW_ERROR(ERR::ReqNonZeroRingElem, "myDeg(rawf, index)");
    if (myNumIndets() == 1)  return deg(PP(myBeginIter(rawf))); // same as deg(LPP(f)), all term-orders are deg compat
    // ???MAYBE: check if 1st row of OrdMat is non-zero only in posn index???
    long res = 0;
    for (SparsePolyIter i=myBeginIter(rawf); !IsEnded(i); ++i)
      res = max(res, exponent(PP(i), index));
    return res;
  }


  void SparsePolyRingBase::myContent(RawPtr rawcontent, ConstRawPtr rawf) const
  {
    const ring& R = myCoeffRing();
    CoCoA_ASSERT(IsTrueGCDDomain(R));
    // Compute answer in local var to avoid aliasing problems; also exception clean.
    RingElem ans(R);
    for (SparsePolyIter i=myBeginIter(rawf); !IsEnded(i); ++i)
    {
      R->myGcd(raw(ans), raw(ans), raw(coeff(i))); // ans = GCD(ans, coeff(i));
      if (IsOne(ans)) break;
    }
    // Finally, swap answer into rawcontent.
    R->mySwap(rawcontent, raw(ans));
    return;
  }

  void SparsePolyRingBase::myContentFrF(RawPtr rawcontent, ConstRawPtr rawf) const
  {
    CoCoA_ASSERT(IsFractionFieldOfGCDDomain(myCoeffRing()));
    const FractionField R = myCoeffRing();
    const ring& S = BaseRing(R);
    RingElem N(S);
    RingElem D(S,1);
    for (SparsePolyIter i=myBeginIter(rawf); !IsEnded(i); ++i)
    {
      N = gcd(N, num(coeff(i)));
      D = lcm(D, den(coeff(i)));
//      S->myGcd(raw(ans), raw(ans), raw(num(coeff(i)))); // ans = GCD(ans, num(coeff(i)));
//      if (IsOne(ans)) break;
    }
    RingHom phi = EmbeddingHom(R);
    RingElem ans = phi(N)/phi(D);
    // Finally, swap answer into rawcontent.
    R->mySwap(rawcontent, raw(ans));
  }


  void SparsePolyRingBase::myCommonDenom(RawPtr rawcontent, ConstRawPtr rawf) const
  {
    CoCoA_ASSERT(IsFractionFieldOfGCDDomain(myCoeffRing()));
    const ring& R = BaseRing(myCoeffRing());
    // Compute result in local "ans" to avoid aliasing problems; also exception clean.
    RingElem ans = one(R);
    for (SparsePolyIter it=myBeginIter(rawf); !IsEnded(it); ++it)
    {
      const RingElem D = den(coeff(it));
      ans *= D/gcd(ans, D);
    }
    // Finally, swap answer into rawcontent -- cheaper than assignment.
    R->mySwap(rawcontent, raw(ans));
  }


  void SparsePolyRingBase::myClearDenom(RawPtr rawg, ConstRawPtr rawf) const
  {
    CoCoA_ASSERT(IsFractionFieldOfGCDDomain(myCoeffRing()));
    const ring& R = BaseRing(myCoeffRing());
    RingElem c(R);
    myCommonDenom(raw(c), rawf);
    const RingElem coeff = EmbeddingHom(myCoeffRing())(c);
    RingElem ans = RingElemAlias(ring(this),rawf);
    myMulByCoeff(raw(ans), raw(coeff));
    // Finally, swap answer into rawg -- cheaper than assignment.
    mySwap(rawg, raw(ans));
  }


  void SparsePolyRingBase::myRemoveBigContent(RawPtr rawf) const
  {
    CoCoA_ASSERT(IsTrueGCDDomain(myCoeffRing()));
    CoCoA_ASSERT(!myIsZero(rawf));
    RingElem cont(myCoeffRing());
    myContent(raw(cont), rawf);
    myDivByCoeff(rawf, raw(cont));
  }


  void SparsePolyRingBase::myMul(RawPtr rawlhs, ConstRawPtr rawf, ConstRawPtr rawg) const
  {
    if (myIsZero(rawf) || myIsZero(rawg)) { myAssignZero(rawlhs); return; }
    if (myIsConstant(rawf))
    {
      RingElem ans(RingElemAlias(ring(this), rawg));
      myMulByCoeff(raw(ans), raw(myLC(rawf)));  // weak exc guarantee, but not a problem here
      mySwap(rawlhs, raw(ans));
      return;
    }
    if (myIsConstant(rawg) && myCoeffRing()->IamCommutative())  // CoeffRing should be comm
    {
      myMul(rawlhs, rawg, rawf);
      return;
    }
    const long gLen = myNumTerms(rawg);
    if (IamCommutative() && myNumTerms(rawf) > gLen) { myMul(rawlhs, rawg, rawf); return; }
    const SparsePolyRing P(this);
    RingElemAlias g(P, rawg);

    if (myIsMonomial(rawf))
    {
      RingElem ans(P);
      myAddMulLM(raw(ans), rawf, rawg);
      mySwap(raw(ans), rawlhs);
      return;
    }
    geobucket gbk(P);
    for (SparsePolyIter itf=myBeginIter(rawf); !IsEnded(itf); ++itf)
    {
      CheckForInterrupt("polynomial multiplication");
      gbk.myAddMulLM(monomial(P, coeff(itf), PP(itf)), g, gLen);
    }
    RingElem ans(P);
    AddClear(ans, gbk); // this is for exception safety
    mySwap(raw(ans), rawlhs); // really an assignment
  }


  bool SparsePolyRingBase::myIsDivisible(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    if (myIsZero(rawy)) { return false; }
    if (myIsZero(rawx)) { myAssignZero(rawlhs); return true; }
    if (myIsConstant(rawy))
    {
      RingElem ans(RingElemAlias(ring(this), rawx));
      if (!myDivByCoeff(raw(ans), raw(myLC(rawy)))) // exc safe?
        return false;
      mySwap(rawlhs, raw(ans));
      return true;
    }
    if (!IamCommutative()) { CoCoA_THROW_ERROR(ERR::NYI, "SparsePolyRingBase::myDiv non commutative"); }
    const SparsePolyRing P(this);
    RingElem xCopy(RingElemAlias(P, rawx));
    geobucket gbk(P);
    gbk.myAddClear(xCopy, NumTerms(xCopy));
    RingElemAlias y(P, rawy);
    const long yLen = NumTerms(y);
    RingElem ans(P);
    const ring& R = myCoeffRing();
    RingElem coeff(R);
    // if not divisible will throw ERR::BadQuot
    while ( !IsZero(gbk) )
    {
      if (!R->myIsDivisible(raw(coeff), raw(LC(gbk)), raw(myLC(rawy)))) return false;
      if (!IsDivisible(LPP(gbk), LPP(y))) return false;      
      RingElem m(monomial(P, coeff, LPP(gbk)/LPP(y)));
      // xCopy -= m*y;
      gbk.myAddMulLM(-m, y, yLen);
      // ans += m;
      myAppendClear(raw(ans), raw(m));
    }
    mySwap(raw(ans), rawlhs); // really an assignment
    return true;
  }


  void SparsePolyRingBase::myDeriv(RawPtr rawlhs, ConstRawPtr rawf, ConstRawPtr rawx) const
  {
    if (myIsOne(rawx)) { myAssign(rawlhs, rawf); return; }
    const long n = myNumIndets();
    ring R = myCoeffRing();
    const SparsePolyRing P(this);
    const vector<long> expD = exponents(myLPP(rawx));

    RingElem ans(P);
    vector<long> expPP(n);
    for (SparsePolyIter itf=myBeginIter(rawf); !IsEnded(itf); ++itf)
    {
      exponents(expPP, PP(itf));
      BigInt scale(1);
      for (long indet=0; indet < n; ++indet)
        if (expD[indet] != 0)
        {
          const long d = expPP[indet];////exponent(PP(itf), indet);
          if (d < expD[indet]) { scale = 0; break; }
          scale *= RangeFactorial(d-expD[indet]+1, d);
        }
      if (IsZero(scale)) continue;
      RingElem m(monomial(P, scale*coeff(itf), PP(itf)/myLPP(rawx)));
      if (!IsZero(m)) myAppendClear(raw(ans), raw(m));
    }
    mySwap(raw(ans), rawlhs); // really an assignment
  }


  void SparsePolyRingBase::myOutput(std::ostream& out, ConstRawPtr rawf) const
  {
    if (!out) return;  // short-cut for bad ostreams

    if (myIsZero(rawf)) { out << '0'; return; }

    const ring& R = myCoeffRing();
    const PPMonoid PPM = myPPM();
    const bool IsQuotientOfZZ = IsQuotientRing(R) && IsZZ(BaseRing(R));
    //const bool IsQuotientOfZ = IsRingFp(R);
    const bool IsNumberRing = IsZZ(R) || IsQuotientOfZZ || IsRingTwinFloat(R); // || IsQQ(R) 

    bool IsFirstCoeff = true;
    for (SparsePolyIter itf=myBeginIter(rawf); !IsEnded(itf) ; ++itf, IsFirstCoeff = false)
    {
      bool PrintStar = true;
      RingElemAlias c = coeff(itf);
      ConstRefPPMonoidElem pp = PP(itf);
      const bool IsWithMinus = R->myIsPrintedWithMinus(raw(c));

      // ---- coefficient ----
      if (!IsFirstCoeff)
      {
        out << ' ';
        if (!IsWithMinus)  out << '+';
      }
      if (IsOne(pp)) { out << c; continue; }
      // ---------- PP != 1 ----------
      if (IsOne(c))
      { // Do not print "1 * "...
        PrintStar = false;
        goto PrintPP;
      }
      if ( IsWithMinus && IsMinusOne(c) ) // in some Z/(n) "-1" prints as n-1
      { // Do not print "-1 * "...
        out << '-';
        PrintStar = false;
        goto PrintPP;
      }
      // General case: coeff is neither +1 nor -1
      if (IsNumberRing || R->myIsPrintAtom(raw(c)) ||
          (IsWithMinus && R->myIsPrintAtom(raw(-c))) )
      {
        out << c;
        goto PrintPP;
      }
      if (!IsFirstCoeff && IsWithMinus) out << '+'; // not printed before
      out << '(' << c << ')';

    PrintPP:      // ---- PP ----
      if (PrintStar)  out << '*';
      out << pp;
    }
  }


  bool SparsePolyRingBase::myIsPrintAtom(ConstRawPtr rawx) const
  {
    long NoUse;
    if (myIsIndet(NoUse, rawx)) return true;
    if (!myIsConstant(rawx)) return false;
    return myCoeffRing()->myIsPrintAtom(raw(myLC(rawx)));
  }


  bool SparsePolyRingBase::myIsPrintedWithMinus(ConstRawPtr rawx) const
  {
    //    if (IsMinusOne(myLC(rawx)) && IsIndet(myLPP(rawx))) return true;
    //    if (!myIsConstant(rawx)) return false;
    return myCoeffRing()->myIsPrintedWithMinus(raw(myLC(rawx)));
  }


  void SparsePolyRingBase::myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawx) const
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("cocoa", "DMPSummands");
    OMOut << myNumTerms(rawx);

    for (SparsePolyIter itf=myBeginIter(rawx); !IsEnded(itf) ; ++itf)
    {
      OMOut << coeff(itf);
      OMOut << PP(itf);
//      R->myOutput_OM(OMOut, it->myCoeff);
//      ordering(PPM)->myOutput_OM(OMOut, it->myOrdv);
    }
    OMOut->mySendApplyEnd();
  }


  bool SparsePolyRingBase::myIsOne(ConstRawPtr rawf) const
  {
    if (!myIsMonomial(rawf)) return false;
    if (!IsOne(myLPP(rawf))) return false;
    return IsOne(myLC(rawf));
  }


  bool SparsePolyRingBase::myIsMinusOne(ConstRawPtr rawf) const
  {
    if (!myIsMonomial(rawf)) return false;
    if (!IsOne(myLPP(rawf))) return false;
    return IsMinusOne(myLC(rawf));
  }


  bool SparsePolyRingBase::myIsConstant(ConstRawPtr rawf) const
  {
    if (myIsZero(rawf)) return true;
    if (!myIsMonomial(rawf)) return false;
    return (IsOne(myLPP(rawf)));
  }


  bool SparsePolyRingBase::myIsIndet(long& IndetIndex, ConstRawPtr rawf) const
  {
    if (!myIsMonomial(rawf)) return false;
    if (!IsOne(myLC(rawf))) return false;
    return IsIndet(IndetIndex, myLPP(rawf));
  }


  bool SparsePolyRingBase::myIsIndetPosPower(long& index, BigInt& EXP, ConstRawPtr rawf) const
  {
    if (!myIsMonomial(rawf)) return false;
    if (!IsOne(myLC(rawf))) return false;
    return IsIndetPosPower(index, EXP, myLPP(rawf));
  }


  bool SparsePolyRingBase::myIsIndetPosPower(ConstRawPtr rawf) const
  {
    if (!myIsMonomial(rawf)) return false;
    if (!IsOne(myLC(rawf))) return false;
    return IsIndetPosPower(myLPP(rawf));
  }


  bool SparsePolyRingBase::myIsEvenPoly(ConstRawPtr rawf) const
  {
    if (myIsZero(rawf)) { return true; }
    for (SparsePolyIter itf=myBeginIter(rawf); !IsEnded(itf); ++itf)
    {
      if ((deg(PP(itf))&1) != 0) return false;
    }
    return true;
  }
  
  bool SparsePolyRingBase::myIsOddPoly(ConstRawPtr rawf) const
  {
    if (myIsZero(rawf)) { return true; }
    for (SparsePolyIter itf=myBeginIter(rawf); !IsEnded(itf); ++itf)
    {
      if ((deg(PP(itf))&1) != 1) return false;
    }
    return true;
  }


  bool SparsePolyRingBase::myIsInvertible(ConstRawPtr rawx) const
  {
    return !myIsZero(rawx) && myIsConstant(rawx) && IsInvertible(myLC(rawx));
  }


  bool SparsePolyRingBase::myIsZeroDivisor(ConstRawPtr rawx) const
  {
    if (myIsZero(rawx)) return true;
    if (IsTrue3(IamIntegralDomain3(true/*quick*/))) return false;
    if (myIsConstant(rawx)) return myCoeffRing()->myIsZeroDivisor(raw(myLC(rawx)));
    if (!myCoeffRing()->myIsZeroDivisor(raw(myLC(rawx)))) return false; // LC is not zd
    CoCoA_THROW_ERROR(ERR::NYI, "SparsePolyRingBase::myIsZeroDivisor when coeffs are not int dom,and LC(poly) is zd ");
    return true; // just to kep compiler quiet
  }


  namespace // anonymous
  {
    // needed for fn immediately after this namespace
    RingElem GCDBySyz(const RingElem& f, const RingElem& g)
    {
         const vector<ModuleElem> v = gens(SyzOfGens(ideal(f,g)));//syz(vector<RingElem>({f,g})) ???
     if ( len(v) != 1 ) 
       CoCoA_THROW_ERROR("Unable to compute GCD", "SparsePolyRingBase::gcd");
     return f/((v[0])[1]);
    }
  } // end of namespace anonymous
  

  // code for R[x,y]: compute gcd in FractionField(R)[x,y]
  void SparsePolyRingBase::myGcd(RawPtr rawlhs, ConstRawPtr rawf, ConstRawPtr rawg) const
  {
    CoCoA_ASSERT(myNumIndets() != 0);
    if ( myIsInvertible(rawf) || myIsInvertible(rawg) ) { myAssign(rawlhs, raw(myOne())); return; }
    if ( myIsZero(rawf) ) { myAssign(rawlhs, rawg); return; }
    if ( myIsZero(rawg) ) { myAssign(rawlhs, rawf); return; }
    const SparsePolyRing P(this);
    RingElemAlias f(P, rawf);
    RingElemAlias g(P, rawg);
    if (IsZZ(myCoeffRing()) || IsQQ(myCoeffRing()))
    {
      RingElem ans = GCD_DMPZ(f, g);
      mySwap(rawlhs, raw(ans));
      return;
    }
    if (!IsField(myCoeffRing()))
    {
      if (!IsTrueGCDDomain(myCoeffRing()))
        CoCoA_THROW_ERROR("NYI gcd of poly with coeffs not in field or TrueGCDDomain", "myGcd");
      else
      {
        FractionField K = NewFractionField(myCoeffRing());
        SparsePolyRing Kx = NewPolyRing_DMP(K, PPM(P));
        RingHom phi = PolyRingHom(P, Kx, CoeffEmbeddingHom(Kx)(EmbeddingHom(K)), indets(Kx));
        RingElem h = gcd(phi(f), phi(g));
        Kx->myClearDenom(raw(h), raw(h));
        RingElem ans(P);
        for (SparsePolyIter it=BeginIter(h) ; !IsEnded(it) ; ++it )
          ans += monomial(P, num(coeff(it)), PP(it));
        myDivByCoeff(raw(ans), raw(content(ans)));
        myMulByCoeff(raw(ans), raw(gcd(content(f), content(g))));
        P->mySwap(rawlhs, raw(ans));
        return;
      }
    }
    // From here onwards myCoeffRing is a *FIELD*
    const long fx = UnivariateIndetIndex(f);
    const long gx = UnivariateIndetIndex(g);  // pointless if fx < 0
    if (myNumIndets() == 1 || (fx >= 0 && fx == gx))
    {
      if (IsRingFp(myCoeffRing()))
      {
        SmallFpImpl ModP(ConvertTo<long>(characteristic(myCoeffRing())));
        RingElem ans = ConvertFromDUPFp(indet(P,fx), gcd(ConvertToDUPFp(ModP, f), ConvertToDUPFp(ModP, g)));
//        RingElem ans = GCD_DUPFF(f,g);
        P->mySwap(rawlhs, raw(ans));
        return;
      }
      // Univariate, coeffs in a field but not a small finite field ==> compute G-basis of principal ideal
      const ideal I = ideal(f, g);
      const vector<RingElem>& GB = GBasis(I);
      CoCoA_ASSERT(len(GB) == 1);
      myAssign(rawlhs, raw(GB[0]));
      return;
    }
    // Multivariate polynomial ==> use syzygy method.
    const RingElem f_indets = IndetsProd(f);
    const RingElem g_indets = IndetsProd(g);
    const RingElem CommonIndets = monomial(P, gcd(LPP(f_indets), LPP(g_indets)));
    if (IsConstant(CommonIndets))  { P->myAssign(rawlhs, raw(CommonIndets)); return; } // it is just 1
    const long nvars = NumIndets(P);
    RingElem f2(f);
    RingElem g2(g);
    if (CommonIndets != f_indets) { const RingElem tmp = f_indets/CommonIndets; vector<long> vars; for (int i=0; i < nvars; ++i) if (deg(tmp,i)>0)  vars.push_back(i);  f2 = ContentWRT(f, vars); }
    if (IsConstant(f2))  { P->myAssign(rawlhs, raw(one(P))); return; } // it is just 1
    if (CommonIndets != g_indets) { const RingElem tmp = g_indets/CommonIndets; vector<long> vars; for (int i=0; i < nvars; ++i) if (deg(tmp,i)>0)  vars.push_back(i);  g2 = ContentWRT(g, vars); }
    if (IsConstant(g2))  { P->myAssign(rawlhs, raw(one(P))); return; } // it is just 1
    if (deg(CommonIndets) == 1)
    {
      // case: multivariate ring, but only 1 indet in common (so f2, g2 are univariate polys)
        SmallFpImpl ModP(ConvertTo<long>(characteristic(myCoeffRing())));
        P->myAssign(rawlhs, raw(ConvertFromDUPFp(CommonIndets,gcd(ConvertToDUPFp(ModP,f2), ConvertToDUPFp(ModP,g2)))));
///LEAKS MEMORY      P->myAssign(rawlhs, raw(GCD_DUPFF(f2,g2)));
      return;
    }
    factorization<RingElem> F = ContentFreeFactor(f2); // all multiplicities are 1
    factorization<RingElem> G = ContentFreeFactor(g2); // all multiplicities are 1
    const vector<RingElem>& facsf = F.myFactors(); // better to sort by IndetsProd?
    const vector<RingElem>& facsg = G.myFactors(); // better to sort by IndetsProd?
    vector<RingElem> IndetsF;  for (const RingElem& ff: facsf)  IndetsF.push_back(IndetsProd(ff));
    vector<RingElem> IndetsG;  for (const RingElem& gg: facsg)  IndetsG.push_back(IndetsProd(gg));
    // Compute GCDs only of matching pairs (if any) -- loop is simple rather than efficient!
    RingElem ans = one(P);
    for (int i=0; i < len(facsf); ++i)
      for (int j=0; j < len(facsg); ++j)
      {
        if (IndetsF[i] != IndetsG[j])  continue;
        ans *= GCDBySyz(facsf[i], facsg[j]);
      }
    P->myAssign(rawlhs, raw(ans));
    return;
  }


  void SparsePolyRingBase::myNormalizeFracNoGcd(RawPtr rawnum, RawPtr rawden) const
  {
    CoCoA_ASSERT(!myIsZero(rawden));
    // Handle case of 0 specially; later code fails otherwise.
    if (myIsZero(rawnum)) { myAssign(rawden, 1); return; }

    const ring& k = myCoeffRing();
    if (IsTrueGCDDomain(k))
    {
      // Coeff ring is a (true) GCD domain
      if (!IsOrderedDomain(k)) return;
      if (myLC(rawden) > 0) return;
      myNegate(rawnum, rawnum);
      myNegate(rawden, rawden);
      return;
    }
    if (IsFractionFieldOfGCDDomain(k))
    {
      // Coeff ring is a FractionField
      const ring R = BaseRing(k);
      const RingHom embed = EmbeddingHom(k);
      RingElemAlias N(ring(this), rawnum);
      RingElemAlias D(ring(this), rawden);

      const RingElem ContN = content(N);
      const RingElem ContD = content(D);
      const RingElem ContQuot = ContN/ContD;
      myMulByCoeff(rawnum, raw(embed(num(ContQuot))/ContN));
      myMulByCoeff(rawden, raw(embed(den(ContQuot))/ContD));

//       RingElem scale(k);
//       scale = embed(lcm(CommonDenom(N), CommonDenom(D)));
//       if (!IsOne(scale))
//       {
//         myMulByCoeff(rawnum, raw(scale));
//         myMulByCoeff(rawden, raw(scale));
//       }
//       scale = embed(gcd(content(N), content(D)));
//       if (!IsOne(scale))
//       {
//         myDivByCoeff(rawnum, raw(scale));
//         myDivByCoeff(rawden, raw(scale));
//       }
      if (!IsOrderedDomain(k)) return;
      if (myLC(rawden) > 0) return;
      myNegate(rawnum, rawnum);
      myNegate(rawden, rawden);
      return;
    }
    // This must come after the case handling FractionFields!
    if (IsField(k))
    {
      if (IsOne(myLC(rawden))) return;
      myDivByCoeff(rawnum, raw(myLC(rawden))); // exc safe?
      myDivByCoeff(rawden, raw(myLC(rawden)));
      return;
    }
    if (IsOrderedDomain(k))
    {
      if (myLC(rawden) > 0) return;
      myNegate(rawnum, rawnum);
      myNegate(rawden, rawden);
      return;
    }
    CoCoA_THROW_ERROR(ERR::NYI, "SparsePolyRing::myNormalizeFracNoGcd");
  }


  bool SparsePolyRingBase::myIsInteger(BigInt& N, ConstRawPtr rawx) const
  {
    if (myIsZero(rawx)) { N = 0; return true; }
    if (myNumTerms(rawx) > 1) return false;
    if (!IsOne(myLPP(rawx))) return false;
    return IsInteger(N, myLC(rawx));
  }


  bool SparsePolyRingBase::myIsRational(BigRat& Q, ConstRawPtr rawx) const
  {
    if (myIsZero(rawx)) { Q = 0; return true; }
    if (myNumTerms(rawx) > 1) return false;
    if (!IsOne(myLPP(rawx))) return false;
    return IsRational(Q, myLC(rawx));
  }


  void SparsePolyRingBase::mySquare(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    if (myIsZero(rawx)) { myAssignZero(rawlhs); return; }
    PowerOverflowCheck(myLPP(rawx),2); // check whether LPP(f)^2 would overflow expv; if not, assume no expv-overflow will occur
    if (IamCommutative())
    {
      const long NT = myNumTerms(rawx);
      if (NT > 4)
      {
        RingElem f(ring(this));
        RingElem g(ring(this));
        long count=0;
        for (SparsePolyIter it=myBeginIter(rawx) ; !IsEnded(it) ; ++it )
        {
          if (++count < NT/2) myPushBack(raw(f), raw(coeff(it)), raw(PP(it)));
          else myPushBack(raw(g), raw(coeff(it)), raw(PP(it)));
        }
        RingElem h = power(f, 2) + power(g, 2) + 2*f*g;
        mySwap(rawlhs, raw(h));
        return;
      }
    }
    mySequentialPower(rawlhs, rawx, 2); //??? BUG/SLUG myBinaryPower better if univariate or coeffs are finite field
  }


  void SparsePolyRingBase::myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const
  {
    // Assert that we have a genuinely non-trivial case.
    CoCoA_ASSERT(n > 1);
    CoCoA_ASSERT(!myIsZero(rawx) && !myIsOne(rawx) && !myIsMinusOne(rawx));

    PowerOverflowCheck(myLPP(rawx),n); // check whether LPP(f)^n would overflow expv; if not, assume no expv-overflow will occur
    // anna 3 apr 2008: under testing
    if (myIsMonomial(rawx) && IamCommutative())
    {
      RingElem m = monomial(SparsePolyRing(this), power(myLC(rawx), n), power(myLPP(rawx), n));
      mySwap(rawlhs, raw(m));
      return;
    }
    if (myIsMonomial(rawx)) { myBinaryPower(rawlhs, rawx, n); return; }
    if (n==2)
    {
      mySquare(rawlhs, rawx);
      return;
    }
    mySequentialPower(rawlhs, rawx, n); //??? BUG/SLUG myBinaryPower better if univariate or coeffs are finite field
  }


  //---- Special functions on RingElem owned by SparsePolyRing

  long UnivariateIndetIndex(ConstRefRingElem f)  // returns < 0 if not univariate
  {
    const SparsePolyRing P = owner(f);
    const long nvars = NumIndets(P);
    if (IsConstant(f))
      CoCoA_THROW_ERROR("Poly must be non-constant","UnivariateIndetIndex"); // f is constant
    if (nvars == 1) return 0; // there is only 1 indet, and f is non-const.
    vector<long> expv = exponents(LPP(f));
    long ans = 0;
    while (expv[ans] == 0)  ++ans;  // f is not constant, so loop is safe!
    for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it)
    {
      exponents(expv, PP(it));
      for (int i=0; i < nvars; ++i)
        if (i != ans && expv[i] != 0) return -1;
    }
    return ans;
  }


  ConstRefPPMonoidElem PPContent(ConstRefRingElem f)     // gcd of all PPs in support
  {
    const SparsePolyRing P = owner(f);
    if (IsZero(f)) return LPP(one(P));
    const int nvars = NumIndets(P);
    vector<long> gcd_exps(nvars);
    vector<long> exps(nvars);
    for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it)
    {
      exponents(exps, PP(it));  // BUG BUG BUG will throw if exps are too big!!!
      bool AllExpsZero = true;
      for (int i=0; i < nvars; ++i)
      {
        if (gcd_exps[i] == 0) continue;
        gcd_exps[i] = min(gcd_exps[i], exps[i]);
        if (gcd_exps[i] != 0) AllExpsZero = false;
      }
      if (AllExpsZero) return LPP(one(P));
    }
    return PPMonoidElem(PPM(P), gcd_exps);
  }


  RingElem IndetsProd(ConstRefRingElem f) // (monomial) prod of indets appearing in f
  {
    const ring& P = owner(f);
    if (!IsSparsePolyRing(P)) CoCoA_THROW_ERROR(ERR::ReqSparsePolyRing, "IndetsProd");
    const vector<long> VarIndices = IndetsIn(f);
    RingElem ans = one(P);
    for (auto k: VarIndices)
      ans *= indet(P,k);
    return ans;
  }


  std::vector<long> IndetsIn(ConstRefRingElem f)
  {
    if (!IsSparsePolyRing(owner(f))) CoCoA_THROW_ERROR(ERR::ReqSparsePolyRing, "IndetsIn");
    const int nvars = NumIndets(owner(f));

    // This approach is cleaner but about 10% slower :-(
    // vector<bool> seen(nvars);
    // const PPMonoid& M = PPM(owner(f));
    // for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it)
    //   M->myIndetsIn(seen, raw(PP(it)));

    long NumIndetsSeen = 0;
    vector<bool> seen(nvars);
    vector<long> exps(nvars);
    for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it)
    {
      exponents(exps, PP(it));  // BUG BUG BUG will throw if exps are too big!!!
      for (int i=0; i < nvars; ++i)
      {
        if (exps[i] == 0 || seen[i]) continue;
        seen[i] = true;
        ++NumIndetsSeen;
      }
      if (NumIndetsSeen == nvars) break;
    }
    // Now convert answer to "list of indices"
    vector<long> ans;
    for (int i=0; i < nvars; ++i)
      if (seen[i]) ans.push_back(i);
    return ans;
  }


  std::ostream& operator<<(std::ostream& out, const CoeffPP& term)
  {
    if (!out) return out;  // short-cut for bad ostreams
    out << "[coeff:=" << term.myCoeff << ", PP:=" << term.myPP << "]" ;
    return out;
  }


  RingElem ConstantCoeff(ConstRefRingElem f)
  {
    // SLUG??? simple impl via iterator; might be able to access last term directly via a pointer???
    const ring& P = owner(f);
    if (!IsSparsePolyRing(P))
      CoCoA_THROW_ERROR(ERR::ReqSparsePolyRing, "ConstantCoeff(f)");
//???    if (IsZero(f)) return zero(CoeffRing(P));
    for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
    {
      if (IsOne(PP(it))) return coeff(it);
    }
    return zero(CoeffRing(P));
  }


  RingElem ContentWRT(ConstRefRingElem f, const std::vector<long>& indets)
  {
    if (!IsSparsePolyRing(owner(f)))
      CoCoA_THROW_ERROR(ERR::ReqSparsePolyRing, "ContentWRT(f,indets)");
    const SparsePolyRing P = owner(f);
    for (long i=0; i < len(indets); ++i)
      if (indets[i] < 0 || indets[i] >= NumIndets(P))
        CoCoA_THROW_ERROR(ERR::BadIndex, "ContentWRT(f,indets)");

    const vector<CoeffPP> CoeffList = CoefficientsWRT(f, indets);
    RingElem ans(P);
    const long n = len(CoeffList);
    for (long i=0; i < n; ++i)
      ans = gcd(ans, CoeffList[i].myCoeff);
    return ans;
  }


  RingElem ContentWRT(ConstRefRingElem f, ConstRefRingElem x)
  {
    if (owner(f) != owner(x))
      CoCoA_THROW_ERROR(ERR::MixedRings, "ContentWRT(f,x)");
    vector<long> indices(1);
    if (!IsIndet(indices[0], x))
      CoCoA_THROW_ERROR(ERR::ReqIndet, "ContentWRT(f,x)");
    return ContentWRT(f, indices);
  }


namespace // anonymous for file local defns
{

  BigRat RatReconstruct(BigInt X, BigInt M)
  {
    RatReconstructByContFrac reconstructor; // default ctor arg
    reconstructor.myAddInfo(X, M);
    if (!IsConvincing(reconstructor))
      CoCoA_THROW_ERROR(ERR::CannotReconstruct, "RatReconstruct");
    // std::cout << " rat = " << ReconstructedRat(reconstructor) << std::endl;
    return ReconstructedRat(reconstructor);
  }
  
  
  BigInt AsINT(ConstRefRingElem c)
  {
    BigInt i;
    if (!IsInteger(i, c)) CoCoA_THROW_ERROR("not integer", "AsINT");
    return i;
  }


  BigInt CRT_modulus(BigInt M1, BigInt M2)
  {
    CRTMill CRT;
    CRT.myAddInfo(BigInt(0), M1);
    CRT.myAddInfo(BigInt(0), M2);
    return CombinedModulus(CRT);
  }

  BigInt CRT4(BigInt R1, BigInt M1, BigInt R2, BigInt M2)
  {
    CRTMill CRT;
    CRT.myAddInfo(R1, M1);
    CRT.myAddInfo(R2, M2);
    return CombinedResidue(CRT);
  }

} // end of anonymous namespace

  void CRTPoly(RingElem& f, BigInt& M,
               ConstRefRingElem f1, const BigInt& M1,
               ConstRefRingElem f2, const BigInt& M2)
  {
    const SparsePolyRing P = owner(f1);
    RingElem ans(P);
    if (!IsQQ(CoeffRing(P)))
      CoCoA_THROW_ERROR("reconstruction only into QQ","RatReconstructPoly");
    SparsePolyIter it1 = BeginIter(f1);
    SparsePolyIter it2 = BeginIter(f2);
    const BigInt Zero;
    while ( !(IsEnded(it1)||IsEnded(it2)) )
      if (PP(it1) == PP(it2))
      {  
        ans += monomial(P, CRT4(AsINT(coeff(it1)),M1, AsINT(coeff(it2)),M2), PP(it1));
        ++it1;
        ++it2;
      }
      else 
      {
        if (PP(it1) < PP(it2))
        {
          ans += monomial(P, CRT4(Zero, M1, AsINT(coeff(it2)), M2), PP(it2));
          ++it2;
        }
        else // LPP(f1) > LPP(f2)
        {
          ans += monomial(P, CRT4(AsINT(coeff(it1)), M1, Zero, M2), PP(it1));
          ++it1;
        }
      }
    for (; !IsEnded(it1); ++it1)
      ans += monomial(P, CRT4(AsINT(coeff(it1)), M1, Zero, M2), PP(it1));
    for (; !IsEnded(it2); ++it2)
      ans += monomial(P, CRT4(Zero, M1, AsINT(coeff(it2)), M2), PP(it2));
    swap(f, ans);
    //    std::cout << " CRT f = " << f << std::endl;
    M = CRT_modulus(M1, M2);   // very likely M1*M2
  }


  RingElem RatReconstructPoly(ConstRefRingElem fCRT, const BigInt& modulus)
  {
    const SparsePolyRing P(owner(fCRT));
    if (!IsQQ(CoeffRing(P)))
      CoCoA_THROW_ERROR("reconstruction only into QQ","RatReconstructPoly");
    RingElem f(P);
    for (SparsePolyIter it = BeginIter(fCRT); !IsEnded(it); ++it)
      f += monomial(P, RatReconstruct(AsINT(coeff(it)), modulus), PP(it));
    return f;
  }


  //--- move to PolyOps-RingElem.C ???
  RingElem CoeffHeight(ConstRefRingElem f)
  {
    const char* const FnName = "CoeffHeight";
    const ring& P = owner(f);
    if (!IsPolyRing(P))
      CoCoA_THROW_ERROR(ERR::ReqPolyRing, FnName);
    if (!IsOrderedDomain(CoeffRing(P)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, FnName);
    RingElem ans = zero(CoeffRing(P));
    // Loop below makes wasteful copies of the coeffs.
    for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
    {
      const RingElem c = abs(coeff(it));
      if (c > ans)
        ans = c;
    }
    return ans;
  }


  // Just for univariate -- requires non-zero const term
  // A bit ugly because I wanted to avoid copying any coeffs.
  bool IsPalindromic(ConstRefRingElem f)
  {
    const ring& Rx = owner(f);
    if (!IsSparsePolyRing(Rx))
      CoCoA_THROW_ERROR(ERR::ReqSparsePolyRing, "IsPalindromic");
    if (IsConstant(f)) return true;
    if (UnivariateIndetIndex(f) < 0)
      CoCoA_THROW_ERROR(ERR::ReqUnivariate, "IsPalindromic");
    const long n = NumTerms(f);
    if (n == 1 || IsZero(ConstantCoeff(f))) return false;
    const long degf = deg(f);
    if (IsOdd(degf) && IsOdd(n)) return false;
    const ring& R = CoeffRing(Rx);
    vector<RingElemConstRawPtr> C; C.reserve(n/2); // integer division!
    vector<long> exp; exp.reserve(n/2);            // integer division!
    long count = 0;
    for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
    {
      ++count;
      if (count <= n/2)
      {
        // first half: we just record exps & coeffs
        C.push_back(raw(coeff(it)));
        const long d = deg(PP(it));
        if (d <= degf/2) return false;
        exp.push_back(d);
        continue;
      }
      // We are in second half of poly, so check against first half
      if (IsOdd(n) && count == (n+1)/2)  // integer division!
      {
        if (deg(PP(it)) != degf/2) return false;
        continue;// no need to check middle coeff
      }
      if (deg(PP(it)) != degf-exp[n-count]) return false;
      if (coeff(it) != RingElemAlias(R, C[n-count])) return false;
    }
    return true;
  }


  // univariate case
  RingElem reverse(ConstRefRingElem f)
  {
    const ring& P = owner(f);
    if (!IsSparsePolyRing(P))
      CoCoA_THROW_ERROR(ERR::ReqSparsePolyRing, "reverse");
    if (IsZero(f)) return f;
    if (UnivariateIndetIndex(f) < 0)
      CoCoA_THROW_ERROR(ERR::ReqUnivariate, "reverse");
    return reverse(f, LPP(f));
  }

  RingElem reverse(ConstRefRingElem f, ConstRefPPMonoidElem t)
  {
    const SparsePolyRing P = owner(f);
    if (!IsSparsePolyRing(P))
      CoCoA_THROW_ERROR(ERR::ReqSparsePolyRing, "reverse");
    if (IsZero(f)) return f;
    RingElem ans = zero(P);
    for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it)
    {
      P->myPushFront(raw(ans), raw(coeff(it)), raw(t/PP(it)));
// 2024-02-16 PREVIOUSLY:      PushFront(ans, coeff(it), t/PP(it));
    }
    return ans;
  }


  // Standard graeffe transformation (squares the roots)
  RingElem graeffe(ConstRefRingElem f)
  {
    const SparsePolyRing P(owner(f), CoCoA_ERROR_CONTEXT);
////    if (!IsSparsePolyRing(P))
////      CoCoA_THROW_ERROR(ERR::ReqSparsePolyRing, "graeffe");
    if (IsZero(f))  return f;
    if (deg(f) == 0)  return one(P); //?????
    const long IndetIndex = UnivariateIndetIndex(f);
    if (IndetIndex < 0)
      CoCoA_THROW_ERROR(ERR::ReqUnivariate, "graeffe");
    PPMonoidElem pp = indet(PPM(P),IndetIndex);
    RingElem EvenPart(P);
    RingElem OddPart(P);
    for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
    {
      const int d = deg(PP(it));
      if (IsEven(d))
        PushBack(EvenPart, coeff(it), power(pp,d/2));
      else
        PushBack(OddPart, coeff(it), power(pp,d/2)); // integer division!
    }
    const RingElem& x = indet(P, IndetIndex);
    // Sign of result taken from Wikipedia entry for "Graeffe's Method"
    if (IsEven(deg(f)))
      return power(EvenPart,2) - x*power(OddPart,2);
    else
      return x*power(OddPart,2) - power(EvenPart,2);
  }


  // Cubic graeffe transformation  (cubes the roots)
  // FORMULA:  f0^3 -3*f0*f1*f2*x +f1^3*x +f2^3*x^2
  RingElem graeffe3(ConstRefRingElem f)
  {
    const ring& P = owner(f);
    if (!IsSparsePolyRing(P))
      CoCoA_THROW_ERROR(ERR::ReqSparsePolyRing, "graeffe");
    if (IsZero(f)) return f;
    const long IndetIndex = UnivariateIndetIndex(f);
    if (IndetIndex < 0)
      CoCoA_THROW_ERROR(ERR::ReqUnivariate, "graeffe");

    PPMonoidElem pp = indet(PPM(P),IndetIndex);
    vector<RingElem> part(3, zero(P));
    for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
    {
      const int d = deg(PP(it));
      PushBack(part[d%3], coeff(it), power(pp, d/3)); // integer division!
    }
    const RingElem& x = indet(P, IndetIndex);
    const RingElem part012 = part[0]*part[1]*part[2];
    RingElem ans = power(part[0],3) + x*(power(part[1],3)-3*part012 + x*power(part[2],3));
    return ans;
  }


} // end of namespace CoCoA

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SparsePolyOps-RingElem.C,v 1.53 2024/08/01 16:51:46 bigatti Exp $
// $Log: SparsePolyOps-RingElem.C,v $
// Revision 1.53  2024/08/01 16:51:46  bigatti
// Summary: error codes: Req(Elem)SparsePolyRing;
// now using CoCoA_THROW_ERROR1/2
//
// Revision 1.52  2024/07/04 06:04:47  bigatti
// Summary: bug fix: missing one ":" in error code "ERR::.."
//
// Revision 1.51  2024/07/03 16:07:44  bigatti
// Summary: Error codes (Redmine #427): changed some string messages into equivalent error codes
//
// Revision 1.50  2024/07/03 09:35:15  bigatti
// Summary: error codes (redmine 92):
// removed BadIntedIndex, BadDegIndex, BadComptIndex
//
// Revision 1.49  2024/07/02 20:42:07  abbott
// Summary: Renaming of errors (redmine 308)
//
// Revision 1.48  2024/07/02 09:57:58  bigatti
// Summary: error codes, first batch:
// ReqUnivariate, ReqNonZero, ReqNonZeroGradingDim, ReqNonZeroModulus,  ReqNonZeroRingElem
//
// Revision 1.47  2024/05/03 12:15:49  abbott
// Summary: Fixed memory leak (redmine 1814)
//
// Revision 1.46  2024/03/25 16:23:13  bigatti
// Summary: removed useless includes, removed code for
// Coefficients... (now in SparsePolyOps-vector.C)
//
// Revision 1.45  2024/03/21 20:09:05  abbott
// Summary: Cleaned up after resolving redmine 1641
//
// Revision 1.44  2024/03/15 18:57:32  abbott
// Summary: Improved use of PolyRing inside fns; employed CoCoA_ERROR_CONTEXT
//
// Revision 1.43  2024/02/23 10:50:36  abbott
// Summary: Made reverse a bit faster
//
// Revision 1.42  2023/10/05 18:02:55  abbott
// Summary: Corrected defn of graeffe (using sign from Wikipedia); improved deg (redmine 1756)
//
// Revision 1.41  2023/08/06 10:14:11  abbott
// Summary: Removed incorrect "sign trick" (redmine 1758)
//
// Revision 1.40  2023/06/16 07:55:57  abbott
// Summary: Added short-cut for deg of poly in univariate ring
//
// Revision 1.39  2023/03/14 22:35:22  abbott
// Summary: New iterator for loop syntax (redmine 1436)
//
// Revision 1.38  2023/03/09 22:41:07  abbott
// Summary: Tidied use of exponents (redmine 1694)
//
// Revision 1.37  2022/11/30 15:24:11  abbott
// Summary: Minor update (better readability)
//
// Revision 1.36  2022/11/26 13:59:47  abbott
// Summary: Fixed IsPalindromic (was badly buggy: wrong ring in RingElemAlias)
//
// Revision 1.35  2022/10/11 15:44:00  bigatti
// Summary: using new interface of exponents returning vector
//
// Revision 1.34  2022/08/08 19:02:35  abbott
// Summary: Cleaning & better comments
//
// Revision 1.33  2022/03/09 17:17:10  abbott
// Summary: Changed args to const ref for CRTPoly and RatReconstructPoly
//
// Revision 1.32  2022/03/09 07:57:07  bigatti
// Summary: added TmpDivAlg
//
// Revision 1.31  2022/03/07 10:57:52  bigatti
// Summary: commented out MaxDegB (unused)
//
// Revision 1.30  2022/02/18 14:11:59  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.29  2022/02/14 14:51:11  bigatti
// Summary: vector of RingElem functions moved to SparsePolyOps-vector
//
// Revision 1.28  2022/02/08 20:18:55  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.27  2022/02/02 09:41:19  abbott
// Summary: Added PPContent
//
// Revision 1.26  2022/01/07 19:16:15  abbott
// Summary: First gcd improvement (redmine 1641)
//
// Revision 1.25  2021/11/28 10:05:03  abbott
// Summary: Added arg checks to NR (redmine 1272, comment 2)
//
// Revision 1.24  2021/11/19 16:33:24  abbott
// Summary: Made poly moultiplication interruptible (redmine 1633)
//
// Revision 1.23  2021/10/30 19:35:59  abbott
// Summary: Commented out unused fn CheckCompatible
//
// Revision 1.22  2021/10/20 18:51:04  abbott
// Summary: Revised semantics of UnivariateIndetIndex (redmine 1617)
//
// Revision 1.21  2021/10/14 12:11:22  abbott
// Summary: Corrected nvars check in UnivariateIndetIndex
//
// Revision 1.20  2020/12/04 10:45:07  abbott
// Summary: Made reduction interruptible
//
// Revision 1.19  2020/10/27 10:01:50  abbott
// Summary: Changed comment
//
// Revision 1.18  2020/06/17 15:49:28  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.17  2020/04/19 17:15:07  abbott
// Summary: Added ConstantCoeff
//
// Revision 1.16  2020/03/11 17:00:27  abbott
// Summary: Added new fn HomogCompt; split of fns to do with homog polys into new file.  Cleaned includes.
//
// Revision 1.15  2020/03/04 20:05:40  abbott
// Summary: Changed name CoeffVecWRTBasis to CoeffVecWRTSupport
//
// Revision 1.14  2020/02/27 14:00:28  abbott
// Summary: Added CoeffVecWRTBasis (or CoeffListWRTBasis in CoCoA-5)
//
// Revision 1.13  2020/02/14 12:22:37  abbott
// Summary: Removed old undocumented fn indets(RingElem); merged impl into IndetsIn
//
// Revision 1.12  2020/02/13 13:54:57  abbott
// Summary: Added IndetsProd
//
// Revision 1.11  2020/02/12 15:38:54  bigatti
// -- made new file SparsePolyIter.C for BeginIter and EndIter
//
// Revision 1.10  2020/02/11 16:56:42  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.9  2020/02/11 16:12:19  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.8  2020/01/26 14:37:24  abbott
// Summary: Removed useless include (of NumTheory)
//
// Revision 1.7  2019/12/28 17:53:53  abbott
// Summary: Added myIsIndetPosPower (with 3 args)
//
// Revision 1.6  2019/10/29 11:38:03  abbott
// Summary: Added IndetsIn also for vector<RingElem>
//
// Revision 1.5  2019/09/16 17:38:19  abbott
// Summary: Impl myIsEvenPoly, myIsOddPoly
//
// Revision 1.4  2019/03/18 11:23:21  abbott
// Summary: Added new include after splitting NumTheory
//
// Revision 1.3  2018/08/06 13:38:20  bigatti
// -- added functins from SparsePolyOps.C
//
// Revision 1.2  2018/07/26 15:25:38  abbott
// Summary: Added comment about SLUG (not using geobucket in RandomLinearForm)
//
// Revision 1.1  2018/05/18 16:35:52  bigatti
// -- split SparsePolyOps-RingElem from SparsePolyRing
//
