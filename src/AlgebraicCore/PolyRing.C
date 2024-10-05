//   Copyright (c)  2003-2018  John Abbott and Anna M. Bigatti

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


// Source code for abstract class PolyRing and friends

#include "CoCoA/PolyRing.H"

#include "CoCoA/CanonicalHom.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/NumTheory-gcd.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/SparsePolyOps-eval.H"
//#include "CoCoA/SparsePolyRing.H" // from SparsePolyOps-RingElem.H
#include "CoCoA/VectorOps.H"    // for HasUniqueOwner
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/factor.H"  // for IsIrredPoly
#include "CoCoA/interrupt.H"
#include "CoCoA/symbol.H"  // for myIndetsCalled
#include "CoCoA/utils.H"  // for len

#include <functional>
using std::not1;    // for AreMonomials
using std::ptr_fun; // for AreMonomials
//#include <vector>
using std::vector;


namespace CoCoA
{

  BigInt PolyRingBase::myCharacteristic() const
  { return characteristic(myCoeffRing()); }  // not inline as o/w requires more includes in .H file

  // void PolyRingBase::myCheckIndetIndex(long i, const char* where) const
  // {
  //   if (i < 0 || i >= myNumIndets())
  //     CoCoA_THROW_ERROR(ERR::BadIndetIndex, where);
  // }

  void PolyRingBase::myCheckIndetIndex(long i, const ErrorContext& ErrCtx) const
  {
    if (i < 0 || i >= myNumIndets())
      CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::BadIndex, ErrCtx);
  }

  bool PolyRingBase::myIsIrred(ConstRawPtr rawf) const
  {
    if (myIsZero(rawf)) CoCoA_THROW_ERROR(ERR::ReqNonZeroRingElem, "RingBase::myIsIrred(rawx)");
    if (myIsInvertible(rawf)) CoCoA_THROW_ERROR(ERR::InvertibleRingElem, "RingBase::myIsIrred(rawx)");
    return IsIrredPoly(RingElemAlias(ring(this), rawf));
  }


  void PolyRingBase::myDiv(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    if (myIsZero(rawy)) // or CoCoA_ASSERT???
      CoCoA_THROW_ERROR(ERR::DivByZero,"PolyRingBase::myDiv");
    if (!myIsDivisible(rawlhs, rawx, rawy))
      CoCoA_THROW_ERROR(ERR::BadQuot, "PolyRingBase::myDiv");
  }  
  

  void PolyRingBase::myMonic(RawPtr rawmonic, ConstRawPtr rawf) const
  {
    if (myIsZero(rawf)) // or CoCoA_ASSERT???
      CoCoA_THROW_ERROR(ERR::ReqNonZeroRingElem,"PolyRingBase::myMonic");
    RingElem ans = RingElemAlias(ring(this), rawf);
    if (!IsOne(myLC(rawf)) && !myDivByCoeff(raw(ans), raw(myLC(rawf))))
      CoCoA_THROW_ERROR(ERR::BadQuot, "PolyRingBase::myDiv");
    mySwap(rawmonic, raw(ans));
  }
  

  bool IsSqFree(ConstRefRingElem f)
  {
    const ring& P = owner(f);
    if (!IsPolyRing(P)) CoCoA_THROW_ERROR(ERR::ReqPolyRing, "IsSqFree");
    if (!IsZZ(CoeffRing(P)) && !IsField(CoeffRing(P)))
      CoCoA_THROW_ERROR("CoeffRing must be ZZ or a field", "IsSqFree");
    if (IsZero(f)) CoCoA_THROW_ERROR(ERR::ReqNonZero, "IsSqFree");
    if (IsConstant(f)) return true;

    // Special case if poly is actually a single term
    if (IsMonomial(f))
      return IsSqFree(LPP(f));

    // Distinguish univariate and multivariate...
    const RingElem xyz = IndetsProd(f);
    const bool multivariate = (deg(xyz) > 1);
    if (!multivariate)
    {
      const long x = UnivariateIndetIndex(f); // f *is* univariate!
      return IsConstant(gcd(f, deriv(f,x))); // NOT IsCoprime in case coeffring is not field (e.g. ZZ) [redmine 1710]
    }
    
    // General multivariate case: pick an indet x, test content wrt x, and then the rest.
    const vector<long> expv = exponents(LPP(f));
    long x=0;
    while (expv[x] == 0)
      ++x;
    
    const RingElem content = ContentWRT(f, indet(P,x));
    if (!IsSqFree(content)) return false;
    const RingElem g = f/content;
    return IsConstant(gcd(g, deriv(g,x))); // NOT IsCoprime in case coeffring is not field (e.g. ZZ) [redmine 1710]
  }


  bool PolyRingBase::myImageLiesInSubfield(const RingHom& /*phi*/) const
  {
    CoCoA_THROW_ERROR(ERR::NYI, "PolyRingBase::myImageLiesInSubfield");
    return false; // just to keep compiler quiet
  }


  std::vector<RingElem> PolyRingBase::myIndets(const std::string& s) const
  {
    std::vector<RingElem> inds;
    for (long i=0; i < myNumIndets(); ++i)
      if (head(myIndetSymbol(i)) == s)
        inds.push_back(myIndets()[i]);
    return inds;
  }
  

  void PolyRingBase::myOutputSelf(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams

    out << "RingWithID(" << myID << ", \"";
    myCoeffRing()->myOutputSelfShort(out);
    out << "[" << myIndets()[0];
    for (long i=1; i<myNumIndets(); ++i)  out << "," << myIndets()[i];
    out <<"]\")";
  }


  void PolyRingBase::myOutputSelfLong(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams

    out << "RingWithID(" << myID
        << ", \"RingWithID(" << RingID(myCoeffRing()) << ")[" << myIndets()[0];
    for (long i=1; i<myNumIndets(); ++i)  out << "," << myIndets()[i];
    out <<"] -- " << myImplDetails() << "\")\n  with CoeffRing ";
    myCoeffRing()->myOutputSelfLong(out);
  }


  namespace // anonymous
  {
    void RingQQtDtor(void* ptr)
    {
      delete static_cast<PolyRing*>(ptr);
    }
  } // end of namespace anonymous

  const PolyRing& RingQQt(const MachineInt& NumIndets)
  {
    static vector<PolyRing*> QQtTable;
    if (IsNegative(NumIndets) || IsZero(NumIndets)) CoCoA_THROW_ERROR(ERR::ReqPositive, "RingQQt");
    const long n = AsSignedLong(NumIndets);
    if (n >= len(QQtTable)) QQtTable.resize(n+1); // will fill with NULL ptrs
    if (QQtTable[n] == nullptr)
    {
      vector<symbol> IndetNames;
      if (n == 1) IndetNames = symbols("t"); else IndetNames = SymbolRange("t",1,n);
      QQtTable[n] = new SparsePolyRing(NewPolyRing(RingQQ(), IndetNames)); // wasteful copy!!
      RegisterDtorForGlobal(&RingQQtDtor, QQtTable[n]);
    }
    return *QQtTable[n];
  }



  const RingElem& indet(const PolyRing& P, long var)
  {
    P->myCheckIndetIndex(var, CoCoA_ERROR_CONTEXT);
    return P->myIndets()[var];
  }


  RingElem IndetPower(const PolyRing& P, long var, long exp)  // error if exp < 0
  {
    P->myCheckIndetIndex(var, CoCoA_ERROR_CONTEXT);
    if (exp < 0) CoCoA_THROW_ERROR(ERR::NegExp, "IndetPower(P, var, exp)");
    RingElem ans(P, 1);
    P->myIndetPower(raw(ans), var, exp);
    return ans;
  }


  RingElem IndetPower(const PolyRing& P, long var, const BigInt& EXP)  // error if EXP < 0
  {
    return IndetPower(P, var, ConvertTo<long>(EXP));
  }


  // long NumTerms(ConstRefRingElem f)
  // {
  //   if (!IsPolyRing(owner(f)))
  //     CoCoA_THROW_ERROR(ERR::NotElemPolyRing, "NumTerms(f)");
  //   return PolyRingPtr(owner(f))->myNumTerms(raw(f));
  // }

  long NumTerms(ConstRefRingElem f)
  { return PolyRingPtr(owner(f), CoCoA_ERROR_CONTEXT)->myNumTerms(raw(f)); }


  bool IsMonomial(ConstRefRingElem f)
  {
    if (!IsPolyRing(owner(f)))
      CoCoA_THROW_ERROR(ERR::ReqElemPolyRing, "IsMonomial(f)");
    return PolyRingPtr(owner(f))->myIsMonomial(raw(f));
  }


  bool AreMonomials(const std::vector<RingElem>& v)
  {
    // morally:  return find_if(v.begin(), v.end(), not1(IsMonomial)) == v.end();
    if (!HasUniqueOwner(v)) CoCoA_THROW_ERROR(ERR::MixedRings, "AreMonomials(v)");
    const long n = len(v);
    for (long i=0; i < n; ++i)
      if (!IsMonomial(v[i])) return false;
    return true;
//  We *DO NOT USE* STL algorithm because ptr_fun fails when args are references.
//     return find_if(v.begin(), v.end(),
//                    not1(ptr_fun(static_cast<bool(*)(const RingElemAlias&)>(CoCoA::IsMonomial))))
//       == v.end(); 
  }


  bool IsConstant(ConstRefRingElem f)
  {
    if (!IsPolyRing(owner(f)))
      CoCoA_THROW_ERROR(ERR::ReqElemPolyRing, "IsConstant(f)");
    return PolyRingPtr(owner(f))->myIsConstant(raw(f));
  }


  bool IsIndet(ConstRefRingElem f)
  {
    if (!IsPolyRing(owner(f)))
      CoCoA_THROW_ERROR(ERR::ReqElemPolyRing, "IsIndet(f)");
    long junk;
    return PolyRingPtr(owner(f))->myIsIndet(junk, raw(f));
  }


  bool IsIndet(long& index, ConstRefRingElem f)
  {
    if (!IsPolyRing(owner(f)))
      CoCoA_THROW_ERROR(ERR::ReqElemPolyRing, "IsIndet(index,f)");
    return PolyRingPtr(owner(f))->myIsIndet(index, raw(f));
  }


  bool IsIndetPosPower(ConstRefRingElem f)
  {
    if (!IsPolyRing(owner(f)))
      CoCoA_THROW_ERROR(ERR::ReqElemPolyRing, "IsIndetPosPower(f)");
    return PolyRingPtr(owner(f))->myIsIndetPosPower(raw(f));
  }


  bool IsEvenPoly(ConstRefRingElem f)
  {
    if (!IsPolyRing(owner(f)))
      CoCoA_THROW_ERROR(ERR::ReqElemPolyRing, "IsEvenPoly(f)");
    return PolyRingPtr(owner(f))->myIsEvenPoly(raw(f));
  }
  
  bool IsOddPoly(ConstRefRingElem f)
  {
    if (!IsPolyRing(owner(f)))
      CoCoA_THROW_ERROR(ERR::ReqElemPolyRing, "IsOddPoly(f)");
    return PolyRingPtr(owner(f))->myIsOddPoly(raw(f));
  }


  long deg(ConstRefRingElem f, long var)
  {
    const char* const FnName = "deg(f,var)";
    if (!IsPolyRing(owner(f)))
      CoCoA_THROW_ERROR(ERR::ReqElemPolyRing, FnName);
    const PolyRingBase* P = PolyRingPtr(owner(f));
    P->myCheckIndetIndex(var, CoCoA_ERROR_CONTEXT);
    if (IsZero(f))
      CoCoA_THROW_ERROR(ERR::ReqNonZeroRingElem, FnName);
    return P->myDeg(raw(f), var);
  }


  // This impl allows content to be greater than 1 (but requires integer coeffs)
  RingElem FixedDivisor(ConstRefRingElem f)
  {
    // hint: could be clever if f is even or odd (need do only half as many evals)
    const char* const FnName = "FixedDivisor";
///    const ring& P = owner(f);
    const PolyRing P(owner(f), CoCoA_ERROR_CONTEXT);
//    if (!IsPolyRing(P))  CoCoA_THROW_ERROR(ERR::ReqPolyRing, FnName);
    if (!IsZero(characteristic(P)))  CoCoA_THROW_ERROR(ERR::BadArg, "FixedDivisor: characteristic must be 0");
    if (IsZero(f))  CoCoA_THROW_ERROR(ERR::ReqNonZeroRingElem, "FixedDivisor");
    if (deg(f) == 0)  return abs(LC(f));
    const long x_ind = UnivariateIndetIndex(f);
    if (x_ind < 0)  CoCoA_THROW_ERROR(ERR::ReqUnivariate, "FixedDivisor");
    const ring& k = CoeffRing(P);
    const BigInt cont = ConvertTo<BigInt>(content(f)); // requires coeffs to be integers!
    EvalUPoly F(f/cont);
    vector<RingElem> EvalPt(NumIndets(P), zero(k)); // initially all zeroes
    long x_lo = 0;
    BigInt val = F(0);
    if (IsOne(val))  return RingElem(k,cont);
    BigInt IC = gcd(val, factorial(deg(f)));  // wasteful if val0 is small!
    if (IsOne(IC))  return RingElem(k,cont);

    long x_hi = 0;
    BigInt val_lo = val;
    BigInt val_hi = val;
    long width=1;
    BigInt fact(1); // factorial(width)
    while (true)
    {
      CheckForInterrupt(FnName);
      ++width;
      fact *= width;
      const bool MoveHi = (val_lo >= val_hi);
      if (MoveHi)
      {
        ++x_hi;
        val = F(x_hi);
      }
      else
      {
        --x_lo;
        val = F(x_lo);
      }
      IC = gcd(IC, val);
      if (IsOne(IC))
        return  RingElem(k,cont);
      if (IC <= fact  &&  gcd(IC,fact) == IC)
        return RingElem(k, IC*cont);
      if (MoveHi)
      {
        val_hi = abs(val);
      }
      else
      {
        val_lo = abs(val);
      }
    }
  }


  namespace // anonymous for file local fn
  {
    RingElem DerivFrF(ConstRefRingElem f, ConstRefRingElem x)
    {
      const FractionField FrF = owner(f);
      if (!IsPolyRing(BaseRing(FrF))) CoCoA_THROW_ERROR(ERR::ReqPolyRing, "deriv(f,x)");
      if (!IsOne(den(x))) CoCoA_THROW_ERROR("Nonsensical derivative", "deriv(f,x)");
      RingElem ans(FrF);
      FrF->myDeriv(raw(ans), raw(f), raw(x));
      return ans;
    }
  } // end of anonymous namespace

  RingElem deriv(ConstRefRingElem f, ConstRefRingElem x)
  {
    if (owner(x) != owner(f)) CoCoA_THROW_ERROR(ERR::MixedRings, "deriv(f, x)");
    if (IsFractionField(owner(f))) return DerivFrF(f,x);
    // From here on we are in the "polynomial" case.
    if (!IsIndet(x)) CoCoA_THROW_ERROR(ERR::ReqIndet, "deriv(f,x)");
    const PolyRing Rx = owner(f);
    RingElem ans(Rx);
    Rx->myDeriv(raw(ans), raw(f), raw(x));
    return ans;
  }


  RingElem deriv(ConstRefRingElem f, long x) // here x is the index of the variable
  {
    const PolyRing Rx = owner(f);
    Rx->myCheckIndetIndex(x, CoCoA_ERROR_CONTEXT);
    return deriv(f, indet(Rx, x));
  }


  RingHom CoeffEmbeddingHom(const PolyRing& Rx)
  {
    return Rx->myCoeffEmbeddingHomCtor();
  }


  // Rx is the domain, S is the codomain
  RingHom PolyRingHom(const PolyRing& Rx, const ring& S, RingHom CoeffHom, const std::vector<RingElem>& IndetImages)
  {
    const std::string FnName = "PolyRingHom(Rx,S,CoeffHom,IndetImages): ";
    if (domain(CoeffHom) != CoeffRing(Rx))
      CoCoA_THROW_ERROR(ERR::BadDomain, FnName + "argument CoeffHom");
    if (IsPolyRing(S) && codomain(CoeffHom) == CoeffRing(S))
      CoeffHom = CoeffEmbeddingHom(S)(CoeffHom);
    if (codomain(CoeffHom) != S)
      CoCoA_THROW_ERROR(ERR::BadCodomain, FnName + "argument CoeffHom");
    if (NumIndets(Rx) != len(IndetImages))
      CoCoA_THROW_ERROR(ERR::BadArraySize, FnName + "arguments Rx, IndetImages");
    for (long i=0; i < NumIndets(Rx); ++i)
      if (owner(IndetImages[i]) != S)
        CoCoA_THROW_ERROR(ERR::BadPolyRingHomImages, FnName + "argument IndetImages");

    return Rx->myHomCtor(S, CoeffHom, IndetImages);
  }


  // Rx is the domain, S is the codomain
  RingHom PolyRingHom(const PolyRing& Rx, const ring& S, RingHom CoeffHom, const std::string& IndetImages)
  { return PolyRingHom(Rx, S, CoeffHom, RingElems(S, IndetImages)); }
  

  RingHom EvalHom(const PolyRing& Rx, const std::vector<RingElem>& IndetImages)
  {
    const char* const FnName = "EvalHom(Rx,IndetImages)";
    const ring& R = CoeffRing(Rx);
    if (NumIndets(Rx) != len(IndetImages))
      CoCoA_THROW_ERROR(ERR::BadArraySize, FnName);
    for (long i=0; i < NumIndets(Rx); ++i)
      if (owner(IndetImages[i]) != R)
        CoCoA_THROW_ERROR(ERR::BadPolyRingHomImages, FnName);

    return Rx->myHomCtor(R, IdentityHom(R), IndetImages);
  }

  RingHom EvalHom(const PolyRing& Rx, const MachineInt& n) // Maps f in R[x] into f(n) in R
  {
    if (NumIndets(Rx) != 1) CoCoA_THROW_ERROR(ERR::BadArg, "EvalHom(Rx,n)");
    const ring& R = CoeffRing(Rx);
    const vector<RingElem> IndetImage(1, RingElem(R,n));
    return Rx->myHomCtor(R, IdentityHom(R), IndetImage);
  }

  RingHom EvalHom(const PolyRing& Rx, const BigInt& N)     // Maps f in R[x] into f(N) in R
  {
    if (NumIndets(Rx) != 1) CoCoA_THROW_ERROR(ERR::BadArg, "EvalHom(Rx,N)");
    const ring& R = CoeffRing(Rx);
    const vector<RingElem> IndetImage(1, RingElem(R,N));
    return Rx->myHomCtor(R, IdentityHom(R), IndetImage);
  }

  RingHom EvalHom(const PolyRing& Rx, const BigRat& q)     // Maps f in R[x] into f(q) in R
  {
    if (NumIndets(Rx) != 1) CoCoA_THROW_ERROR(ERR::BadArg, "EvalHom(Rx,N)");
    const ring& R = CoeffRing(Rx);
    const vector<RingElem> IndetImage(1, RingElem(R,q));
    return Rx->myHomCtor(R, IdentityHom(R), IndetImage);
  }

  RingHom EvalHom(const PolyRing& Rx, ConstRefRingElem r)  // Maps f in R[x] into f(r) in R
  {
    if (NumIndets(Rx) != 1) CoCoA_THROW_ERROR(ERR::BadArg, "EvalHom(Rx,r)");
    const ring& R = CoeffRing(Rx);
    if (owner(r) != R) CoCoA_THROW_ERROR(ERR::MixedRings, "EvalHom(Rx,r");
    const vector<RingElem> IndetImage(1, r);
    return Rx->myHomCtor(R, IdentityHom(R), IndetImage);
  }


  RingHom PolyAlgebraHom(const PolyRing& Rx, const ring& Ry, const std::vector<RingElem>& IndetImages)
  {
    const char* const FnName = "PolyAlgebraHom(Rx,Ry,IndetImages)";
    // Check that IndetImages are sensible...
    if (NumIndets(Rx) != len(IndetImages))
      CoCoA_THROW_ERROR(ERR::BadArraySize, FnName);
    for (long i=0; i < NumIndets(Rx); ++i)
      if (owner(IndetImages[i]) != Ry)
        CoCoA_THROW_ERROR(ERR::BadPolyRingHomImages, FnName);
//     // Special case: codomain is coeff ring.
//     if (Ry == CoeffRing(Rx))
//       return Rx->myHomCtor(Ry, IdentityHom(Ry), IndetImages);
//     // General case: codomain must be a poly ring with same coeffs
//     if (!IsPolyRing(Ry))
//       CoCoA_THROW_ERROR(ERR::BadCodomain, FnName);
//     if (CoeffRing(Rx) != CoeffRing(Ry))
//       CoCoA_THROW_ERROR(ERR::MixedCoeffRings, FnName);
//    return Rx->myHomCtor(Ry, CoeffEmbeddingHom(Ry), IndetImages);
    return Rx->myHomCtor(Ry, CanonicalHom(CoeffRing(Rx),Ry), IndetImages);
  }


  RingHom PolyAlgebraHom(const PolyRing& Rx, const ring& Ry, const std::string& IndetImages)
  { return PolyAlgebraHom(Rx, Ry, RingElems(Ry, IndetImages)); }
  

} // end of namespace CoCoA

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/PolyRing.C,v 1.90 2024/07/22 15:43:58 abbott Exp $
// $Log: PolyRing.C,v $
// Revision 1.90  2024/07/22 15:43:58  abbott
// Summary: Renamed macros CoCoA_THROW_ERROR_WITH_CONTEXT so that last char is num of args
//
// Revision 1.89  2024/07/03 15:35:47  bigatti
// Summary: Specialized some ERR:BadArg with explicit error codes
//
// Revision 1.88  2024/07/03 09:35:15  bigatti
// Summary: error codes (redmine 92):
// removed BadIntedIndex, BadDegIndex, BadComptIndex
//
// Revision 1.87  2024/07/02 20:42:07  abbott
// Summary: Renaming of errors (redmine 308)
//
// Revision 1.86  2024/07/02 15:36:14  bigatti
// Summary:Changed error codes: LogZero into ReqNonZero
// and Not... into ReqNonNegative, ReqNonNegativeGrading, ReqPositive, ReqPositiveGrading
//
// Revision 1.85  2024/07/02 09:57:57  bigatti
// Summary: error codes, first batch:
// ReqUnivariate, ReqNonZero, ReqNonZeroGradingDim, ReqNonZeroModulus,  ReqNonZeroRingElem
//
// Revision 1.84  2024/03/26 13:57:40  bigatti
// Summary: applied ErrorContext to NumTerms (Redmine #582)
//
// Revision 1.83  2024/03/16 10:15:12  abbott
// Summary: Use ErrorContext when calling arg-check fns (redmine 1793)
//
// Revision 1.82  2024/03/15 18:57:32  abbott
// Summary: Improved use of PolyRing inside fns; employed CoCoA_ERROR_CONTEXT
//
// Revision 1.81  2023/12/09 20:29:49  abbott
// Summary: Split "eval" fns out of SparsePolyOps-cyclotomic, into SparsePolyOps-eval; added example & doc
//
// Revision 1.80  2023/12/01 20:27:09  abbott
// Summary: Improved spacing
//
// Revision 1.79  2023/11/21 22:08:31  abbott
// Summary: Partial update to FixedDivisor (using new TEMPORARY EvalUPoly class)
//
// Revision 1.78  2023/11/20 20:44:56  abbott
// Summary: Made loop interruptible; improved indentation
//
// Revision 1.77  2023/03/08 19:43:56  abbott
// Summary: Made myCharacteristic a function (instead of proc with OUT param) redmine 536
//
// Revision 1.76  2022/11/16 19:44:39  abbott
// Summary: Fixed bug when coeffs are ZZ: redmine 1710
//
// Revision 1.75  2022/10/11 15:32:30  bigatti
// Summary: using new interface of exponents returning vector
//
// Revision 1.74  2022/06/08 15:22:36  abbott
// Summary: IsSqFree now works over ZZ (redmine 1683)
//
// Revision 1.73  2022/02/18 14:11:56  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.72  2021/03/01 19:40:24  abbott
// Summary: Moved defns of monic, CommonDenom, ClearDenom to new file PolyRing-content.C
//
// Revision 1.71  2020/06/17 15:49:25  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.70  2020/02/11 16:56:41  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.69  2020/02/11 16:12:18  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.68  2020/01/26 14:41:59  abbott
// Summary: Revised includes after splitting NumTheory (redmine 1161)
//
// Revision 1.67  2020/01/08 15:33:37  abbott
// Summary: Added check for empty list to CommonDenom (redmine 1382)
//
// Revision 1.66  2019/10/18 14:16:01  bigatti
// -- renamed ReadExprVector --> RingElems
//
// Revision 1.65  2019/10/10 06:46:25  bigatti
// -- added PolyRingHom PolyAlgebraHom with string (using ReadExprVector)
//
// Revision 1.64  2019/09/16 17:22:27  abbott
// Summary: Added IsEvenPoly, IsOddPoly and FixedDivisor
//
// Revision 1.63  2019/03/19 11:07:07  abbott
// Summary: Replaced 0 by nullptr where appropriate
//
// Revision 1.62  2018/06/15 08:47:24  abbott
// Summary: Replaced IsOne(gcd...)) by IsCoprime(...)
//
// Revision 1.61  2018/05/18 16:42:11  bigatti
// -- added include SparsePolyOps-RingElem.H
//
// Revision 1.60  2018/05/17 15:39:19  bigatti
// -- renamed VectorOperations --> VectorOps
// -- sorted includes
//
// Revision 1.59  2018/04/20 18:51:25  abbott
// Summary: Changed ctors for BigInt/BigRat from string or from MPZ/MPQ
//
// Revision 1.58  2018/03/20 11:40:43  bigatti
// -- added CommonDenom for vector of polynomials over QQ
//
// Revision 1.57  2017/04/21 18:07:53  bigatti
// -- improved error messages for PolyRingHom
//
// Revision 1.56  2017/04/21 15:46:23  bigatti
// -- improved error messages for PolyRingHom
//
// Revision 1.55  2016/11/03 12:25:25  abbott
// Summary: Changed IsRadical (for PPMonoidElem) into IsSqFree
//
// Revision 1.54  2016/10/25 20:54:09  abbott
// Summary: Added new fn IsSqFree (for BigInt and ringelem of PolyRing over field)
//
// Revision 1.53  2015/11/04 10:45:02  abbott
// Summary: Consequential change: changed arg order to RegisterDtorForGlobal
//
// Revision 1.52  2015/07/29 13:50:45  bigatti
// -- renamed ID into RingID
//
// Revision 1.51  2015/07/29 11:04:55  bigatti
// -- added space after comma in printing rings
//
// Revision 1.50  2015/03/04 11:22:52  bigatti
// -- changed ERR::BadPolyRingHomImages ERR::BadArraySize
//
// Revision 1.49  2014/07/31 14:45:18  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.48  2014/07/28 15:44:30  abbott
// Summary: Added defn of CoeffEmbeddingHom (was inline)
// Author: JAA
//
// Revision 1.47  2014/07/14 15:07:42  abbott
// Summary: Changed include of tmp.H into UtilsTemplate.H
// Author: JAA
//
// Revision 1.46  2014/07/14 11:39:32  abbott
// Summary: Added convenient "univariate" EvalHom signatures
// Author: JAA
//
// Revision 1.45  2014/07/11 15:44:42  bigatti
// -- default impl for myOutputSelf
//
// Revision 1.44  2014/07/09 11:42:11  abbott
// Summary: Removed AsPolyRing from commented out code
// Author: JAA
//
// Revision 1.43  2014/07/08 08:36:07  abbott
// Summary: Removed AsFractionField
// Author: JAA
//
// Revision 1.42  2014/07/07 17:11:20  abbott
// Summary: [MAJOR CHANGE] Removed AsPolyRing; added PolyRingPtr
// Author: JAA
//
// Revision 1.41  2014/05/06 13:20:41  abbott
// Summary: Changed names (my)MaxExponents into (my)Deg
// Author: JAA
//
// Revision 1.40  2014/04/11 15:07:21  abbott
// Summary: Renamed TmpFactor to factor in include
// Author: JAA
//
// Revision 1.39  2014/03/14 11:08:20  abbott
// Summary: Changed non-standard errors into std ones (with code ERR::BadQuot)
// Author: JAA
//
// Revision 1.38  2013/06/18 12:27:37  abbott
// Renamed HibertSeriesPolyRing to RingQQt.
//
// Revision 1.37  2013/05/20 16:20:33  abbott
// Improved arg sanity checks to deriv.
//
// Revision 1.36  2013/05/14 14:21:31  abbott
// Revised/improved impl of derivative of ratfns.
//
// Revision 1.35  2013/03/26 14:58:08  abbott
// Added IndetPower allowing exp to be a BigInt (handy elsewhere).
//
// Revision 1.34  2013/01/17 13:46:00  abbott
// Changed name:  IndetsCalled -->  indets
// Added virt mem fn  myImageLiesInSubfield
//
// Revision 1.33  2012/10/24 13:34:42  abbott
// Changed ConstRefRingElem into RingElemAlias in ctor calls.
// Removed use of std::ptr_fun  (used explicit loop instead)
//
// Revision 1.32  2012/05/22 10:02:37  abbott
// Removed IsGCDDomain; substituted by IsTrueGCDDomain.
// Added IsFractionFieldOfGCDDomain.
//
// Revision 1.31  2012/05/20 09:45:10  abbott
// Corrected error code (to NotTrueGCDDomain) in content.
//
// Revision 1.30  2012/04/04 13:29:50  abbott
// Removed more commented out old code.
//
// Revision 1.29  2012/04/04 13:28:35  abbott
// Removed old commented out code.
//
// Revision 1.28  2012/04/04 11:03:25  bigatti
// -- fixed content for myContentFrF: ans is in FrF, not in BaseRing(FrF)
//
// Revision 1.27  2012/01/30 23:26:25  abbott
// Added check for monic imput poly in myMonic.
//
// Revision 1.26  2011/11/07 10:30:30  bigatti
// -- added AreMonomials
//
// Revision 1.25  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.24  2011/06/27 13:43:57  bigatti
// -- allowing CoeffRing(Rx) -> CoeffRing(S) as CoeffHom in PolyRingHom
//
// Revision 1.23  2011/06/23 16:04:47  abbott
// Added IamExact mem fn for rings.
// Added myRecvTwinFloat mem fn for rings.
// Added first imple of RingHom from RingTwinFloat to other rings.
//
// Revision 1.22  2011/05/20 16:10:58  bigatti
// -- small extension to PolyAlgebraHom
//
// Revision 1.21  2011/03/16 15:40:22  bigatti
// -- added myIsIndetPosPower(f), IsIndetPosPower(f)
//
// Revision 1.20  2011/03/16 13:17:21  abbott
// Removed indet(PolyRing,ZZ).
// Added commented out IndetPower with ZZ exponent.
//
// Revision 1.19  2011/03/10 16:39:34  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.18  2011/03/01 14:10:47  bigatti
// -- added ClearDenom/myClearDenom
//
// Revision 1.17  2011/01/24 17:16:25  bigatti
// -- added monic/myMonic
//
// Revision 1.16  2011/01/18 14:38:43  bigatti
// -- moved **_forC5 functions into CoCoA-5/CoCoALibSupplement:
//    myMonomials_forC5, mySupport_forC5, monomials_forC5, support_forC5,
//    LPP_forC5, LT_forC5, LM_forC5
//
// Revision 1.15  2010/11/30 11:19:50  bigatti
// -- unique implementation of IndetsCalled
// -- added virtual myIndetSymbol
//
// Revision 1.14  2010/11/25 09:33:04  bigatti
// -- sorted includes (cstddef and vector)
//
// Revision 1.13  2010/11/05 16:03:50  bigatti
// -- added #include <ZZ> for CoCoA_ASSERT
//
// Revision 1.12  2010/11/05 15:56:13  bigatti
// -- added myMonomials_forC5, mySupport_forC5
// -- added monomials_forC5, support_forC5
// -- added LPP_forC5, LM_forC5
//
// Revision 1.11  2010/11/02 16:03:12  bigatti
// -- added indet(*, ZZ) [for CoCoA-5]
// -- indet(PPM, i) now returns "const PPMonoidElem&"
//
// Revision 1.10  2009/10/02 13:27:26  bigatti
// -- unique implementation of myDiv in PolyRing.C
//
// Revision 1.9  2009/07/24 12:26:43  abbott
// Added CommonDenom function for polynomials.
//
// Revision 1.8  2008/04/21 12:50:09  abbott
// Added EvalHom, an easy way to create evaluation homomorphisms for poly rings.
//
// Revision 1.7  2008/03/12 16:02:47  bigatti
// -- added: IsIrred
//
// Revision 1.6  2007/12/05 11:06:24  bigatti
// -- changed "size_t StdDeg/myStdDeg(f)" into "long"  (and related functions)
// -- changed "log/myLog(f, i)" into "MaxExponent/myMaxExponent(f, i)"
// -- fixed bug in "IsOne(ideal)" in SparsePolyRing.C
//
// Revision 1.5  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.4  2007/09/24 14:24:08  abbott
// Changed content to work if input is zero or if coeffs are in a field.
// Cleaned impl of PolyAlgebraHom to keep gcc-4.3 prerelease quiet.
//
// Revision 1.3  2007/05/31 15:26:23  bigatti
// -- default implementation for IamCommutative, IamIntegralDomain,
//    IamGCDDomain, IamField, myCharacteristic
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.7  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.6  2007/03/07 13:59:45  bigatti
// -- added PolyAlgebraHom(Rx, Ry, IndetImages)
//
// Revision 1.5  2007/03/05 21:06:07  cocoa
// New names for homomorphism pseudo-ctors: removed the "New" prefix.
//
// Revision 1.4  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.3  2007/02/28 13:51:59  bigatti
// -- added function IsMonomial
//
// Revision 1.2  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.2  2006/02/20 22:41:19  cocoa
// All forms of the log function for power products now return SmallExponent_t
// (instead of int).  exponents now resizes the vector rather than requiring
// the user to pass in the correct size.
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.3  2005/08/08 16:36:32  cocoa
// Just checking in before going on holiday.
// Don't really recall what changes have been made.
// Added IsIndet function for RingElem, PPMonoidElem,
// and a member function of OrdvArith.
// Improved the way failed assertions are handled.
//
// Revision 1.2  2005/07/01 16:08:15  cocoa
// Friday check-in.  Major change to structure under PolyRing:
// now SparsePolyRing and DUPolyRing are separated (in preparation
// for implementing iterators).
//
// A number of other relatively minor changes had to be chased through
// (e.g. IndetPower).
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.7  2005/04/29 15:42:02  cocoa
// Improved documentation for GMPAllocator.
// Added example program for GMPAllocator.
// Added example program for simple ops on polynomials.
// Added two new ctors for (principal) ideals (from long, and from ZZ).
// Added (crude) printing for PPMonoids.
// Updated library.H (#included GMPAllocator.H).
//
// Revision 1.6  2005/04/27 16:14:56  cocoa
// Cleaned up example programs -- added "free use" permit.
// Changed a couple of ErrorInfo object names, and added
// ERR::NotTrueGCDDomain.
//
// Revision 1.5  2005/04/20 15:40:48  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.4  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.3  2005/03/11 18:36:36  cocoa
// -- new: NewPolyRingHom Rx-->Ry ; code is still commented out
//
// Revision 1.2  2005/02/11 14:15:20  cocoa
// New style ring elements and references to ring elements;
// I hope I have finally got it right!
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.10  2004/12/16 13:35:01  cocoa
// Added new monomial functions (with more coherent signatures).
//
// Revision 1.9  2004/11/25 16:14:21  cocoa
// (1) Fixed definition of specialization of std::swap template function
//     so that it compiles with gcc 3.4.3
// (2) Implemented monomial function for polynomial rings.
// (3) Added one(PPM) and PPM->myOne() functions.
//
// Revision 1.8  2004/11/19 16:15:51  cocoa
// (a) Removed unused error message about degree of zero;
//     replaced it by a more general message about needing a
//     non-zero polynomial (for various ops such as LC, LPP).
// (b) Added some missing arg checking in LC, LPP and deg
//     (for elems of a PolyRing).
// (c) Updated some commented out code in GPair and GPoly.
//
// Revision 1.7  2004/11/19 15:44:27  cocoa
// Changed names of "casting" functions which convert a ring into
// one with a more special structure (e.g. FractionField).  These
// functions now have names starting with "As".  There were several
// consequential changes.
//
// Revision 1.6  2004/11/18 18:33:41  cocoa
// Now every ring know its own "one" element (as well as "zero").
// Several consequential changes.
//
// Revision 1.5  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.4  2004/10/29 16:05:03  cocoa
// -- changed PPOrdering::ExpvElem --> SmallExponent_t
//
// Revision 1.3  2004/05/24 15:52:14  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.2  2003/10/09 12:18:42  cocoa
// New coding convention for rings.
//
// Revision 1.5  2003/06/23 17:04:45  abbott
// Minor cleaning prior to public release.
//
// Revision 1.4  2003/05/30 12:25:12  abbott
// Corrected signature of IndetPower: previously allowed signed exponent.
//
// Revision 1.3  2003/05/14 16:53:05  abbott
// Numerous changes reflecting the renaming and restructuring
// mentioned in PolyRing.H (qv.).
//
// Revision 1.2  2002/11/15 15:15:55  abbott
// Revised according to the renaming in ring.H.
// Error messages improved.
//
// Revision 1.1  2002/06/22 16:54:40  abbott
// Initial revision
//
//

