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

  // void PolyRingBase::myCheckIndetIndex(long i, const ErrorContext& where) const
  // {
  //   if (i < 0 || i >= myNumIndets())
  //     CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::BadIndetIndex, where);
  // }

  void PolyRingBase::myCheckIndetIndex(long i, const ErrorContext& ErrCtx) const
  {
    if (i < 0 || i >= myNumIndets())
      CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::BadIndex, ErrCtx);
  }

  bool PolyRingBase::myIsIrred(ConstRawPtr rawf) const
  {
    if (myIsZero(rawf))  CoCoA_THROW_ERROR1(ERR::ReqNonZeroRingElem);
    if (myIsInvertible(rawf))  CoCoA_THROW_ERROR1(ERR::InvertibleRingElem);
    return IsIrredPoly(RingElemAlias(ring(this), rawf));
  }


  void PolyRingBase::myDiv(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    if (myIsZero(rawy)) // or CoCoA_ASSERT???
      CoCoA_THROW_ERROR1(ERR::DivByZero);
    if (!myIsDivisible(rawlhs, rawx, rawy))
      CoCoA_THROW_ERROR1(ERR::BadQuot);
  }  
  

  void PolyRingBase::myMonic(RawPtr rawmonic, ConstRawPtr rawf) const
  {
    if (myIsZero(rawf)) // or CoCoA_ASSERT???
      CoCoA_THROW_ERROR1(ERR::ReqNonZeroRingElem);
    RingElem ans = RingElemAlias(ring(this), rawf);
    if (!IsOne(myLC(rawf)) && !myDivByCoeff(raw(ans), raw(myLC(rawf))))
      CoCoA_THROW_ERROR1(ERR::BadQuot);
    mySwap(rawmonic, raw(ans));
  }
  

  bool IsSqFree(ConstRefRingElem f)
  {
    const ring& P = owner(f);
    if (!IsPolyRing(P))  CoCoA_THROW_ERROR1(ERR::ReqPolyRing);
    if (!IsZZ(CoeffRing(P)) && !IsField(CoeffRing(P)))
      CoCoA_THROW_ERROR1("CoeffRing must be ZZ or a field");
    if (IsZero(f))  CoCoA_THROW_ERROR1(ERR::ReqNonZero);
    if (IsConstant(f))  return true;

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
    CoCoA_THROW_ERROR1(ERR::NYI);
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
    if (IsNegative(NumIndets) || IsZero(NumIndets))  CoCoA_THROW_ERROR1(ERR::ReqPositive);
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
    if (exp < 0)  CoCoA_THROW_ERROR1(ERR::NegExp);
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
  //     CoCoA_THROW_ERROR1(ERR::NotElemPolyRing);
  //   return PolyRingPtr(owner(f))->myNumTerms(raw(f));
  // }

  long NumTerms(ConstRefRingElem f)
  { return PolyRingPtr(owner(f), CoCoA_ERROR_CONTEXT)->myNumTerms(raw(f)); }


  bool IsMonomial(ConstRefRingElem f)
  {
    if (!IsPolyRing(owner(f)))
      CoCoA_THROW_ERROR1(ERR::ReqElemPolyRing);
    return PolyRingPtr(owner(f))->myIsMonomial(raw(f));
  }


  bool AreMonomials(const std::vector<RingElem>& v)
  {
    // morally:  return find_if(v.begin(), v.end(), not1(IsMonomial)) == v.end();
    if (!HasUniqueOwner(v))  CoCoA_THROW_ERROR1(ERR::MixedRings);
    const long n = len(v);
    for (long i=0; i < n; ++i)
      if (!IsMonomial(v[i]))  return false;
    return true;
//  We *DO NOT USE* STL algorithm because ptr_fun fails when args are references.
//     return find_if(v.begin(), v.end(),
//                    not1(ptr_fun(static_cast<bool(*)(const RingElemAlias&)>(CoCoA::IsMonomial))))
//       == v.end(); 
  }


  bool IsConstant(ConstRefRingElem f)
  {
    if (!IsPolyRing(owner(f)))
      CoCoA_THROW_ERROR1(ERR::ReqElemPolyRing);
    return PolyRingPtr(owner(f))->myIsConstant(raw(f));
  }


  bool IsIndet(ConstRefRingElem f)
  {
    if (!IsPolyRing(owner(f)))
      CoCoA_THROW_ERROR1(ERR::ReqElemPolyRing);
    long junk;
    return PolyRingPtr(owner(f))->myIsIndet(junk, raw(f));
  }


  bool IsIndet(long& index, ConstRefRingElem f)
  {
    if (!IsPolyRing(owner(f)))
      CoCoA_THROW_ERROR1(ERR::ReqElemPolyRing);
    return PolyRingPtr(owner(f))->myIsIndet(index, raw(f));
  }


  bool IsIndetPosPower(ConstRefRingElem f)
  {
    if (!IsPolyRing(owner(f)))
      CoCoA_THROW_ERROR1(ERR::ReqElemPolyRing);
    return PolyRingPtr(owner(f))->myIsIndetPosPower(raw(f));
  }


  bool IsEvenPoly(ConstRefRingElem f)
  {
    if (!IsPolyRing(owner(f)))
      CoCoA_THROW_ERROR1(ERR::ReqElemPolyRing);
    return PolyRingPtr(owner(f))->myIsEvenPoly(raw(f));
  }
  
  bool IsOddPoly(ConstRefRingElem f)
  {
    if (!IsPolyRing(owner(f)))
      CoCoA_THROW_ERROR1(ERR::ReqElemPolyRing);
    return PolyRingPtr(owner(f))->myIsOddPoly(raw(f));
  }


  long deg(ConstRefRingElem f, long var)
  {
    if (!IsPolyRing(owner(f)))
      CoCoA_THROW_ERROR1(ERR::ReqElemPolyRing);
    const PolyRingBase* P = PolyRingPtr(owner(f));
    P->myCheckIndetIndex(var, CoCoA_ERROR_CONTEXT);
    if (IsZero(f))
      CoCoA_THROW_ERROR1(ERR::ReqNonZeroRingElem);
    return P->myDeg(raw(f), var);
  }


  // This impl allows content to be greater than 1 (but requires integer coeffs)
  RingElem FixedDivisor(ConstRefRingElem f)
  {
    // hint: could be clever if f is even or odd (need do only half as many evals)
    const PolyRing P(owner(f), CoCoA_ERROR_CONTEXT);
//    if (!IsPolyRing(P))  CoCoA_THROW_ERROR1(ERR::ReqPolyRing);
    if (!IsZero(characteristic(P)))  CoCoA_THROW_ERROR1(ERR::ReqChar0);
    if (IsZero(f))  CoCoA_THROW_ERROR1(ERR::ReqNonZeroRingElem);
    if (deg(f) == 0)  return abs(LC(f));
    const long x_ind = UnivariateIndetIndex(f);
    if (x_ind < 0)  CoCoA_THROW_ERROR1(ERR::ReqUnivariate);
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
      CheckForInterrupt("FixedDivisor");
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
      if (!IsPolyRing(BaseRing(FrF)))  CoCoA_THROW_ERROR1(ERR::ReqPolyRing);
      if (!IsOne(den(x)))  CoCoA_THROW_ERROR1("Nonsensical derivative");  // BUG ??? also check that num(x) IsIndet ????
      RingElem ans(FrF);
      FrF->myDeriv(raw(ans), raw(f), raw(x));
      return ans;
    }
  } // end of anonymous namespace

  RingElem deriv(ConstRefRingElem f, ConstRefRingElem x)
  {
    if (owner(x) != owner(f))  CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (IsFractionField(owner(f))) return DerivFrF(f,x);
    // From here on we are in the "polynomial" case.
    if (!IsIndet(x))  CoCoA_THROW_ERROR1(ERR::ReqIndet);
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


  void MakeMonic(std::vector<RingElem>& v)
  {
    for (RingElem& g: v)  g = monic(g);  // might throw
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
      CoCoA_THROW_ERROR2(ERR::BadDomain, "CoeffHom");
    if (IsPolyRing(S) && codomain(CoeffHom) == CoeffRing(S))
      CoeffHom = CoeffEmbeddingHom(S)(CoeffHom);
    if (codomain(CoeffHom) != S)
      CoCoA_THROW_ERROR2(ERR::BadCodomain, "CoeffHom");
    if (NumIndets(Rx) != len(IndetImages))
      CoCoA_THROW_ERROR1(ERR::BadArraySize);
    for (long i=0; i < NumIndets(Rx); ++i)
      if (owner(IndetImages[i]) != S)
        CoCoA_THROW_ERROR1(ERR::BadPolyRingHomImages);

    return Rx->myHomCtor(S, CoeffHom, IndetImages);
  }


  // Rx is the domain, S is the codomain
  RingHom PolyRingHom(const PolyRing& Rx, const ring& S, RingHom CoeffHom, const std::string& IndetImages)
  { return PolyRingHom(Rx, S, CoeffHom, RingElems(S, IndetImages)); }
  

  RingHom EvalHom(const PolyRing& Rx, const std::vector<RingElem>& IndetImages)
  {
    const ring& R = CoeffRing(Rx);
    if (NumIndets(Rx) != len(IndetImages))
      CoCoA_THROW_ERROR1(ERR::BadArraySize);
    for (long i=0; i < NumIndets(Rx); ++i)
      if (owner(IndetImages[i]) != R)
        CoCoA_THROW_ERROR1(ERR::BadPolyRingHomImages);

    return Rx->myHomCtor(R, IdentityHom(R), IndetImages);
  }

  RingHom EvalHom(const PolyRing& Rx, const MachineInt& n) // Maps f in R[x] into f(n) in R
  {
    if (NumIndets(Rx) != 1)  CoCoA_THROW_ERROR1(ERR::ReqUnivariate);
    const ring& R = CoeffRing(Rx);
    const vector<RingElem> IndetImage(1, RingElem(R,n));
    return Rx->myHomCtor(R, IdentityHom(R), IndetImage);
  }

  RingHom EvalHom(const PolyRing& Rx, const BigInt& N)     // Maps f in R[x] into f(N) in R
  {
    if (NumIndets(Rx) != 1)  CoCoA_THROW_ERROR1(ERR::ReqUnivariate);
    const ring& R = CoeffRing(Rx);
    const vector<RingElem> IndetImage(1, RingElem(R,N));
    return Rx->myHomCtor(R, IdentityHom(R), IndetImage);
  }

  RingHom EvalHom(const PolyRing& Rx, const BigRat& q)     // Maps f in R[x] into f(q) in R
  {
    if (NumIndets(Rx) != 1)  CoCoA_THROW_ERROR1(ERR::ReqUnivariate);
    const ring& R = CoeffRing(Rx);
    const vector<RingElem> IndetImage(1, RingElem(R,q));
    return Rx->myHomCtor(R, IdentityHom(R), IndetImage);
  }

  RingHom EvalHom(const PolyRing& Rx, ConstRefRingElem r)  // Maps f in R[x] into f(r) in R
  {
    if (NumIndets(Rx) != 1)  CoCoA_THROW_ERROR1(ERR::ReqUnivariate);
    const ring& R = CoeffRing(Rx);
    if (owner(r) != R)  CoCoA_THROW_ERROR1(ERR::MixedRings);
    const vector<RingElem> IndetImage(1, r);
    return Rx->myHomCtor(R, IdentityHom(R), IndetImage);
  }


  RingHom PolyAlgebraHom(const PolyRing& Rx, const ring& Ry, const std::vector<RingElem>& IndetImages)
  {
    // Check that IndetImages are sensible...
    if (NumIndets(Rx) != len(IndetImages))
      CoCoA_THROW_ERROR1(ERR::BadArraySize);
    for (long i=0; i < NumIndets(Rx); ++i)
      if (owner(IndetImages[i]) != Ry)
        CoCoA_THROW_ERROR1(ERR::BadPolyRingHomImages);
//     // Special case: codomain is coeff ring.
//     if (Ry == CoeffRing(Rx))
//       return Rx->myHomCtor(Ry, IdentityHom(Ry), IndetImages);
//     // General case: codomain must be a poly ring with same coeffs
//     if (!IsPolyRing(Ry))
//       CoCoA_THROW_ERROR1(ERR::BadCodomain);
//     if (CoeffRing(Rx) != CoeffRing(Ry))
//       CoCoA_THROW_ERROR1(ERR::MixedCoeffRings);
//    return Rx->myHomCtor(Ry, CoeffEmbeddingHom(Ry), IndetImages);
    return Rx->myHomCtor(Ry, CanonicalHom(CoeffRing(Rx),Ry), IndetImages);
  }


  RingHom PolyAlgebraHom(const PolyRing& Rx, const ring& Ry, const std::string& IndetImages)
  { return PolyAlgebraHom(Rx, Ry, RingElems(Ry, IndetImages)); }
  

} // end of namespace CoCoA
