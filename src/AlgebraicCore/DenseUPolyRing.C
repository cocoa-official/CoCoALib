//   Copyright (c)  2007,2009 John Abbott and Anna Bigatti
//   Author:  2007,2009  Anna Bigatti

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

#include "CoCoA/DenseUPolyRing.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/CpuTimeLimit.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/QuotientRing.H"  // for IsQuotientRing
#include "CoCoA/RingDenseUPolyClean.H" // for NewPolyRing_DUP
#include "CoCoA/RingFp.H" // for IsRingFp
#include "CoCoA/RingQQ.H" // for IsQQ
#include "CoCoA/RingTwinFloat.H" // for IsRingTwinFloat
#include "CoCoA/RingZZ.H" // for IsZZ
#include "CoCoA/assert.H"
#include "CoCoA/degree.H"
#include "CoCoA/error.H"
#include "CoCoA/symbol.H"
#include "CoCoA/utils.H"

//#include <vector>
using std::vector;
#include <algorithm>
using std::max;
using std::sort;


namespace CoCoA
{

  DenseUPolyRing NewPolyRing(const ring& CoeffRing)
  {
    return NewPolyRing_DUP(CoeffRing);
  }

  DenseUPolyRing NewPolyRing(const ring& CoeffRing, const symbol& IndetSym)
  {
    return NewPolyRing_DUP(CoeffRing, IndetSym);
  }

  bool IsGoodIndetName(const ring& CoeffRing, const symbol& IndetName)
  {
    vector<symbol> syms = symbols(CoeffRing);
    const long NumSyms = len(syms);
    for (long i=0; i < NumSyms; ++i)
    {
      if (syms[i] == IndetName) return false;
      if (head(syms[i]) == head(IndetName) &&
          NumSubscripts(syms[i]) != NumSubscripts(IndetName))
        return false;
    }
    return true;
  }
  

  RingElem DenseUPolyRingBase::mySymbolValue(const symbol& s) const
  {
    if (s == myIndetSymbol()) return myIndets()[0];
    return myCoeffEmbeddingHomCtor()(myCoeffRing()->mySymbolValue(s));
  }


  //---- Functions for creating/building polynomials

  RingElem monomial(const DenseUPolyRing& P, ConstRefRingElem c, const MachineInt& exp)
  {
    if (owner(c) != CoeffRing(P)) CoCoA_THROW_ERROR(ERR::MixedCoeffRings, "monomial(P,c,d)");
    if (IsNegative(exp))
      CoCoA_THROW_ERROR(ERR::NegExp, "monomial(P,c,d)");
    if (!IsSignedLong(exp))
      CoCoA_THROW_ERROR(ERR::ExpTooBig, "monomial(P,c,d)");
    if (IsZero(c)) return zero(P);
    return P->myMonomial(raw(c), AsSignedLong(exp));
  }
  

  RingElem monomial(const DenseUPolyRing& P, const BigInt& N, const MachineInt& exp)
  {
    // All arg checking made by the call below...
    return monomial(P, RingElem(CoeffRing(P), N), exp);
  }
  

  RingElem monomial(const DenseUPolyRing& P, const BigRat& N, const MachineInt& exp)
  {
    // All arg checking made by the call below...
    return monomial(P, RingElem(CoeffRing(P), N), exp);
  }
  

  RingElem monomial(const DenseUPolyRing& P, const MachineInt& n, const MachineInt& exp)
  {
    // All arg checking made by the call below...
    return monomial(P, RingElem(CoeffRing(P), n), exp);
  }


  /*----------------------------------------------------------------------
    Member functions every concrete DenseUPolyRing implementation
    must have in addition to those of PolyRingBase.
    ----------------------------------------------------------------------*/

  bool DenseUPolyRingBase::myIsValid(ConstRawPtr rawf) const
  {
    if (myDegPlus1(rawf)>mySize(rawf)) return false;
    for (long i=myDegPlus1(rawf); i<mySize(rawf); ++i)
      if (!IsZero(myCoeff(rawf, i))) return false;
    if (myIsZero(rawf)) return true;
    if (IsZero(myLC(rawf))) return false;
    return true;
  }


  long DenseUPolyRingBase::myStdDeg(ConstRawPtr rawf) const
  {
    if (myIsZero(rawf)) CoCoA_THROW_ERROR(ERR::ReqNonZeroRingElem, "myStdDeg(rawf)");
    return myDegPlus1(rawf) - 1;
  }


  long DenseUPolyRingBase::myDeg(ConstRawPtr rawf, long index) const
  {
    if (index!=0) CoCoA_THROW_ERROR(ERR::BadIndex, "myDeg(rawf, index)");
    return myStdDeg(rawf);
  }


  RingElemAlias DenseUPolyRingBase::myLC(ConstRawPtr rawf) const
  {
    if (myIsZero(rawf)) CoCoA_THROW_ERROR(ERR::ReqNonZeroRingElem, "DenseUPolyRingBase::myLC(rawf)");
    return myCoeff(rawf, myStdDeg(rawf));
  }
  

  void DenseUPolyRingBase::myContent(RawPtr rawcontent, ConstRawPtr rawf) const
  {
    const ring& R = myCoeffRing();
    CoCoA_ASSERT(IsTrueGCDDomain(R));
    // Compute answer in local var to avoid aliasing problems; also exception clean.
    RingElem ans(R);
    for (long i=0; i<myDegPlus1(rawf); ++i)
    {
      if (IsZero(myCoeff(rawf,i))) continue;
      R->myGcd(raw(ans), raw(ans), raw(myCoeff(rawf, i))); // ans = GCD(ans, coeff(i));
      if (IsOne(ans)) break;
    }
    // Finally, swap answer into rawcontent.
    R->mySwap(rawcontent, raw(ans));
  }

  void DenseUPolyRingBase::myContentFrF(RawPtr rawcontent, ConstRawPtr rawf) const
  {
    CoCoA_ASSERT(IsFractionFieldOfGCDDomain(myCoeffRing()));
    const FractionField R = myCoeffRing();
    const ring& S = BaseRing(R);
    RingElem N(S);
    RingElem D(S,1);
    for (long i=0; i < myDegPlus1(rawf); ++i)
    {
      if (IsZero(myCoeff(rawf,i))) continue;
//      R->myGcd(raw(ans), raw(ans), raw(num(myCoeff(rawf, i)))); // ans = GCD(ans, num(coeff(i)));
      N = gcd(N, num(myCoeff(rawf, i)));
      D = lcm(D, den(myCoeff(rawf, i)));
//      if (IsOne(ans)) break;
    }
    RingHom phi = EmbeddingHom(R);
    RingElem ans = phi(N)/phi(D);
    // Finally, swap answer into rawcontent.
    R->mySwap(rawcontent, raw(ans));
  }


  void DenseUPolyRingBase::myCommonDenom(RawPtr rawcontent, ConstRawPtr rawf) const
  {
    CoCoA_ASSERT(IsFractionFieldOfGCDDomain(myCoeffRing()));
    const ring& R = BaseRing(myCoeffRing());
    // Compute answer in local var to avoid aliasing problems; also exception clean.
    RingElem ans = one(R);
    for (long i=0; i<myDegPlus1(rawf); ++i)
    {
      ConstRefRingElem C = myCoeff(rawf, i);
      if (!IsZero(C))
      {
        const RingElem D = den(C);
        ans *= D/gcd(D, ans);
      }
    }
    // Finally, swap answer into rawcontent.
    R->mySwap(rawcontent, raw(ans));
  }


  void DenseUPolyRingBase::myClearDenom(RawPtr rawg, ConstRawPtr rawf) const
  {
    CoCoA_ASSERT(IsFractionFieldOfGCDDomain(myCoeffRing()));
    const ring& R = BaseRing(myCoeffRing());
    // Compute answer in local var to avoid aliasing problems; also exception clean.
    RingElem c(R);
    myCommonDenom(raw(c), rawf);
    const RingElem coeff = CoeffEmbeddingHom(DenseUPolyRing(this))(EmbeddingHom(myCoeffRing()))(c);
    myMul(rawg, raw(coeff), rawf);
  }


  void DenseUPolyRingBase::myRemoveBigContent(RawPtr rawf) const
  {
    CoCoA_ASSERT(IsTrueGCDDomain(myCoeffRing()));
    CoCoA_ASSERT(!myIsZero(rawf));
    RingElem cont(myCoeffRing());
    myContent(raw(cont), rawf);
    myDivByCoeff(rawf, raw(cont));
  }


  /*----------------------------------------------------------------------
    Member functions inherited from ring with a single implementation
    for all DenseUPolyRing implementations
    ----------------------------------------------------------------------*/

  namespace // anonymous
  {
    // only if CoeffRing is a field
    void mod(RingElem& a, ConstRefRingElem b)  // a = a%b
    {
      CoCoA_ASSERT(owner(a) == owner(b));
      CoCoA_ASSERT(!IsZero(b));  // this should be an error
      while (!IsZero(a) && deg(a) >= deg(b))
        DenseUPolyRingPtr(owner(a))->myAddMulLM(raw(a), raw(-LC(a)/LC(b)), deg(a)-deg(b), raw(b));
    }
    
    void DivMod(RingElem& q, RingElem& a, ConstRefRingElem b)  // q = a/b;  a = a%b
    {
      // pointers must be different
      CoCoA_ASSERT(owner(a) == owner(b));
      CoCoA_ASSERT(owner(q) == owner(a));
      CoCoA_ASSERT(!IsZero(b));  // this should be an error
      const DenseUPolyRing P = owner(a);
      q = zero(owner(q));
      if (IsZero(a) || (deg(a)<deg(b)) ) return;
      P->myResize(raw(q), deg(a)-deg(b)+1);
      while (!IsZero(a) && deg(a) >= deg(b))
      {
        // q += monomial(P, -LC(a)/LC(b), deg(a)-deg(b));
        P->myAssignNonZeroCoeff(raw(q), raw(LC(a)/LC(b)), deg(a)-deg(b));
        P->myAddMulLM(raw(a), raw(-LC(a)/LC(b)), deg(a)-deg(b), raw(b));
      }
    }
    
    void ExgcdAux(RingElem& acofac, RingElem& bcofac, RingElem& a, RingElem& b) // gcd is acofac*orig_a + bcofac*orig_b;  terminates with a=gcd and b=0
    {
      // all pointers must be different
      CoCoA_ASSERT(owner(a) == owner(b));
      CoCoA_ASSERT(!IsZero(b));  // this should be an error
      acofac = 1;
      bcofac = 0;  // i.e. a = acofac*a + bcofac*b
      RingElem q(owner(b));
      RingElem b_acofac(owner(b));
      RingElem b_bcofac(one(owner(b)));  // i.e. b = b_acofac*a + b_bcofac*b
      while (!IsZero(b))
      {
        // q = div(a, b);  a = a - b*d;
        DivMod(q, a, b);
        acofac -= b_acofac*q;
        bcofac -= b_bcofac*q;
        swap(a, b);
        swap(acofac, b_acofac);
        swap(bcofac, b_bcofac);
      }
  }


  }  // anonymous namespace


  void DenseUPolyRingBase::myMul(RawPtr rawlhs, ConstRawPtr rawf, ConstRawPtr rawg) const
  {
    ring P(this);
    RingElemAlias g(P, rawg);

    RingElem ans(P);
    for (long i=0; i<myDegPlus1(rawf); ++i)
      myAddMulLM(raw(ans), raw(myCoeff(rawf,i)), i, rawg);
    mySwap(raw(ans), rawlhs); // really an assignment
  }


  bool DenseUPolyRingBase::myIsDivisible(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    ring P(this);
    RingElem r = RingElemAlias(P, rawx);
    RingElem ans(P);
    DivMod(ans, r, RingElemAlias(P, rawy));
    mySwap(raw(ans), rawlhs); // really an assignment
    return IsZero(r);
  }


  // code for R[x]: compute gcd in FractionField(R)[x]
  // see also myExgcd
  void DenseUPolyRingBase::myGcd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    RingElem a = RingElemAlias(ring(this), rawx);
    if (IsInvertible(a))  a=1;
    if (myIsZero(rawy))
    {
      mySwap(rawlhs, raw(a));
      return;
    }
    RingElem b = RingElemAlias(ring(this), rawy);
    if (IsInvertible(b))  b=1;
    if (myIsZero(rawx))
    {
      mySwap(rawlhs, raw(b));
      return;
    }
    if (!IsField(myCoeffRing()))
    {
      if (!IsTrueGCDDomain(myCoeffRing()))
        CoCoA_THROW_ERROR(ERR::NYI, "DenseUPolyRingBase::myGcd");
      else
      {
        const FractionField K = NewFractionField(myCoeffRing());
        const DenseUPolyRing Kx = NewPolyRing(K);
        const DenseUPolyRing P(this);
        const RingHom phi = PolyRingHom(P, Kx, CoeffEmbeddingHom(Kx)(EmbeddingHom(K)), indets(Kx));
        RingElem z = gcd(phi(a), phi(b));
        Kx->myClearDenom(raw(z), raw(z));
        RingElem g(P);
        for (long d=myStdDeg(raw(z)) ; d >= 0 ; --d)
          g += monomial(P, num(myCoeff(raw(z), d)), d);
        myDivByCoeff(raw(g), raw(content(g)));
        myMulByCoeff(raw(g), raw(gcd(content(a), content(b))));
        P->mySwap(rawlhs, raw(g));
        return;
      }
    }
    // finally Euclid's Algorithm
    mod(a, b);  // a = remainder(a, b)
    while (!IsZero(a))
    {
      mySwap(raw(a), raw(b));
      mod(a, b);  // a = remainder(a, b)
    }
    if (IsInvertible(b))  b=1;
    mySwap(rawlhs, raw(b));
  }


  // See also myGcd
  void DenseUPolyRingBase::myExgcd(RawPtr rawlhs, RawPtr rawxcofac, RawPtr rawycofac, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    ring P(ring(this));
    RingElem ZeroP(P);
    RingElem OneP(one(P));
    RingElem a(RingElemAlias(P, rawx));
    if (myIsZero(rawy))  // = raw(b)
    {
      mySwap(rawlhs, raw(a));
      mySwap(rawxcofac, raw(OneP));
      mySwap(rawycofac, raw(ZeroP));
      return;
    }
    RingElem b(RingElemAlias(P, rawy));
    if (myIsZero(rawx))  // = raw(a)
    {
      mySwap(rawlhs, raw(b));
      mySwap(rawxcofac, raw(ZeroP));
      mySwap(rawycofac, raw(OneP));
      return;
    }
    if (!IsField(myCoeffRing()))
      CoCoA_THROW_ERROR(ERR::NYI, "DenseUPolyRingBase::myGcd(RawPtr, ConstRawPtr, ConstRawPtr) const");
    // finally Euclid's Algorithm
    RingElem res(P), acofac(P), bcofac(P);
    ExgcdAux(acofac, bcofac, a, b);
    mySwap(rawlhs, raw(a));
    mySwap(rawxcofac, raw(acofac));
    mySwap(rawycofac, raw(bcofac));
  }


  void DenseUPolyRingBase::mySymbols(std::vector<symbol>& SymList) const
  {
    myCoeffRing()->mySymbols(SymList);
    SymList.push_back(myIndetSymbol());
  }


  void DenseUPolyRingBase::myOutput(std::ostream& out, ConstRawPtr rawf) const
  {
    if (!out) return;  // short-cut for bad ostreams

    if (myIsZero(rawf)) { out << "0"; return; }

    const ring& R = myCoeffRing();
    bool PrintStar, IsFirstCoeff=true;
    const bool IsQuotientOfZZ = IsQuotientRing(R) && IsZZ(BaseRing(R));
    //bool IsQuotientOfZZ = IsRingFp(R);
    const bool IsNumberRing = IsZZ(R) || IsQuotientOfZZ || IsQQ(R) || IsRingTwinFloat(R);

    for (long d=myStdDeg(rawf) ; d != 0 ; --d)
    {
      PrintStar = true;
      ConstRefRingElem c = myCoeff(rawf, d);
      if (IsZero(c)) continue;
      //----------------------------------------------------------------------
      // coefficient
      //----------------------------------------------------------------------

      if (!IsNumberRing)
      {
        if (!IsFirstCoeff) out << " +";
        out << '(' << c << ')';
        goto EndCoeffPrint;
      }
      // assuming Q is printed with positive denominator
      if (!IsFirstCoeff)
      {
        if (R->myIsPrintedWithMinus(raw(c))) out << " ";
        else out <<" +";
      }
      if (IsOne(c)) // do not print "1 * "
      {
        PrintStar = false;
        goto EndCoeffPrint;
      }
      if ( !IsQuotientOfZZ && IsMinusOne(c) ) // do not print "-1 * "
      {
        PrintStar = false;
        out << "-";
        goto EndCoeffPrint;
      }
      out << c;
      
    EndCoeffPrint:
      
      //----------------------------------------------------------------------
      // PP
      //----------------------------------------------------------------------
      if (PrintStar)  out << "*";
      out << myIndetSymbol();
      if (d > 1) out << "^" << d;

      IsFirstCoeff = false;
    }      
      //----------------------------------------------------------------------
      // coefficient
      //----------------------------------------------------------------------

    ConstRefRingElem c = myCoeff(rawf, 0);
    if (IsZero(c)) return;
    if (!IsNumberRing)
    {
      if (!IsFirstCoeff) out << " +";
      out << '(' << c << ')';
      return;
    }
      // assuming Q is printed with positive denominator
    if (!IsFirstCoeff)
    {
      if (R->myIsPrintedWithMinus(raw(c))) out << " ";
      else out <<" +";
    }
    out << c;
  }


  void DenseUPolyRingBase::myOutputSelf_OM(OpenMathOutput& OMOut) const
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("setname1", myImplDetails());
    OMOut << myCoeffRing();
    OMOut->mySendApplyEnd();
  }


  bool DenseUPolyRingBase::myIsZero(ConstRawPtr rawx) const
  {
    return myDegPlus1(rawx)==0;
  }


  bool DenseUPolyRingBase::myIsPrintAtom(ConstRawPtr rawx) const
  {
    if (!myIsConstant(rawx)) return false;
    // same as in "output"
    const ring& R = myCoeffRing();

    if (IsZZ(R) || IsQQ(R)) return R->myIsPrintAtom(raw(myLC(rawx)));
    return true;
  }


  void DenseUPolyRingBase::myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawx) const
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("cocoa", "DUPCoeffs");
    // ??? if (IsZero(rawx)) .....
    OMOut << myStdDeg(rawx);

    for (long i=0 ; i <= myStdDeg(rawx) ; ++i)
      OMOut << myCoeff(rawx, i);
    OMOut->mySendApplyEnd();
  }
  

  bool DenseUPolyRingBase::myIsOne(ConstRawPtr rawf) const
  {
    if (myDegPlus1(rawf)!=1) return false;
    return IsOne(myCoeff(rawf, 0));
  }


  bool DenseUPolyRingBase::myIsMinusOne(ConstRawPtr rawf) const
  {
    if (myDegPlus1(rawf)!=1) return false;
    return IsMinusOne(myCoeff(rawf, 0));
  }


  bool DenseUPolyRingBase::myIsConstant(ConstRawPtr rawf) const
  {
    if (myDegPlus1(rawf)>1) return false;
    return true;
  }


  bool DenseUPolyRingBase::myIsIndet(long& IndetIndex, ConstRawPtr rawf) const
  {
    if (myDegPlus1(rawf)!=2) return false;
    if (!(IsOne(myCoeff(rawf,1)) && IsZero(myCoeff(rawf,0))))
      return false;
    IndetIndex = 0;
    return true;
  }


  bool DenseUPolyRingBase::myIsIndetPosPower(ConstRawPtr rawf) const
  {
    long d = -1;
    for (long i=0; i<myDegPlus1(rawf); ++i)
      if (!IsZero(myCoeff(rawf, i)))
      {
        if (!IsOne(myCoeff(rawf, i))) return false;
        if (d != -1) return false; else d = i;
      }
    return true;
  }


  bool DenseUPolyRingBase::myIsEvenPoly(ConstRawPtr rawf) const
  {
    const long d1 = myDegPlus1(rawf);
    if (d1 == 0) return true;
    for (long i=1; i<d1; i+=2)
      if (!IsZero(myCoeff(rawf, i))) return false;
    return true;
  }
  
  bool DenseUPolyRingBase::myIsOddPoly(ConstRawPtr rawf) const
  {
    const long d1 = myDegPlus1(rawf);
    if (d1 == 0) return true;
    for (long i=0; i<d1; i+=2)
      if (!IsZero(myCoeff(rawf, i))) return false;
    return true;
  }


  bool DenseUPolyRingBase::myIsMonomial(ConstRawPtr rawf) const
  {
    if (myIsZero(rawf)) return true;
    for (long i=0; i<myStdDeg(rawf); ++i)
      if (!myIsZero(raw(myCoeff(rawf,i)))) return false;
    return true;
  }


  void DenseUPolyRingBase::myIndetPower(RawPtr rawf, long var, long exp) const
  {
    (void)(var); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(var==0);
    RingElem ans(ring(this));
    myResize(raw(ans), exp+1);
    myAssignNonZeroCoeff(raw(ans), raw(one(myCoeffRing())), exp);
    mySwap(raw(ans), rawf); // do it this way to be exception clean
  }


  long DenseUPolyRingBase::myNumTerms(ConstRawPtr rawf) const
  {
    long nt = 0;
    for (long i=0; i<myDegPlus1(rawf); ++i) 
      if (!myCoeffRing()->myIsZero(raw(myCoeff(rawf,i)))) ++nt;
    return nt;
  }
  
  RingElem DenseUPolyRingBase::myMonomial(ConstRawPtr rawc, unsigned long exp) const
  {
    CoCoA_ASSERT(!myCoeffRing()->myIsZero(rawc));
    RingElem ans(ring(this));
    myResize(raw(ans), exp+1);
    myAssignNonZeroCoeff(raw(ans), rawc, exp);
    return ans;
  }


  bool DenseUPolyRingBase::myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    if (myIsZero(rawx))
    {
      if (myIsZero(rawy)) return true;
      else return false;
    }
    if (myDegPlus1(rawx) != myDegPlus1(rawy)) return false;
    for (long i=0; i<myDegPlus1(rawx); ++i) 
      if (myCoeff(rawx,i) != myCoeff(rawy,i))
        return false;
    return true;
  }


  bool DenseUPolyRingBase::myIsInteger(BigInt& N, ConstRawPtr rawx) const
  {
    if (myIsZero(rawx))
    {
      N = 0;
      return true;
    }
    if (myIsConstant(rawx))
      return myCoeffRing()->myIsInteger(N, raw(myCoeff(rawx, 0)));
    return false;
  }


  bool DenseUPolyRingBase::myIsRational(BigRat& Q, ConstRawPtr rawx) const
  {
    if (myIsZero(rawx)) { Q = 0; return true; }
    if (!myIsConstant(rawx)) return false;
    return IsRational(Q, myCoeff(rawx,0));
  }


//   bool DenseUPolyRingBase::myIsHomog(ConstRawPtr rawf) const
//   {
//     if (myIsConstant(rawf)) { return true; }
//     CoCoA_THROW_ERROR(ERR::NYI, "DenseUPolyRingBase::myIsHomog(ConstRawPtr)");
//     return true;
//   }


  bool DenseUPolyRingBase::myIsInvertible(ConstRawPtr rawx) const
  {
    return (!myIsZero(rawx)) && myIsConstant(rawx) && IsInvertible(myLC(rawx));
  }


  void DenseUPolyRingBase::myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const
  {
    // Assert that we have a genuinely non-trivial case.
    CoCoA_ASSERT(n > 1);
    CoCoA_ASSERT(!myIsZero(rawx) && !myIsOne(rawx) && !myIsMinusOne(rawx));
    long NoUse;
    if (myIsIndet(NoUse, rawx))  //??? is this the right place for this check??
      myIndetPower(rawlhs, 0, n);
    else
      myBinaryPower(rawlhs, rawx, n);
    //    mySequentialPower(rawlhs, rawx, n); //??? BUG/SLUG myBinaryPower better if univariate or coeffs are finite field
  }


  const symbol& DenseUPolyRingBase::myIndetSymbol(long idx) const
  {
    if (idx!=0)
      CoCoA_THROW_ERROR(ERR::BadIndex, "DenseUPolyRingBase::myIndetSymbol(idx)");
    return myIndetSymbol();
  }


  //---- Special functions on RingElem owned by DenseUPolyRing


  //-- IdealImpl ----------------------------------------

  ideal DenseUPolyRingBase::myIdealCtor(const std::vector<RingElem>& gens) const
  {
    return ideal(new IdealImpl(DenseUPolyRing(this), gens)); //??? ugly ???
  }


  DenseUPolyRingBase::IdealImpl::IdealImpl(const DenseUPolyRing& P, const std::vector<RingElem>& gens):
    myP(P),
    myGensValue(gens)
  {
    myTidyGensIsValid = false;
    if (!IsField(CoeffRing(P)))
      CoCoA_THROW_ERROR("NYI ideal of polynomials with coeffs not in a field", "ideal(DenseUPolyRing, gens)");//???
    const long ngens = len(gens);
    for (long i=0; i < ngens; ++i)
      if (owner(gens[i]) != myP) CoCoA_THROW_ERROR(ERR::MixedRings, "DenseUPolyRingBase::IdealImpl ctor");
  }


  IdealBase* DenseUPolyRingBase::IdealImpl::myClone() const
  {
    return new IdealImpl(*this);
  }


  bool DenseUPolyRingBase::IdealImpl::IamZero() const
  {
    if (myTidyGensIsValid) return myTidyGensValue.empty();
    const vector<RingElem>& g = myGensValue;
    const long ngens = len(g);
    for (long i=0; i < ngens; ++i)
      if (!IsZero(g[i])) return false;
    return true;
  }


  void DenseUPolyRingBase::IdealImpl::myTestIsMaximal() const
  {
    if (!IsField(CoeffRing(myP)))
      CoCoA_THROW_ERROR(ERR::NYI, "DenseUPolyRingBase::IdealImpl::myTestIsMaximal() non-field CoeffRing");
    myTestIsPrime();
  }


  void DenseUPolyRingBase::IdealImpl::myTestIsPrimary() const
  {
    CoCoA_THROW_ERROR(ERR::NYI, "DenseUPolyRingBase::IdealImpl::myTestIsPrimary()");
  }


  void DenseUPolyRingBase::IdealImpl::myTestIsPrime() const
  {
    if (IamZero())
    { myAssignPrimeFlag(IsIntegralDomain(CoeffRing(myRing()))); return; }
    if (!IsField(CoeffRing(myP)))
      CoCoA_THROW_ERROR(ERR::NYI, "DenseUPolyRingBase::IdealImpl::myTestIsPrime() non-field CoeffRing");
    // The ring is K[x], thus a PID; hence ideal has a single generator.
    CoCoA_ASSERT(len(myTidyGens(NoCpuTimeLimit()))==1);
    const RingElem& f = myTidyGens(NoCpuTimeLimit())[0];
    myAssignPrimeFlag(!IsConstant(f) && IsIrred(f));
  }

  // redefine myAssignMaximalFlag etc, for PID
  void DenseUPolyRingBase::IdealImpl::myTestIsRadical() const
  {
    if (!IsField(CoeffRing(myP)))
      CoCoA_THROW_ERROR(ERR::NYI, "DenseUPolyRingBase::IdealImpl::myTestIsRadical() non-field CoeffRing");
    // The ring is K[x], thus a PID; hence ideal has a single generator.
    CoCoA_ASSERT(len(myTidyGens(NoCpuTimeLimit()))==1);
    myAssignRadicalFlag(IsSqFree(myTidyGens(NoCpuTimeLimit())[0]));
  }


  void DenseUPolyRingBase::IdealImpl::myReduceMod(RingElemRawPtr rawf) const
  {
    if (myP->myIsZero(rawf)) return;
    if (IamZero()) return;
    RingElem g = RingElemAlias(myP, rawf);
    mod(g, (myTidyGens(NoCpuTimeLimit()))[0]);
    myP->mySwap(rawf, raw(g));
  }


  bool DenseUPolyRingBase::IdealImpl::IhaveElem(RingElemConstRawPtr rawf) const
  {
    RingElem g = RingElemAlias(myP, rawf);
    myReduceMod(raw(g));
    return IsZero(g);
  }


  inline const DenseUPolyRingBase::IdealImpl* DenseUPolyRingBase::IdealImpl::GetPtr(const ideal& I)
  {
    return dynamic_cast<const DenseUPolyRingBase::IdealImpl*>(I.myIdealPtr());
  }


  void DenseUPolyRingBase::IdealImpl::myAdd(const ideal& Jin)
  {
    // ANNA: clever insert (skip 0s, deal with 1s)
    myGensValue.insert(myGensValue.end(), gens(Jin).begin(), gens(Jin).end());
    myTidyGensIsValid = false;
    myTidyGensValue.clear();
  }


  void DenseUPolyRingBase::IdealImpl::myMul(const ideal& Jin)
  {
    CoCoA_THROW_ERROR(ERR::NYI, "DenseUPolyRingBase::IdealImpl::mul");
  }


  void DenseUPolyRingBase::IdealImpl::myIntersect(const ideal& J)
  {
    if (IamZero()) return;
    if (IsZero(J))
    {
      myGensValue.clear();
      myTidyGensValue.clear();
      myTidyGensIsValid = true;
      return;
    }
    if (!IsField(CoeffRing(myP)))
      CoCoA_THROW_ERROR(ERR::NYI, "DenseUPolyRingBase::IdealImpl::intersect non-field CoeffRing");
    RingElem f = myTidyGens(NoCpuTimeLimit())[0];
    RingElem g = TidyGens(J)[0];
    //    ideal I(myP, (g*f)/gcd(f,g));
    RingElem h = (g*f)/gcd(f,g);
    myGensValue.clear();
    myGensValue.push_back(h);
    myTidyGensValue.clear();
    myTidyGensIsValid = false;
  }


  void DenseUPolyRingBase::IdealImpl::myColon(const ideal& J)
  {
    if (!IsField(CoeffRing(myP)))
      CoCoA_THROW_ERROR(ERR::NYI, "DenseUPolyRingBase::IdealImpl::colon non-field CoeffRing");
    if (IamZero()) return;
    if (IsZero(J))
    {
      myGensValue.clear();
      myGensValue.push_back(one(myP));
      myTidyGensValue.clear();
      myTidyGensIsValid = false;
      return;
    }
    RingElem f = myTidyGens(NoCpuTimeLimit())[0];
    RingElem g = TidyGens(J)[0];
    RingElem h = f/gcd(f,g);
    myGensValue.clear();
    myGensValue.push_back(h);
    myTidyGensValue.clear();
    myTidyGensIsValid = false;
  }


  bool DenseUPolyRingBase::IdealImpl::myDivMod(RingElemRawPtr /*rawlhs*/, RingElemConstRawPtr /*rawnum*/, RingElemConstRawPtr /*rawden*/) const
  {
    CoCoA_THROW_ERROR(ERR::NYI, "DenseUPolyRingBase::IdealImpl::myDivMod");
    return false; // just to keep compiler quiet!!
  }


  const std::vector<RingElem>& DenseUPolyRingBase::IdealImpl::myTidyGens(const CpuTimeLimit& CheckForTimeout) const
  {
    if (!IsField(CoeffRing(myP)))
      CoCoA_THROW_ERROR(ERR::NYI, "DenseUPolyRingBase::IdealImpl::myTidyGens non-field CoeffRing");
    if (myTidyGensIsValid) return myTidyGensValue;
    CoCoA_ASSERT(myTidyGensValue.empty());
    RingElem g(myP);
    const long ngens = len(myGensValue);
    CheckForTimeout.myReset(IterationVariability::high);
    for (long i=0; i < ngens; ++i)
    {
      CheckForTimeout("DenseUPolyRingBase::IdealImpl::myTidyGens"); // could this ever be useful?
      myP->myGcd(raw(g), raw(g), raw(myGensValue[i]));
      if (IsInvertible(g)) break; // gcd is 1
    }
    if (!IsZero(g))  myTidyGensValue.push_back(g);
    myTidyGensIsValid = true;
    return myTidyGensValue;
  }


  RingElemAlias coeff(ConstRefRingElem f, long d)
  {
    DenseUPolyRing DUP(owner(f), CoCoA_ERROR_CONTEXT);
//    if (!IsDenseUPolyRing(owner(f)))
//      CoCoA_THROW_ERROR(ERR::BadRing, "coeff(f,d)");
    if (d < 0)
      CoCoA_THROW_ERROR(ERR::NegExp, "coeff(f,d);");
    return DUP->myCoeff(raw(f), AsSignedLong(d));
//    return DenseUPolyRingPtr(owner(f))->myCoeff(raw(f), AsSignedLong(d));
  }


  //-- HomImpl ----------------------------------------

  DenseUPolyRingBase::HomImpl::HomImpl(const DenseUPolyRing& domain, const ring& codomain, const RingHom& CoeffHom, ConstRefRingElem IndetImage):
      RingHomBase(domain, codomain),
      myCoeffHom(CoeffHom),
      myIndetImage(IndetImage)
  {
    // No need to check anything: checks already made when CoeffHom was built.
  }

namespace
{
  // ??? this is essentially the same as ApplyGeneral: it should be improved
  // assume image==0
  void ApplyDPRCodomain(RingElem& image, ConstRefRingElem arg, const RingHom CoeffHom, ConstRefRingElem IndetImage)
  {
    const ring DPR = owner(image);
    RingElem g(DPR);
    if (!IsZero(arg))
      for (long i=0; i<=deg(arg); ++i) 
      {
        RingElem SummandImage = CoeffHom(coeff(arg,i));
        CoCoA_ASSERT(owner(SummandImage) == DPR);
        if (IsZero(SummandImage)) continue; // efficiency hack????
        if (i != 0) 
          SummandImage *= power(IndetImage, i);
        g += SummandImage;
      }
    swap(image, g);
  }


  void ApplyGeneral(RingElem& image, ConstRefRingElem arg, const RingHom CoeffHom, ConstRefRingElem IndetImage)
  {
    ring R = owner(image);
    RingElem g(R);
    if (!IsZero(arg))
      for (long i=0; i<=deg(arg); ++i) 
      {
        RingElem SummandImage = CoeffHom(coeff(arg,i));
        CoCoA_ASSERT(owner(SummandImage) == R);
        if (IsZero(SummandImage)) continue; // efficiency hack????
        if (i != 0) 
          SummandImage *= power(IndetImage, i);
        g += SummandImage;
      }
    swap(image, g);
  }
}  // end of anonymous namespace

  void DenseUPolyRingBase::HomImpl::myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const
  {
    RingElem ans(myCodomain);  // Putting result into ans is exception safe and avoids aliasing problems.
    if ( IsDenseUPolyRing(myCodomain) )
      ApplyDPRCodomain(ans, RingElemAlias(myDomain, rawarg), myCoeffHom, myIndetImage);
    else
      ApplyGeneral(ans, RingElemAlias(myDomain, rawarg), myCoeffHom, myIndetImage);
    myCodomain->mySwap(rawimage, raw(ans));
  }


  void DenseUPolyRingBase::HomImpl::myOutputSelfDetails(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    if (NumIndets(myDomain) == 0) return;
    out << " sending "
        << indet(myDomain, 0) << " to " << myIndetImage;    
  }


  RingHom DenseUPolyRingBase::myHomCtor(const ring& codomain, const RingHom& CoeffHom, const std::vector<RingElem>& IndetImages) const
  {
    // Args already sanity checked by PolyRingHom (see PolyRing.C)
// DON'T KNOW IF I REALLY WANT TO MAKE THIS CHECK...
//       // Check to see if we're building an identity homomorphism
//       if (ring(this) == codomain && IsIdentity(CoeffHom))
//       {
//         bool IndetsFixed = true;
//         for (long i=0; i < myNumIndetsValue; ++i)
//           IndetsFixed &= (myIndetVector[i] == IndetImages[i]);
//         if (IndetsFixed) return IdentityHom(ring(this));
//       }
      // General case
    return RingHom(new HomImpl(DenseUPolyRing(this), codomain, CoeffHom, IndetImages[0]));
  }


  RingHom DenseUPolyRingBase::myCompose(const RingHom& phi, const RingHom& theta) const
  {
    vector<RingElem> IndetImages;
    IndetImages.push_back(phi(theta(myIndets()[0])));

    return myHomCtor(codomain(phi), phi(theta(myCoeffEmbeddingHomCtor())), IndetImages);
  }


  //-- CoeffEmbeddingHomImpl ----------------------------------------

  //---------------------------------------------------------------------------
  // Functions for the class DenseUPolyRingBase::CoeffEmbeddingHomImpl


  DenseUPolyRingBase::CoeffEmbeddingHomImpl::CoeffEmbeddingHomImpl(const DenseUPolyRing& P):
    RingHomEmbeddingBase(CoeffRing(P), P)
  {}


  void DenseUPolyRingBase::CoeffEmbeddingHomImpl::myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const
  {
    const DenseUPolyRing P = myCodomain;
    RingElem ans(P);  // don't use image here for aliasing
    // ??? ANNA profile this:  (probably better to have myMonomial)
    if (!myDomain->myIsZero(rawarg))
      ans = monomial(P, RingElemAlias(myDomain, rawarg), 0);
    P->mySwap(rawimage, raw(ans));
  }


} // end of namespace CoCoA

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/DenseUPolyRing.C,v 1.85 2024/07/03 15:40:27 bigatti Exp $
// $Log: DenseUPolyRing.C,v $
// Revision 1.85  2024/07/03 15:40:27  bigatti
// Summary: Error codes: changed BadIndetIndex into BadIndex
//
// Revision 1.84  2024/07/02 09:57:57  bigatti
// Summary: error codes, first batch:
// ReqUnivariate, ReqNonZero, ReqNonZeroGradingDim, ReqNonZeroModulus,  ReqNonZeroRingElem
//
// Revision 1.83  2024/06/09 20:10:20  abbott
// Summary: Renamed myPrepareForNewLoop to myReset (ProgressReporter & CpuTimeLimit); remine 1705
//
// Revision 1.82  2024/03/16 10:15:12  abbott
// Summary: Use ErrorContext when calling arg-check fns (redmine 1793)
//
// Revision 1.81  2024/01/19 15:44:29  bigatti
// Summary: Redmine #1623: modified printing for RingHom
//
// Revision 1.80  2023/07/04 09:20:12  abbott
// Summary: Changed TimeOut to Timeout
//
// Revision 1.79  2022/02/18 14:11:53  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.78  2022/02/08 20:18:54  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.77  2021/03/03 22:09:32  abbott
// Summary: New enum class (redmine 894)
//
// Revision 1.76  2020/06/17 15:49:22  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.75  2020/02/12 10:48:59  bigatti
// -- changed calls to myAssignBLAHFlagPID
//
// Revision 1.74  2020/02/12 09:01:47  bigatti
// -- changed myTestIsMaximal etc to return void (and consequences)
//
// Revision 1.73  2020/02/11 16:56:40  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.72  2020/02/11 16:12:17  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.71  2019/09/16 17:38:48  abbott
// Summary: Impl myIsEvenPoly, myIsOddPoly
//
// Revision 1.70  2018/05/25 09:24:46  abbott
// Summary: Major redesign of CpuTimeLimit (many consequences)
//
// Revision 1.69  2018/05/18 12:13:36  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.68  2018/04/18 14:25:59  abbott
// Summary: Added missing return
//
// Revision 1.67  2018/03/29 09:36:40  bigatti
// -- added member functions myTestIsRadical, myTestIsPrimary and flags
//
// Revision 1.66  2018/03/20 11:45:54  bigatti
// -- changed myPrimeTest --> myTestIsPrime (etc)
//
// Revision 1.65  2017/03/28 12:56:51  bigatti
// -- changed copyright
//
// Revision 1.64  2016/11/08 13:03:31  bigatti
// -- removed useless comment
//
// Revision 1.63  2016/09/08 14:12:32  bigatti
// -- mySetMaximalFlag --> myAssignMaximalFlag
// -- mySetPrimeFlag --> myAssignPrimeFlag
// -- updated the related code
// -- (still "old design", but better aligned)
//
// Revision 1.62  2015/05/04 13:14:10  bigatti
// -- removed and fixed some comments
//
// Revision 1.61  2015/04/24 15:40:59  bigatti
// -- renamed: myAddMul --> myAddMulLM
// -- renamed: myMoveLM --> myMoveLMToFront
// -- new myMoveLMToBack (used in ReductionCog --> bug in test-TmpMorseGraph??)
//
// Revision 1.60  2014/07/31 14:45:17  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.59  2014/07/30 14:03:34  abbott
// Summary: Changed myAmbientRing into myRing
// Author: JAA
//
// Revision 1.58  2014/07/28 15:42:26  abbott
// Summary: Changed myCoeffEmbeddingHom into myCoeffEmbeddingHomCtor
// Author: JAA
//
// Revision 1.57  2014/07/14 15:05:58  abbott
// Summary: Added include of utils.H
// Author: JAA
//
// Revision 1.56  2014/07/11 15:35:25  bigatti
// -- default implementation of myOutputSelf(OpenMathOutput& OMOut)
//
// Revision 1.55  2014/07/09 13:01:17  abbott
// Summary: Removed AsDenseUPolyRing
// Author: JAA
//
// Revision 1.54  2014/07/08 15:21:56  abbott
// Summary: Updated comment
// Author: JAA
//
// Revision 1.53  2014/07/08 13:14:40  abbott
// Summary: Removed AsQuotientRing; added new defn of BaseRing
// Author: JAA
//
// Revision 1.52  2014/07/08 08:34:41  abbott
// Summary: Removed AsFractionField
// Author: JAA
//
// Revision 1.51  2014/06/17 10:23:29  abbott
// Summary: Added (void)(var) to avoid compiler warning about unused param in myIndetPower
// Author: JAA
//
// Revision 1.50  2014/05/06 13:20:41  abbott
// Summary: Changed names (my)MaxExponents into (my)Deg
// Author: JAA
//
// Revision 1.49  2014/04/30 16:04:56  abbott
// Summary: Replaced X.size() by len(X)
// Author: JAA
//
// Revision 1.48  2014/04/15 14:24:21  abbott
// Summary: Improved/cleaned IdealImpl:myPrimeTest
// Author: JAA
//
// Revision 1.47  2013/06/28 17:03:51  abbott
// Modified semantics of IdealBase::myDivMod;
// it now returns a boolean.
// Several consequential changes.
//
// Revision 1.46  2013/05/31 09:15:21  abbott
// Changed arg type of fn "coeff" from MachineInt to long because it is an index.
//
// Revision 1.45  2012/10/24 12:12:36  abbott
// Changed return type of coeff and myLC.
// Replaced several ctor calls so that they build RingElemAlias.
//
// Revision 1.44  2012/10/17 12:03:19  abbott
// Replaced  RefRingElem  by RingElem&
//
// Revision 1.43  2012/05/30 16:04:55  bigatti
// -- applied "3" convention on bool3 functions and member fields
//
// Revision 1.42  2012/05/29 07:45:23  abbott
// Implemented simplification change to bool3:
//  changed names of the constants,
//  changed names of the testing fns.
//
// Revision 1.41  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.40  2012/05/24 14:53:35  bigatti
// -- changed symbol "index" into "subscripts"
//
// Revision 1.39  2012/05/22 10:02:37  abbott
// Removed IsGCDDomain; substituted by IsTrueGCDDomain.
// Added IsFractionFieldOfGCDDomain.
//
// Revision 1.38  2012/04/27 15:08:58  abbott
// Corrected myContentFrF
//
// Revision 1.37  2012/02/10 10:26:40  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.36  2012/02/08 17:12:37  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.35  2012/01/26 16:47:00  bigatti
// -- changed back_inserter into insert
//
// Revision 1.34  2011/11/09 14:03:40  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.33  2011/08/24 10:24:17  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.32  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.31  2011/06/23 16:04:47  abbott
// Added IamExact mem fn for rings.
// Added myRecvTwinFloat mem fn for rings.
// Added first imple of RingHom from RingTwinFloat to other rings.
//
// Revision 1.30  2011/04/27 08:21:48  bigatti
// -- added gcd with coefficients in GCDDomain
//
// Revision 1.29  2011/03/16 15:38:31  bigatti
// -- added myIsIndetPosPower(f)
//
// Revision 1.28  2011/03/10 16:39:34  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.27  2011/03/01 14:10:47  bigatti
// -- added ClearDenom/myClearDenom
//
// Revision 1.26  2011/02/23 16:04:50  bigatti
// -- more fixes to get gcd==1 when coprime
//
// Revision 1.25  2011/01/18 14:38:43  bigatti
// -- moved **_forC5 functions into CoCoA-5/CoCoALibSupplement:
//    myMonomials_forC5, mySupport_forC5, monomials_forC5, support_forC5,
//    LPP_forC5, LT_forC5, LM_forC5
//
// Revision 1.24  2010/11/30 11:30:49  bigatti
// -- moved IndetsCalled into unique implementation in PolyRing
// -- renamed IndetName --> IndetSymbol
// -- added myIndetSymbol
//
// Revision 1.23  2010/11/25 12:31:22  bigatti
// -- added myIndetsCalled
//
// Revision 1.22  2010/11/05 15:59:49  bigatti
// -- added myMonomials_forC5, mySupport_forC5
//
// Revision 1.21  2010/10/01 15:45:17  bigatti
// -- added mySymbolValue
//
// Revision 1.20  2010/06/10 08:00:02  bigatti
// -- fixed naming conventions
//
// Revision 1.19  2010/02/04 10:12:04  bigatti
// -- added "mul" for ideal (implemented only for SparsePolyRing)
//
// Revision 1.18  2009/10/02 13:27:26  bigatti
// -- unique implementation of myDiv in PolyRing.C
//
// Revision 1.17  2009/09/24 16:24:19  abbott
// Added include directive (after removing it from RingFp.H).
//
// Revision 1.16  2009/07/24 12:26:43  abbott
// Added CommonDenom function for polynomials.
//
// Revision 1.15  2009/07/02 16:32:11  abbott
// Consequential changes stemming from new class QQ, and modified interface to the member
// function RingBase::myIsRational.  Also some new conversion functions.
//
// Revision 1.14  2009/05/22 10:24:12  bigatti
// -- updated def of IsQuotientOfZ in myOutput
//
// Revision 1.13  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.12  2008/12/12 15:00:53  bigatti
// -- added { } to disambiguate "if"
//
// Revision 1.11  2008/04/21 12:55:39  abbott
// Fixed some minor bugs, and improved code layout slightly.
//
// Revision 1.10  2008/03/12 16:47:14  bigatti
// -- added: myExgcd, myPrimeTest, myMaximalTest
//
// Revision 1.9  2007/12/21 12:29:08  bigatti
// -- abstract implementation in DenseUPolyRing of myDiv, myIsDivisible, myIsInvertible, myGcd
// -- abstract implementation in DenseUPolyRing of ideal operations
// -- some cleaning
//
// Revision 1.8  2007/12/07 15:27:01  bigatti
// -- default implementation of "IamOne" in ideal.C
//
// Revision 1.7  2007/12/05 11:06:24  bigatti
// -- changed "size_t StdDeg/myStdDeg(f)" into "long"  (and related functions)
// -- changed "log/myLog(f, i)" into "MaxExponent/myMaxExponent(f, i)"
// -- fixed bug in "IsOne(ideal)" in SparsePolyRing.C
//
// Revision 1.6  2007/12/04 14:13:23  bigatti
// -- added forgotten declaration.... sorry
//
// Revision 1.5  2007/12/04 13:44:50  bigatti
// -- fixed problem with myMinCapacity in constructor
// -- commented out myPowerSmallExp in RingDenseUPolyCleanImpl
// -- refined myPowerSmallExp in DenseUPolyRing
//
// Revision 1.4  2007/10/30 17:14:08  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.3  2007/10/15 12:40:31  bigatti
// -- new implementation for myAddMul
// -- new: coeffs between deg+1 and size are now guaranteed to be 0
//
// Revision 1.2  2007/10/10 14:02:37  bigatti
// -- added myMulBy1MinusXExp
// -- fixed a few little bugs
//
// Revision 1.1  2007/10/05 15:28:56  bigatti
// -- added abstract class DenseUPolyRing for representing dense
//    univariate polynomials
// -- added concrete class RingDenseUPolyClean, the cleanest
//    implementation
// -- just for testing, still horribly inefficient and incomplete
//
