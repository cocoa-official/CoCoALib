//   Copyright (c)  2005-2009,2016  John Abbott and Anna M. Bigatti

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


// Source code for classes RingBase(abstract), RingElem, ConstRefRingElem.
// Also most operations on ring elements.

#include "CoCoA/ring.H"
#include "CoCoA/ring-AutomaticConversion.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/CanonicalHom.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/NumTheory-prime.H"
#include "CoCoA/NumTheory-factor.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingElemInput.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyIter.H" // for IsPthPower
#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/bool3.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/factor.H"
#include "CoCoA/interrupt.H"
#include "CoCoA/symbol.H"

#include <atomic>
//#include <vector>
using std::vector;


namespace CoCoA
{

  // default ctor sets ring to RingZZ()
  ring::ring(): 
      mySmartPtr(RingZZ().myRawPtr())
  {}


  long RingBase::NewRingID()
  {
    static std::atomic<long> RingCounter;
    return RingCounter++; // post-incr, so ZZ has ID==0
  }


  BigInt characteristic(const ring& R)
  {
    return R->myCharacteristic();
  }

  long LogCardinality(const ring& R)
  {
    return R->myLogCardinality();
  }

  // This makes a wasteful copy.
  vector<symbol> symbols(const ring& R)
  {
    vector<symbol> ans;
    R->mySymbols(ans);
    return ans;
  }


  bool3 IsPID3(const ring& R)
  {
    if (IsZZ(R) || IsField(R)) return true3;
    if (IsPolyRing(R))
    {
      return bool3((NumIndets(R) == 1) && IsField(CoeffRing(R)));
    }
    if (IsQuotientRing(R)) return uncertain3;
    return false3;
  }

  bool IsPID(const ring& R)
  {
    const bool3 ans = IsPID3(R);
    if (!IsUncertain3(ans)) return IsTrue3(ans);
    // At this point we know R is a QuotientRing.
    CoCoA_THROW_ERROR1(ERR::NYI);
    return false;
  }


  //----------------------------------------------------------------------
  // Default definitions of some virtual functions.

  void RingBase::myOutputSelfShort(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "RingWithID(" << myID << ")";
  }

  void RingBase::myOutputSelfLong(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    myOutputSelf(out);
  }

  long RingBase::myLogCardinality() const
  {
    CoCoA_THROW_ERROR1("Defined only for finite fields");
    return 0; // just to keep compiler quiet
  }

  bool RingBase::IamTrueGCDDomain() const
  {
    return !IamField();
  }


  bool RingBase::IamOrderedDomain() const
  {
    return false;
  }


  RingElem RingBase::myFieldGen() const /* default throws error -- def'd only for finite fields */
  { CoCoA_THROW_ERROR1("Only for finite fields"); }


  RingElemRawPtr RingBase::myNew(const BigRat& Q) const
  {
    AutoRingElem n(ring(this), myNew(num(Q)));
    AutoRingElem d(ring(this), myNew(den(Q)));
//    RingElemRawPtr n = myNew(num(Q));
//    RingElemRawPtr d = myNew(den(Q));
//    const bool OK = myIsDivisible(n, n, d);
    const bool OK = myIsDivisible(raw(n), raw(n), raw(d));
//    myDelete(d);
    if (!OK)
    {
//      myDelete(n);
      CoCoA_THROW_ERROR1(ERR::EmbedBigRatFailed);
    }
    return release(n);
  }


  RingElemRawPtr RingBase::myNew(const symbol& s) const
  {
    vector<symbol> syms = symbols(ring(this));   // makes a copy of vec
    syms.push_back(s);
    if ( AreDistinct(syms) )  // should be IsInSymbols(...)
      CoCoA_THROW_ERROR1("symbol not in ring");
    return myNew(raw(mySymbolValue(s)));
  }


  RingElemRawPtr RingBase::myNew(ConstRefRingElem rhs) const
  {
    if (owner(rhs)==ring(this)) return myNew(raw(rhs));
    RingHom phi = CanonicalHom(owner(rhs), ring(this));
    return myNew(raw(phi(rhs)));
  }


  void RingBase::myAssign(RawPtr rawlhs, const BigRat& Q) const
  {
    RingElemRawPtr n = myNew(num(Q));
    RingElemRawPtr d = myNew(den(Q));
    const bool OK = myIsDivisible(n, n, d);
    myDelete(d);
    if (!OK)
    {
      myDelete(n);
      CoCoA_THROW_ERROR1(ERR::EmbedBigRatFailed);
    }
    mySwap(rawlhs, n); // really an assignment
    myDelete(n);
  }


  bool RingBase::myIsInvertible(ConstRawPtr rawx) const
  {
    RingElem junk(ring(this));
    return myIsDivisible(raw(junk), raw(myOne()), rawx); // ??? discard quotient???
  }


  bool RingBase::myIsZeroDivisor(ConstRawPtr rawx) const // default impl
  {
    if (myIsZero(rawx)) return true;
    if (IsTrue3(IamIntegralDomain3(true))) return false;
    return !IsZero(colon(ideal(myZero()), ideal(RingElemAlias(ring(this),rawx))));
  }


  bool RingBase::myIsIrred(ConstRawPtr /*rawx*/) const
  {
//     if (myIsZero(rawx)) return CoCoA_THROW_ERROR1(ERR::ReqNonZeroRingElem);
//    if (IamField()) return false;  // subsumed by line below
//     if (myIsInvertible(rawx)) return false;
//     if (IamTrueGCDDomain()) return CoCoA_THROW_ERROR2(ERR::NYI, "RingBase::myIsIrred(rawx) for GCDDomain");
//     CoCoA_THROW_ERROR1(ERR::NotTrueGCDDomain);
    CoCoA_THROW_ERROR1(ERR::ShouldNeverGetHere);
    return true; // NEVER REACHED; just to keep compiler quiet.
  }
  

  // Default defn for myGcd: throws NotTrueGCDDomain error.
  void RingBase::myGcd(RawPtr /*rawlhs*/, ConstRawPtr /*rawx*/, ConstRawPtr /*rawy*/) const
  {
    CoCoA_THROW_ERROR1(ERR::NotTrueGCDDomain);
  }


  void RingBase::myLcm(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    CoCoA_ASSERT(IamTrueGCDDomain()); // if assert fails, some important code is missing.
    if (myIsZero(rawx) || myIsZero(rawy)) { myAssignZero(rawlhs); return;}
    // if (IamField()) { myAssign(rawlhs, raw(myOne())); return; }
    RingElem g(ring(this));
    myGcd(raw(g), rawx, rawy);
    myDiv(raw(g), rawx, raw(g));
    myMul(raw(g), rawy, raw(g));
    // do not set it to 1 if invertible otherwise a*b!=lcm(a,b)*gcd(a,b)
    mySwap(rawlhs, raw(g));
  }


  void RingBase::myGcdQuot(RawPtr rawlhs, RawPtr rawxquot, RawPtr rawyquot, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    CoCoA_ASSERT(IamTrueGCDDomain());
    myGcd(rawlhs, rawx, rawy);
    myDiv(rawxquot, rawx, rawlhs);
    myDiv(rawyquot, rawy, rawlhs);
  }


  // void RingBase::myExgcd(RawPtr /*rawlhs*/, RawPtr /*rawxcofac*/, RawPtr /*rawycofac*/, ConstRawPtr /*rawx*/, ConstRawPtr /*rawy*/) const
  // {
  //   if (IamTrueGCDDomain()) CoCoA_THROW_ERROR1(ERR::NYI);
  //   CoCoA_THROW_ERROR1(ERR::NotTrueGCDDomain);
  // }


  void RingBase::myNormalizeFrac(RawPtr rawnum, RawPtr rawden) const
  {
    CoCoA_ASSERT(IamTrueGCDDomain());
    {
      RingElem g(ring(this));
      myGcd(raw(g), rawnum, rawden);
      if (!myIsOne(raw(g)))
      {
        myDiv(rawnum, rawnum, raw(g));
        myDiv(rawden, rawden, raw(g));
      }
    }
    return myNormalizeFracNoGcd(rawnum, rawden);
  }

  void RingBase::myNormalizeFracNoGcd(RawPtr rawnum, RawPtr rawden) const
  {
    CoCoA_ASSERT(IamTrueGCDDomain());
    if (!IamOrderedDomain()) return;
    if (mySign(rawden) > 0) return;
    myNegate(rawnum, rawnum);  // not exception safe :-(
    myNegate(rawden, rawden);
  }


  // Deal with all trivial cases; pass non-trivial cases to myPowerSmallExp.
  // Assume inputs args are mathematically sensible.
  void RingBase::myPower(RawPtr rawlhs, ConstRawPtr rawx, long n) const  // assumes n >= 0
  {
    CoCoA_ASSERT(n >= 0);
    if (n==0) // zeroth power is trivial, note that 0^0 gives 1
    {
      myAssign(rawlhs, 1);
      return;
    }
    if (n==1 || myIsOne(rawx) || myIsZero(rawx)) { myAssign(rawlhs, rawx); return; }
    if (myIsMinusOne(rawx))
    {
      if ((n&1) != 0) myAssign(rawlhs, rawx); // -1 to odd power
      else myAssign(rawlhs, 1); // -1 to even power
      return;
    }
    // Non-trivial case: x is not -1, 0, or 1, and n > 1.
    // Handle squaring specially:
    if (n == 2) { mySquare(rawlhs, rawx); return; }
    myPowerSmallExp(rawlhs, rawx, n);
  }

  // Deal with all trivial cases; pass non-trivial cases to myPowerSmallExp or myPowerBigExp
  void RingBase::myPower(RawPtr rawlhs, ConstRawPtr rawx, const BigInt& N) const  // assumes N >= 0
  {
    CoCoA_ASSERT(N >= 0);
    if (IsZero(N))
    {
      myAssign(rawlhs, 1);  // note that 0^0 gives 1
      return;
    }
    if (N==1 || myIsOne(rawx) || myIsZero(rawx)) { myAssign(rawlhs, rawx); return; }
    if (myIsMinusOne(rawx))
    {
      if (IsOdd(N)) myAssign(rawlhs, rawx); // -1 to odd power
      else myAssign(rawlhs, 1);             // -1 to even power
      return;
    }
    // Non-trivial case: x is not -1, 0, or 1, and N > 1.
    // Handle squaring specially:
    if (N == 2) { mySquare(rawlhs, rawx); return; }
    // Call myPowerSmallExp or myPowerBigExp depending on value of exponent.
    long n;
    if (IsConvertible(n, N))
      myPowerSmallExp(rawlhs, rawx, n);
    else
      myPowerBigExp(rawlhs, rawx, N);
  }


  // Default implementation
  void RingBase::mySquare(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    myMul(rawlhs, rawx, rawx);
  }


  // Default implementation does nothing.
  void RingBase::mySymbols(vector<symbol>& /*SymList*/) const
  {}


  bool RingBase::myIsPrintAtom(ConstRawPtr /*rawx*/) const
  {
    return false;
  }


//   bool RingBase::myIsPrintedWithMinus(ConstRawPtr /*rawx*/) const
//   {
//     //    CoCoA_THROW_ERROR1(ERR::ShouldNeverGetHere);
//     return false;
//   }


  bool RingBase::myIsMinusOne(ConstRawPtr rawx) const
  {
    RingElem tmp(ring(this));
    myNegate(raw(tmp), rawx);
    return myIsOne(raw(tmp));
  }


  // Default impl: not very efficient.
  bool RingBase::myIsDouble(double& d, ConstRawPtr rawx) const
  {
    if (myIsZero(rawx)) { d = 0; return true; }
    BigRat Q;
    if (!myIsRational(Q, rawx)) return false;
    return IsConvertible(d, Q);
  }


  bool RingBase::myIsZeroAddMul(RawPtr rawlhs, ConstRawPtr rawfact1, ConstRawPtr rawfact2) const
  {
// return IsZero(lhs += x*y);
//return myIsZero(myAssign(lhs,myAdd(lhs,myMul(x,y))));
    RingElem tmp(ring(this));
    myMul(raw(tmp), rawfact1, rawfact2);
    myAdd(rawlhs, rawlhs, raw(tmp));
    return myIsZero(rawlhs);
  }


  bool RingBase::myIsZeroAddMul(RawPtr rawlhs, RawPtr rawtmp, ConstRawPtr rawfact1, ConstRawPtr rawfact2) const
  {
    myMul(rawtmp, rawfact1, rawfact2);
    myAdd(rawlhs, rawlhs, rawtmp);
    return myIsZero(rawlhs);
  }


  // This function should never be called.  It could be called if you forget
  // to define myCmp for an ordered ring.
  int RingBase::myCmp(ConstRawPtr /*rawx*/, ConstRawPtr /*rawy*/) const
  {
    CoCoA_THROW_ERROR1(ERR::ShouldNeverGetHere);
    return 0; // NEVER REACHED; just to keep compiler quiet.
  }


  // Default defn; obviously correct, probably slow
  int RingBase::myCmpAbs(ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    CoCoA_ASSERT(IamOrderedDomain());
    return cmp(abs(RingElemAlias(ring(this),rawx)),
               abs(RingElemAlias(ring(this),rawy)));
  }


  int RingBase::mySign(ConstRawPtr rawx) const
  {
    CoCoA_ASSERT(IamOrderedDomain());
    const int tmp = myCmp(rawx, raw(myZero()));
    if (tmp < 0) return -1;
    if (tmp == 0) return 0;
    return 1;
  }


  // This function should never be called.  It could be called if you forget
  // to define myFloor for an ordered ring.
  BigInt RingBase::myFloor(ConstRawPtr /*rawx*/) const
  {
    CoCoA_THROW_ERROR1(ERR::ShouldNeverGetHere);
    return BigInt(0); // NEVER REACHED; just to keep compiler quiet.
  }

  // This function should never be called.  It could be called if you forget
  // to define myCeil for an ordered ring.
  BigInt RingBase::myCeil(ConstRawPtr /*rawx*/) const
  {
    CoCoA_THROW_ERROR1(ERR::ShouldNeverGetHere);
    return BigInt(0); // NEVER REACHED; just to keep compiler quiet.
  }

  // This function should never be called.  It could be called if you forget
  // to define myNearestInt for an ordered ring.
  BigInt RingBase::myNearestInt(ConstRawPtr /*rawx*/) const
  {
    CoCoA_THROW_ERROR1(ERR::ShouldNeverGetHere);
    return BigInt(0); // NEVER REACHED; just to keep compiler quiet.
  }


  // In myPowerBigExp we are guaranteed that the exponent N is too
  // large to fit into a long (i.e. it is genuinely large).
  // We also know that x is not -1,0,1; i.e. a non-trivial case.
  // By default we simply complain that the exponent is too big.
  void RingBase::myPowerBigExp(RawPtr /*rawlhs*/, ConstRawPtr /*rawx*/, const BigInt& /*N*/) const
  {
    CoCoA_THROW_ERROR1(ERR::ExpTooBig);
  }

  // Direct iteration for computing powers -- good for multivariate polys.
  void RingBase::mySequentialPower(RawPtr rawlhs, ConstRawPtr rawx, long n) const  // assumes n >= 0
  {
    CoCoA_ASSERT(n >= 0);
    // Dispose of the trivial cases; note that 0^0 gives 1
    if (n == 0 || myIsOne(rawx)) { myAssign(rawlhs, 1); return; }
    if (myIsZero(rawx)) { myAssignZero(rawlhs); return; }

    RingElem ans(myOne());
    for (long i = 0; i < n; ++i)
    {
      CheckForInterrupt("RingBase::mySequentialPower");
      myMul(raw(ans), rawx, raw(ans));
    }
    mySwap(rawlhs, raw(ans));
  }


  // This function is private because of severe restrictions on its args
  // n must be strictly positive, and no aliasing between lhs and x.
  void RingBase::myBinaryPowerLoop(RawPtr rawlhs, ConstRawPtr rawx, long n) const  // assumes n >= 1
  {
    CoCoA_ASSERT(n >= 1);
    if (n == 1) { myAssign(rawlhs, rawx); return; }
    if (n == 2) { mySquare(rawlhs, rawx); return; }
    myBinaryPowerLoop(rawlhs, rawx, n/2);  // floor division!
    mySquare(rawlhs, rawlhs);
    if (n&1) myMul(rawlhs, rawlhs, rawx);
  }

  void RingBase::myBinaryPower(RawPtr rawlhs, ConstRawPtr rawx, long n) const  // assumes n >= 0
  {
    CoCoA_ASSERT(n >= 0);
    // Dispose of the trivial cases; note that 0^0 gives 1
    if (n == 0 || myIsOne(rawx)) { myAssign(rawlhs, raw(myOne())); return; }
    if (myIsZero(rawx)) { myAssignZero(rawlhs); return; }
    if (n == 1) { myAssign(rawlhs, rawx); return; }

    RingElem ans(ring(this));
    myBinaryPowerLoop(raw(ans), rawx, n);
    mySwap(rawlhs, raw(ans));
  }


  void RingBase::myBinaryPower(RawPtr rawlhs, ConstRawPtr rawx, const BigInt& N) const
  {
    CoCoA_ASSERT("NEGATIVE EXPONENT" && N >= 0);
    // Dispose of the trivial cases; note that 0^0 gives 1
    if (IsZero(N) || myIsOne(rawx)) { myAssign(rawlhs, 1); return; }
    if (myIsZero(rawx)) { myAssignZero(rawlhs); return; }
    if (N == 1) { myAssign(rawlhs, rawx); return; }

    RingElemAlias X(ring(this), rawx);
    RingElem ans(X);
    const long NumBits = 1+FloorLog2(N);
    for (long BitPos = NumBits-1; BitPos > 0; --BitPos)
    {
      ans *= ans;
      if (mpz_tstbit(mpzref(N), BitPos-1))
        ans *= X;
    }
    mySwap(rawlhs, raw(ans));
  }


  void RingBase::myGcdInField(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    CoCoA_THROW_ERROR1(ERR::ShouldNeverGetHere);
    CoCoA_ASSERT(IamField());
    // In a field gcd always gives 1, unless both args are 0 when it gives 0.
    if (myIsZero(rawx) && myIsZero(rawy))
      myAssignZero(rawlhs);
    else
      myAssign(rawlhs, raw(myOne()));
  }


  /////////////////////////////////////////////////////////////////////////////

  RingElem::RingElem():
      RingElemAlias(RingZZ(), RingZZ()->myNew())
  {}

  RingElem::RingElem(const ring& R, const MachineInt& n):
    RingElemAlias(R, R->myNew(n))
  {}


  RingElem::RingElem(const ring& R, const mpz_t N):
      RingElemAlias(R, R->myNew(BigIntFromMPZ(N)))
  {}


  RingElem::RingElem(const ring& R, const BigInt& N):
    RingElemAlias(R, R->myNew(N))
  {}


  RingElem::RingElem(const ring& R, const mpq_t Q):
      RingElemAlias(R, R->myNew(BigRatFromMPQ(Q)))
  {}

  RingElem::RingElem(const ring& R, const BigRat& Q):
    RingElemAlias(R, R->myNew(Q))
  {}


  RingElem::RingElem(const ring& R, const symbol& s):
    RingElemAlias(R, R->myNew(s))
  {}


  RingElem::RingElem(const ring& R, const std::string& s):
    RingElemAlias(R, R->myNew(ReadExpr(R,s)))
  {}


  RingElem::RingElem(const ring& R, ConstRefRingElem rhs):
    RingElemAlias(R, R->myNew(rhs))
  {}
  

  RingElem::~RingElem()
  {
    myOwner()->myDelete(myRawPtr());
  }


  RingElem& RingElem::operator=(const RingElem& rhs)
  {
    return operator=(static_cast<const RingElemAlias&>(rhs));
  }


  RingElem& RingElem::operator=(RingElem&& rhs)
  {
    swap(*this, rhs);
    return *this;
  }


  RingElem& RingElem::operator=(ConstRefRingElem rhs)
  {
    if (myOwner() == owner(rhs))
    {
      // Case: assignment IN SAME RING
      // self-assignment handled gracefully
      myOwner()->myAssign(myRawPtr(), raw(rhs));
      return *this;
    }

    // Case: assignment CHANGES RING
    // For exception safety: make copy first then delete old value.
    const RingElemRawPtr copy = owner(rhs)->myNew(raw(rhs));
    myOwner()->myDelete(myRawPtr());
    myR = owner(rhs);
    myValuePtr = copy;
    return *this;
  }


  RingElem& RingElem::operator=(const MachineInt& n)
  {
    myOwner()->myAssign(myRawPtr(), n);
    return *this;
  }


  RingElem& RingElem::operator=(const BigInt& N)
  {
    myOwner()->myAssign(myRawPtr(), N);
    return *this;
  }


  RingElem& RingElem::operator=(const BigRat& Q)
  {
    myOwner()->myAssign(myRawPtr(), Q);
    return *this;
  }


  bool operator==(ConstRefRingElem x, ConstRefRingElem y)
  {
    if (owner(x) != owner(y))
      CoCoA_THROW_ERROR1(ERR::MixedRings);
    return owner(x)->myIsEqual(raw(x), raw(y));
  }


  RingElem operator-(ConstRefRingElem x)
  {
    const ring& R = owner(x);
    RingElem ans(R);
    R->myNegate(raw(ans), raw(x));
    return ans;
  }


  namespace // anonymous
  { // June 2020 -- dyadic arith on elems belonging to the same ring

    RingElem add_SameRing(ConstRefRingElem x, ConstRefRingElem y)
    {
      const ring& Rx = owner(x);
      CoCoA_ASSERT(Rx == owner(y));
      RingElem ans(Rx);
      Rx->myAdd(raw(ans), raw(x), raw(y));
      return ans;
    }
    
    RingElem sub_SameRing(ConstRefRingElem x, ConstRefRingElem y)
    {
      const ring& Rx = owner(x);
      CoCoA_ASSERT(Rx == owner(y));
      RingElem ans(Rx);
      Rx->mySub(raw(ans), raw(x), raw(y));
      return ans;
    }
    
    RingElem mul_SameRing(ConstRefRingElem x, ConstRefRingElem y)
    {
      const ring& Rx = owner(x);
      CoCoA_ASSERT(Rx == owner(y));
      RingElem ans(Rx);
      Rx->myMul(raw(ans), raw(x), raw(y));
      return ans;
    }
    
    RingElem div_SameRing(ConstRefRingElem x, ConstRefRingElem y)
    {
      const ring& Rx = owner(x);
      CoCoA_ASSERT(Rx == owner(y));
      if (IsZeroDivisor(y))  CoCoA_THROW_ERROR1(ERR::DivByZero);
      RingElem ans(Rx);
      if (!Rx->myIsDivisible(raw(ans), raw(x), raw(y)))
        CoCoA_THROW_ERROR1(ERR::BadQuot);
      return ans;
    }
    
  }
  

  RingElem operator+(ConstRefRingElem x, ConstRefRingElem y)
  {
    CoCoA_STATIC_ERROR_MESG(ErrMixed, ERR::MixedRings,"RingElem+RingElem");
    const ring& Rx = owner(x);
    const ring& Ry = owner(y);
    if (Rx == Ry)  return add_SameRing(x, y);
    const RingHom promote = AutomaticConversionHom(Rx,Ry,ErrMixed); // throws ErrMixed if not permitted 
    if (codomain(promote) == Ry)
      return add_SameRing(promote(x), y);
    return add_SameRing(x, promote(y));
  }


  RingElem operator-(ConstRefRingElem x, ConstRefRingElem y)
  {
    CoCoA_STATIC_ERROR_MESG(ErrMixed, ERR::MixedRings,"RingElem-RingElem");
    const ring& Rx = owner(x);
    const ring& Ry = owner(y);
    if (Rx == Ry)  return sub_SameRing(x, y);
    const RingHom promote = AutomaticConversionHom(Rx,Ry,ErrMixed); // throws ErrMixed if not permitted 
    if (codomain(promote) == Ry)
      return sub_SameRing(promote(x), y);
    return sub_SameRing(x, promote(y));
  }


  RingElem operator*(ConstRefRingElem x, ConstRefRingElem y)
  {
    CoCoA_STATIC_ERROR_MESG(ErrMixed, ERR::MixedRings,"RingElem*RingElem");
    const ring& Rx = owner(x);
    const ring& Ry = owner(y);
    if (Rx == Ry)  return mul_SameRing(x, y);
    const RingHom promote = AutomaticConversionHom(Rx,Ry,ErrMixed); // throws ErrMixed if not permitted 
    if (codomain(promote) == Ry)
      return mul_SameRing(promote(x), y);
    return mul_SameRing(x, promote(y));
  }


  RingElem operator/(ConstRefRingElem x, ConstRefRingElem y)
  {
    CoCoA_STATIC_ERROR_MESG(ErrMixed, ERR::MixedRings,"RingElem/RingElem");
    const ring& Rx = owner(x);
    const ring& Ry = owner(y);
    if (Rx == Ry)  return div_SameRing(x, y);
    const RingHom promote = AutomaticConversionHom(Rx,Ry,ErrMixed); // throws ErrMixed if not permitted 
    if (codomain(promote) == Ry)
      return div_SameRing(promote(x), y);
    return div_SameRing(x, promote(y));
  }


  RingElem gcd(ConstRefRingElem x, ConstRefRingElem y)
  {
    const ring& Rx = owner(x);
    const ring& Ry = owner(y);
    if (Rx != Ry)
      CoCoA_THROW_ERROR1(ERR::MixedRings);

    if (!IsTrueGCDDomain(Rx))  CoCoA_THROW_ERROR1(ERR::NotTrueGCDDomain);
    RingElem ans(Rx);
    Rx->myGcd(raw(ans), raw(x), raw(y));
    return ans;
  }


  RingElem lcm(ConstRefRingElem x, ConstRefRingElem y)
  {
    const ring& Rx = owner(x);
    const ring& Ry = owner(y);
    if (Rx != Ry)
      CoCoA_THROW_ERROR1(ERR::MixedRings);

    if (!IsTrueGCDDomain(Rx))  CoCoA_THROW_ERROR1(ERR::NotTrueGCDDomain);
    RingElem ans(Rx);
    Rx->myLcm(raw(ans), raw(x), raw(y));
    return ans;
  }


  bool IsCoprime(ConstRefRingElem x, ConstRefRingElem y)
  {
    const ring& Rx = owner(x);
    const ring& Ry = owner(y);
    if (Rx != Ry)
      CoCoA_THROW_ERROR1(ERR::MixedRings);

    if (!IsTrueGCDDomain(Rx))  CoCoA_THROW_ERROR1(ERR::NotTrueGCDDomain);
    const bool RingIsGood = (IsZZ(Rx) || (IsPolyRing(Rx) && IsField(CoeffRing(Rx))));
    if (!RingIsGood)  CoCoA_THROW_ERROR1("Only over ZZ or k[x,y,...]");
    RingElem ans(Rx);
    Rx->myGcd(raw(ans), raw(x), raw(y));
    return IsInvertible(ans);
  }


  // Naive impl -- taken from Fabio Rossi's old CoCoA-3 "galois.pkg"
    RingElem GcdList(const std::vector<RingElem> &L)
    {
      const size_t n = L.size();
      if (n == 0)  CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);

      RingElem g = L[0];
      for (size_t i = 1; i < n; ++i)
      {
        CheckForInterrupt("gcdList");
        g = gcd(g, L[i]);
        // @JA: In coclib.cpkg5, IsOne could maybe be included?
        if (IsOne(g))
          break;
      }
      return g;
    }


  // RingElem exgcd(RingElem& xcoeff, RingElem& ycoeff, ConstRefRingElem x, ConstRefRingElem y)
  // {
  //   const char* const FnName = "exgcd(c1,c2,r1,r2)";
  //   const ring& R = owner(x);
  //   if (owner(y) != R)
  //     CoCoA_THROW_ERROR1(ERR::MixedRings);
  //   if (!IsTrueGCDDomain(R))
  //     CoCoA_THROW_ERROR1(ERR::NotTrueGCDDomain);
  //   RingElem gcd(R);
  //   // BUG: clobbers xcoeff and ycoeff always; use locals and swap into xcoeff, ycoeff???
  //   xcoeff = gcd; // must force into same ring
  //   ycoeff = gcd; // must force into same ring
  //   R->myExgcd(raw(gcd), raw(xcoeff), raw(ycoeff), raw(x), raw(y));
  //   return gcd;
  // }



  void GcdQuot(RingElem& gcd, RingElem& quot1, RingElem& quot2, ConstRefRingElem x, ConstRefRingElem y)
  {
    const ring& R = owner(x);
    if (owner(y) != R)
      CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (!IsTrueGCDDomain(R))
      CoCoA_THROW_ERROR1(ERR::NotTrueGCDDomain);
    // Next 3 lines are conditional in case there is aliasing between the args
    if (owner(gcd) != R)  gcd = zero(R); 
    if (owner(quot1) != R)  quot1 = zero(R); 
    if (owner(quot2) != R)  quot2 = zero(R); 
    R->myGcdQuot(raw(gcd), raw(quot1), raw(quot2), raw(x), raw(y));
  }


  RingElem& operator+=(RingElem& x, ConstRefRingElem y)
  {
    const ring& Rx = owner(x);
    const ring& Ry = owner(y);
    if (Rx != Ry)
    {
      CoCoA_STATIC_ERROR_MESG(ErrMixed, ERR::MixedRings,"RingElem += RingElem");
      const RingHom promote = AutomaticConversionHom(Rx,Ry,ErrMixed); // throws ErrMixed if not permitted 
      if (codomain(promote) == Ry)
        CoCoA_THROW_ERROR1(ERR::MixedRings);
      return x += promote(y); //     Rx->myAdd(raw(x), raw(x), raw(promote(y)));
    }

    Rx->myAdd(raw(x), raw(x), raw(y));
    return x;
  }


  RingElem& operator-=(RingElem& x, ConstRefRingElem y)
  {
    const ring& Rx = owner(x);
    const ring& Ry = owner(y);
    if (Rx != Ry)
    {
      CoCoA_STATIC_ERROR_MESG(ErrMixed, ERR::MixedRings,"RingElem -= RingElem");
      const RingHom promote = AutomaticConversionHom(Rx,Ry,ErrMixed); // throws ErrMixed if not permitted 
      if (codomain(promote) == Ry)
        CoCoA_THROW_ERROR1(ERR::MixedRings);
      return x -= promote(y); //     Rx->mySub(raw(x), raw(x), raw(promote(y)));
    }

    Rx->mySub(raw(x), raw(x), raw(y));
    return x;
  }


  RingElem& operator*=(RingElem& x, ConstRefRingElem y)
  {
    const ring& Rx = owner(x);
    const ring& Ry = owner(y);
    if (Rx != Ry)
    {
      CoCoA_STATIC_ERROR_MESG(ErrMixed, ERR::MixedRings,"RingElem *= RingElem");
      if (IsPolyRing(Rx) && CoeffRing(Rx) == Ry)  { PolyRing(Rx)->myMulByCoeff(raw(x),raw(y)); return x; }
      const RingHom promote = AutomaticConversionHom(Rx,Ry,ErrMixed); // throws ErrMixed if not permitted 
      if (codomain(promote) == Ry)
        CoCoA_THROW_ERROR1(ERR::MixedRings);
      return x *= promote(y); //     Rx->myMul(raw(x), raw(x), raw(promote(y)));
    }

    Rx->myMul(raw(x), raw(x), raw(y));
    return x;
  }


  RingElem& operator/=(RingElem& x, ConstRefRingElem y)
  {
    const ring& Rx = owner(x);
    const ring& Ry = owner(y);
    if (Rx != Ry)
    {
      CoCoA_STATIC_ERROR_MESG(ErrMixed, ERR::MixedRings,"RingElem /= RingElem");
      if (IsPolyRing(Rx) && CoeffRing(Rx) == Ry)  { if (IsZeroDivisor(y))  CoCoA_THROW_ERROR1(ERR::DivByZero); PolyRing(Rx)->myDivByCoeff(raw(x),raw(y)); return x; }
      const RingHom promote = AutomaticConversionHom(Rx,Ry,ErrMixed); // throws ErrMixed if not permitted 
      if (codomain(promote) == Ry)
        CoCoA_THROW_ERROR1(ERR::MixedRings);
      return x /= promote(y);
    }

    if (IsZeroDivisor(y))  CoCoA_THROW_ERROR1(ERR::DivByZero);
    if (!Rx->myIsDivisible(raw(x), raw(x), raw(y)))
      CoCoA_THROW_ERROR1(ERR::BadQuot);
    return x;
  }


  std::ostream& operator<<(std::ostream& out, ConstRefRingElem x)
  {
    if (!out) return out;  // short-cut for bad ostreams
    owner(x)->myOutput(out, raw(x));
    return out;
  }


  OpenMathOutput& operator<<(OpenMathOutput& OMOut, ConstRefRingElem x)
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("cocoa", "RingElem");
    OMOut << owner(x);
    owner(x)->myOutput_OM(OMOut, raw(x));
    OMOut->mySendApplyEnd();
    return OMOut;
  }


  bool IsZero(ConstRefRingElem x)
  {
    return owner(x)->myIsZero(raw(x));
  }


  bool IsOne(ConstRefRingElem x)
  {
    return owner(x)->myIsOne(raw(x));
  }


  bool IsMinusOne(ConstRefRingElem x)
  {
    return owner(x)->myIsMinusOne(raw(x));
  }


  bool IsInteger(BigInt& N, ConstRefRingElem x)
  {
    return owner(x)->myIsInteger(N, raw(x));
  }


  bool IsRational(BigRat& Q, ConstRefRingElem x)
  {
    return owner(x)->myIsRational(Q, raw(x));
  }


  bool IsDouble(double& d, ConstRefRingElem x)
  {
    return owner(x)->myIsDouble(d, raw(x));
  }


  bool IsInvertible(ConstRefRingElem x)
  {
    return owner(x)->myIsInvertible(raw(x));
  }


  bool IsZeroDivisor(ConstRefRingElem x)
  { return owner(x)->myIsZeroDivisor(raw(x)); }
  

  bool IsIrred(ConstRefRingElem x)
  {
    if (IsZero(x))  CoCoA_THROW_ERROR1(ERR::ReqNonZeroRingElem);
    if (IsField(owner(x)))  CoCoA_THROW_ERROR1(ERR::NotElemGCDDomain);
    if (!IsTrueGCDDomain(owner(x)))  CoCoA_THROW_ERROR1(ERR::NotElemGCDDomain);
    if (IsInvertible(x))  CoCoA_THROW_ERROR1(ERR::InvertibleRingElem);
    return owner(x)->myIsIrred(raw(x));
  }


  //-------------------------------------------------------
  // IsDivisible with 2 args (next section has IsDivisible with 3 args)

  namespace // anonymous
  {
    constexpr bool AllowFields = true;
    constexpr bool DoNotAllowFields = false;
    
    bool IsDivisible_SameRing(ConstRefRingElem num, ConstRefRingElem den, bool AllowFieldFlag)
    {
      const ring& Rnum = owner(num);
      CoCoA_ASSERT(Rnum == owner(den));
      if (!AllowFieldFlag && IsField(Rnum))  CoCoA_THROW_ERROR1("Not allowed in a field");
      if (IsZero(den)) return false;
// NO, because den may be zero-div!!!    if (IsZero(num)) return true;
      RingElem quot(Rnum);
      return Rnum->myIsDivisible(raw(quot), raw(num), raw(den));
    }

    bool IsDivisible(ConstRefRingElem num, ConstRefRingElem den, bool AllowFieldFlag)
    {
      static const char* const FnName = "IsDivisible(RingElem,RingElem)";
      CoCoA_STATIC_ERROR_MESG(ErrMixed, ERR::MixedRings, FnName);
      const ring& Rnum = owner(num);
      const ring& Rden = owner(den);
      if (Rnum == Rden)  return IsDivisible_SameRing(num, den, AllowFieldFlag);
      // Here Rnum != Rden
      const RingHom promote = AutomaticConversionHom(Rnum,Rden,ErrMixed); // throws ErrMixed if not permitted 
      if (codomain(promote) == Rnum)
        return IsDivisible_SameRing(num, promote(den), AllowFieldFlag);
      return IsDivisible_SameRing(promote(num), den, AllowFieldFlag);
    }

  } // end of namespace anonymous


  bool IsDivisible(ConstRefRingElem num, ConstRefRingElem den)
  {
    return IsDivisible(num,den, DoNotAllowFields);
  }

  bool IsDivisible_AllowFields(ConstRefRingElem num, ConstRefRingElem den)
  {
    return IsDivisible(num,den, AllowFields);
  }


  bool IsDivisible(const MachineInt& n, ConstRefRingElem r)
  {
    RingElem tmp(owner(r), n);
    return IsDivisible_SameRing(tmp, r, DoNotAllowFields);
  }

  bool IsDivisible_AllowFields(const MachineInt& n, ConstRefRingElem r)
  {
    RingElem tmp(owner(r), n);
    return IsDivisible_SameRing(tmp, r, AllowFields);
  }

  bool IsDivisible(ConstRefRingElem r, const MachineInt& n)
  {
    RingElem tmp(owner(r), n);
    return IsDivisible_SameRing(r, tmp, DoNotAllowFields);
  }

  bool IsDivisible_AllowFields(ConstRefRingElem r, const MachineInt& n)
  {
    RingElem tmp(owner(r), n);
    return IsDivisible_SameRing(r, tmp, AllowFields);
  }

  bool IsDivisible(const BigInt& N, ConstRefRingElem r)
  {
    RingElem tmp(owner(r), N);
    return IsDivisible_SameRing(tmp, r, DoNotAllowFields);
  }

  bool IsDivisible_AllowFields(const BigInt& N, ConstRefRingElem r)
  {
    RingElem tmp(owner(r), N);
    return IsDivisible_SameRing(tmp, r, AllowFields);
  }

  bool IsDivisible(ConstRefRingElem r, const BigInt& N)
  {
    RingElem tmp(owner(r), N);
    return IsDivisible_SameRing(r, tmp, DoNotAllowFields);
  }

  bool IsDivisible_AllowFields(ConstRefRingElem r, const BigInt& N)
  {
    RingElem tmp(owner(r), N);
    return IsDivisible_SameRing(r, tmp, AllowFields);
  }


  //-------------------------------------------------------
  // IsDivisible with 3 args (prev section has IsDivisible with 2 args)


  namespace // anonymous
  {

    bool IsDivisible_SameRing(RingElem& lhs, ConstRefRingElem num, ConstRefRingElem den, bool AllowFieldFlag)
    {
      const ring& Rnum = owner(num);
      CoCoA_ASSERT(Rnum == owner(den));
      if (!AllowFieldFlag && IsField(Rnum))  CoCoA_THROW_ERROR1("Not allowed in a field");
      if (IsZero(den)) return false;
      /// NO, because den may be a zero-div!!   if (IsZero(num)) return true;
      RingElem quot = zero(Rnum);
      if (!Rnum->myIsDivisible(raw(quot), raw(num), raw(den))) return false;
      swap(lhs,quot) /*faster than lhs=quot*/;
      return true;
    }

    bool IsDivisible(RingElem& lhs, ConstRefRingElem num, ConstRefRingElem den, bool AllowFieldFlag)
    {
      static const char* const FnName = "IsDivisible(3args)";
      CoCoA_STATIC_ERROR_MESG(ErrMixed, ERR::MixedRings, FnName);
      const ring& Rnum = owner(num);
      const ring& Rden = owner(den);
      if (Rnum == Rden) return IsDivisible_SameRing(lhs, num, den, AllowFieldFlag);
      // Here Rnum != Rden
      const RingHom promote = AutomaticConversionHom(Rnum,Rden,ErrMixed); // throws ErrMixed if not permitted 
      if (codomain(promote) == Rnum)
        return IsDivisible_SameRing(lhs, num, promote(den), AllowFieldFlag);
      return IsDivisible_SameRing(lhs, promote(num), den, AllowFieldFlag);
    }

  } // end of namespace anonymous

  bool IsDivisible(RingElem& lhs, ConstRefRingElem num, ConstRefRingElem den)
  {
    return IsDivisible(lhs, num, den, DoNotAllowFields);
  }

  bool IsDivisible_AllowFields(RingElem& lhs, ConstRefRingElem num, ConstRefRingElem den)
  {
    return IsDivisible(lhs, num, den, AllowFields);
  }

  bool IsDivisible(RingElem& lhs, const MachineInt& n, ConstRefRingElem r)
  {
    RingElem tmp(owner(r), n);
    return IsDivisible_SameRing(lhs, tmp, r, DoNotAllowFields);
  }

  bool IsDivisible_AllowFields(RingElem& lhs, const MachineInt& n, ConstRefRingElem r)
  {
    RingElem tmp(owner(r), n);
    return IsDivisible_SameRing(lhs, tmp, r, AllowFields);
  }

  bool IsDivisible(RingElem& lhs, ConstRefRingElem r, const MachineInt& n)
  {
    RingElem tmp(owner(r), n);
    return IsDivisible_SameRing(lhs, r, tmp, DoNotAllowFields);
  }

  bool IsDivisible_AllowFields(RingElem& lhs, ConstRefRingElem r, const MachineInt& n)
  {
    RingElem tmp(owner(r), n);
    return IsDivisible_SameRing(lhs, r, tmp, AllowFields);
  }

  bool IsDivisible(RingElem& lhs, const BigInt& N, ConstRefRingElem r)
  {
    RingElem tmp(owner(r), N);
    return IsDivisible_SameRing(lhs, tmp, r, DoNotAllowFields);
  }

  bool IsDivisible_AllowFields(RingElem& lhs, const BigInt& N, ConstRefRingElem r)
  {
    RingElem tmp(owner(r), N);
    return IsDivisible_SameRing(lhs, tmp, r, AllowFields);
  }

  bool IsDivisible(RingElem& lhs, ConstRefRingElem r, const BigInt& N)
  {
    RingElem tmp(owner(r), N);
    return IsDivisible_SameRing(lhs, r, tmp, DoNotAllowFields);
  }

  bool IsDivisible_AllowFields(RingElem& lhs, ConstRefRingElem r, const BigInt& N)
  {
    RingElem tmp(owner(r), N);
    return IsDivisible_SameRing(lhs, r, tmp, AllowFields);
  }


  //-------------------------------------------------------

  // Check the args are mathematically sensible; use myPower to do the work.
  RingElem power(ConstRefRingElem x, const MachineInt& n) // deliberately allow negative exponents
  {
    if (IsZero(x) && IsNegative(n))
      CoCoA_THROW_ERROR1(ERR::BadPwrZero);
    if (IsNegative(n) && !IsInvertible(x))
      CoCoA_THROW_ERROR1(ERR::NotUnit);
    if (!IsSignedLong(n)) return power(x, BigInt(n));
// ??? special cases for x=1,0,-1 ???
    if (!IsNegative(n))
    {
      RingElem ans(owner(x));
      owner(x)->myPower(raw(ans), raw(x), AsSignedLong(n));
      return ans;
    }
    // Here we know n < 0 and x is invertible
    RingElem invx(owner(x), 1);
    invx /= x;
    RingElem ans(owner(x));
    owner(x)->myPower(raw(ans), raw(invx), -AsSignedLong(n));
    return ans;
  }


  // Check the args are mathematically sensible; use myPower to do the work.
  RingElem power(ConstRefRingElem x, const BigInt& N) // deliberately allow negative exponents
  {
    if (IsZero(x) && N < 0)
      CoCoA_THROW_ERROR1(ERR::BadPwrZero);
    if (N < 0 && !IsInvertible(x))
      CoCoA_THROW_ERROR1(ERR::NotUnit);

    // General case for N a large integer
    if (N >= 0)
    {
      RingElem ans(owner(x));
      owner(x)->myPower(raw(ans), raw(x), N);
      return ans;
    }
    // Here we know N < 0 and x is invertible
    RingElem invx(owner(x), 1);
    invx /= x;
    RingElem ans(owner(x));
    owner(x)->myPower(raw(ans), raw(invx), -N);
    return ans;
  }


  RingElem FieldGen(const ring& Fq)    ///< mult generator, only for finite fields
  {
    if (!IsFiniteField(Fq))
      CoCoA_THROW_ERROR1("Arg must be finite field");
    return Fq->myFieldGen();
  }


  bool IsPthPower(ConstRefRingElem x)  ///< only in Fp and Fp[x,y,...]
  {
    const BigInt P = characteristic(owner(x));
    long p;
    if (!IsConvertible(p, P) || !IsPrime(p))
      CoCoA_THROW_ERROR1(ERR::BadArg);
    if (IsFiniteField(owner(x)))  return true;
    if (IsFractionField(owner(x)))  return IsPthPower(num(x)) && IsPthPower(den(x));
    if (!IsPolyRing(owner(x)))
      CoCoA_THROW_ERROR1(ERR::BadArg);
    if (!IsSparsePolyRing(owner(x)))  CoCoA_THROW_ERROR1(ERR::NYI);
    for (SparsePolyIter it=BeginIter(x); !IsEnded(it); ++it)
    {
      if (!IsPower(PP(it), p)) return false;
      if (!IsPthPower(coeff(it))) return false;
    }
    return  true;
  }

  RingElem PthRoot(ConstRefRingElem x) ///< only in Fp and Fp[x,y,...]
  {
    const BigInt P = characteristic(owner(x));
    if (IsFiniteField(owner(x)))
      return power(x, power(P, LogCardinality(owner(x))-1));
    long p;
    if (!IsConvertible(p, P) || !IsPrime(p))
      CoCoA_THROW_ERROR1(ERR::ArgTooBig);  // BUG -- weak impl!
    if (IsFractionField(owner(x)))
    { // BUG BUG IMPERFECT IMPL -- could fail if num & den have non-trivial invertible factors without PthRoots
      const RingElem N = PthRoot(num(x));
      const RingElem D = PthRoot(den(x));
      const RingHom phi = EmbeddingHom(owner(x));
      return phi(N)/phi(D);
    }
    if (!IsPolyRing(owner(x)))
      CoCoA_THROW_ERROR1(ERR::BadArg);
    if (!IsSparsePolyRing(owner(x)))  CoCoA_THROW_ERROR1(ERR::NYI);
    const SparsePolyRing Rx = owner(x);
    RingElem ans(Rx);
    for (SparsePolyIter it=BeginIter(x); !IsEnded(it); ++it)
    {
      ans += monomial(Rx, PthRoot(coeff(it)), root(PP(it), p));
    }
    return ans;
  }


  RingElem radical(ConstRefRingElem x) /// squarefree "factor support"
  {
    const ring& R = owner(x);
    if (!IsZZ(R) && !IsPolyRing(R))
      CoCoA_THROW_ERROR1(ERR::BadRing);
    if (IsZZ(R))
      return RingElem(R, radical(ConvertTo<BigInt>(x)));
    return RadicalOfPoly(x);
  }



  RingElem binomial(RingElem x, const MachineInt& n)
  {
    if (IsNegative(n))  CoCoA_THROW_ERROR2(ERR::ReqNonNegative, "2nd arg");
    ring R = owner(x);
    if (IsZero(n)) return one(R);
    //??? special check for 0 < char(R) <= n  should give ERR:DivByZero
    if (characteristic(R)==0)
    { // special handling if x happens to be an integer (& char(R) = 0)
      BigInt X;
      if (IsInteger(X,x))
        return RingElem(R, binomial(X,n));
    }
    const long N = AsSignedLong(n);
    RingElem ans = x;
    for (long i=1; i < N; ++i)
      ans *= (x-i);
    return ans/factorial(N);
  }

  RingElem binomial(RingElem x, const BigInt& N)
  {
    const ErrorInfo ErrMesg(ERR::ArgTooBig, "binomial(RingElem,BigInt)");
    return binomial(x, ConvertTo<long>(N, ErrMesg));
  }

  //---------------------------------------------------------------------------
  // More syntactic sugar: arithmetic between RingElems and MachineInts

  bool operator==(ConstRefRingElem r, const MachineInt& n)
  {
    return owner(r)->myIsEqual(raw(r), raw(RingElem(owner(r), n)));
  }


  RingElem operator+(ConstRefRingElem r, const MachineInt& n)
  {
    RingElem ans(owner(r));
    owner(r)->myAdd(raw(ans), raw(r), raw(RingElem(owner(r), n)));
    return ans;
  }


  RingElem operator-(ConstRefRingElem r, const MachineInt& n)
  {
    RingElem ans(owner(r));
    owner(r)->mySub(raw(ans), raw(r), raw(RingElem(owner(r), n)));
    return ans;
  }


  RingElem operator*(ConstRefRingElem r, const MachineInt& n)
  {
    RingElem ans(owner(r));
    owner(r)->myMul(raw(ans), raw(r), raw(RingElem(owner(r), n)));
    return ans;
  }


  RingElem operator/(ConstRefRingElem r, const MachineInt& n)
  {
    RingElem ans(owner(r), n); // also use ans to convert n to RingElem
    if (IsZeroDivisor(ans))  CoCoA_THROW_ERROR1(ERR::DivByZero);
    if (!owner(r)->myIsDivisible(raw(ans), raw(r), raw(ans)))
      CoCoA_THROW_ERROR1(ERR::BadQuot);
    return ans;
  }


  RingElem gcd(ConstRefRingElem r, const MachineInt& n)
  {
    RingElem ans(owner(r));
    owner(r)->myGcd(raw(ans), raw(r), raw(RingElem(owner(r), n)));
    return ans;
  }

  RingElem lcm(ConstRefRingElem r, const MachineInt& n)
  {
    RingElem ans(owner(r));
    owner(r)->myLcm(raw(ans), raw(r), raw(RingElem(owner(r), n)));
    return ans;
  }

  bool IsCoprime(ConstRefRingElem r, const MachineInt& n)
  {
    RingElem ans(owner(r));
    owner(r)->myGcd(raw(ans), raw(r), raw(RingElem(owner(r), n)));
    return IsInvertible(ans);
  }


  //---------------------------------------------------------------------------
  // Operations between MachineInt and RingElem.

  bool operator==(const MachineInt& n, ConstRefRingElem r)
  {
    return owner(r)->myIsEqual(raw(r), raw(RingElem(owner(r), n)));
  }


  RingElem operator+(const MachineInt& n, ConstRefRingElem r)
  {
    RingElem ans(owner(r));
    owner(r)->myAdd(raw(ans), raw(RingElem(owner(r), n)), raw(r));
    return ans;
  }


  RingElem operator-(const MachineInt& n, ConstRefRingElem r)
  {
    RingElem ans(owner(r));
    owner(r)->mySub(raw(ans), raw(RingElem(owner(r), n)), raw(r));
    return ans;
  }


  RingElem operator*(const MachineInt& n, ConstRefRingElem r)
  {
    RingElem ans(owner(r));
    owner(r)->myMul(raw(ans), raw(RingElem(owner(r), n)), raw(r));
    return ans;
  }


  RingElem operator/(const MachineInt& n, ConstRefRingElem r)
  {
    if (IsZeroDivisor(r))  CoCoA_THROW_ERROR1(ERR::DivByZero);
    RingElem ans(owner(r), n); // use ans to convert n to a RingElem
    if (!owner(r)->myIsDivisible(raw(ans), raw(ans), raw(r)))
      CoCoA_THROW_ERROR1(ERR::BadQuot);
    return ans;
  }


  RingElem gcd(const MachineInt& n, ConstRefRingElem r)
  {
    RingElem ans(owner(r));
    owner(r)->myGcd(raw(ans), raw(RingElem(owner(r), n)), raw(r));
    return ans;
  }

  RingElem lcm(const MachineInt& n, ConstRefRingElem r)
  {
    RingElem ans(owner(r));
    owner(r)->myLcm(raw(ans), raw(RingElem(owner(r), n)), raw(r));
    return ans;
  }

  bool IsCoprime(const MachineInt& n, ConstRefRingElem r)
  {
    RingElem ans(owner(r));
    owner(r)->myGcd(raw(ans), raw(RingElem(owner(r), n)), raw(r));
    return IsInvertible(ans);
  }


  RingElem& operator+=(RingElem& r, const MachineInt& n)
  {
    owner(r)->myAdd(raw(r), raw(r), raw(RingElem(owner(r), n)));
    return r;
  }


  RingElem& operator-=(RingElem& r, const MachineInt& n)
  {
    owner(r)->mySub(raw(r), raw(r), raw(RingElem(owner(r), n)));
    return r;
  }


  RingElem& operator*=(RingElem& r, const MachineInt& n)
  {
    owner(r)->myMul(raw(r), raw(r), raw(RingElem(owner(r), n)));
    return r;
  }


  RingElem& operator/=(RingElem& r, const MachineInt& n)
  {
    RingElem den(owner(r), n);
    if (IsZeroDivisor(den))  CoCoA_THROW_ERROR1(ERR::DivByZero);
    if (!owner(r)->myIsDivisible(raw(r), raw(r), raw(den)))
      CoCoA_THROW_ERROR1(ERR::BadQuot);
    return r;
  }


  /////////////////////////////////////////////////////////////////////////////
  // More syntactic sugar: arithmetic between RingElems and BigInts

  bool operator==(ConstRefRingElem r, const BigInt& N)
  {
    return owner(r)->myIsEqual(raw(r), raw(RingElem(owner(r), N)));
  }


  RingElem operator+(ConstRefRingElem r, const BigInt& N)
  {
    RingElem ans(owner(r));
    owner(r)->myAdd(raw(ans), raw(r), raw(RingElem(owner(r), N)));
    return ans;
  }


  RingElem operator-(ConstRefRingElem r, const BigInt& N)
  {
    RingElem ans(owner(r));
    owner(r)->mySub(raw(ans), raw(r), raw(RingElem(owner(r), N)));
    return ans;
  }


  RingElem operator*(ConstRefRingElem r, const BigInt& N)
  {
    RingElem ans(owner(r));
    owner(r)->myMul(raw(ans), raw(r), raw(RingElem(owner(r), N)));
    return ans;
  }


  RingElem operator/(ConstRefRingElem r, const BigInt& N)
  {
    RingElem ans(owner(r), N); // also use ans to convert N to RingElem
    if (IsZeroDivisor(ans))  CoCoA_THROW_ERROR1(ERR::DivByZero);
    if (!owner(r)->myIsDivisible(raw(ans), raw(r), raw(ans)))
      CoCoA_THROW_ERROR1(ERR::BadQuot);
    return ans;
  }


  RingElem gcd(ConstRefRingElem r, const BigInt& N)
  {
    RingElem ans(owner(r));
    owner(r)->myGcd(raw(ans), raw(r), raw(RingElem(owner(r), N)));
    return ans;
  }

  RingElem lcm(ConstRefRingElem r, const BigInt& N)
  {
    RingElem ans(owner(r));
    owner(r)->myLcm(raw(ans), raw(r), raw(RingElem(owner(r), N)));
    return ans;
  }

  bool IsCoprime(ConstRefRingElem r, const BigInt& N)
  {
    RingElem ans(owner(r));
    owner(r)->myGcd(raw(ans), raw(r), raw(RingElem(owner(r), N)));
    return IsInvertible(ans);
  }


  bool operator==(const BigInt& N, ConstRefRingElem r)
  {
    return owner(r)->myIsEqual(raw(r), raw(RingElem(owner(r), N)));
  }


  RingElem operator+(const BigInt& N, ConstRefRingElem r)
  {
    RingElem ans(owner(r), N);
    owner(r)->myAdd(raw(ans), raw(ans), raw(r));
    return ans;
  }


  RingElem operator-(const BigInt& N, ConstRefRingElem r)
  {
    RingElem ans(owner(r), N);
    owner(r)->mySub(raw(ans), raw(ans), raw(r));
    return ans;
  }


  RingElem operator*(const BigInt& N, ConstRefRingElem r)
  {
    RingElem ans(owner(r), N);
    owner(r)->myMul(raw(ans), raw(ans), raw(r));
    return ans;
  }


  RingElem operator/(const BigInt& N, ConstRefRingElem r)
  {
    if (IsZeroDivisor(r))  CoCoA_THROW_ERROR1(ERR::DivByZero);
    RingElem ans(owner(r), N);
    if (!owner(r)->myIsDivisible(raw(ans), raw(ans), raw(r)))
      CoCoA_THROW_ERROR1(ERR::BadQuot);
    return ans;
  }


  RingElem gcd(const BigInt& N, ConstRefRingElem r)
  {
    RingElem ans(owner(r));
    owner(r)->myGcd(raw(ans), raw(RingElem(owner(r), N)), raw(r));
    return ans;
  }

  RingElem lcm(const BigInt& N, ConstRefRingElem r)
  {
    RingElem ans(owner(r));
    owner(r)->myLcm(raw(ans), raw(RingElem(owner(r), N)), raw(r));
    return ans;
  }

  bool IsCoprime(const BigInt& N, ConstRefRingElem r)
  {
    RingElem ans(owner(r));
    owner(r)->myGcd(raw(ans), raw(RingElem(owner(r), N)), raw(r));
    return IsInvertible(ans);
  }


  RingElem& operator+=(RingElem& r, const BigInt& N)
  {
    owner(r)->myAdd(raw(r), raw(r), raw(RingElem(owner(r), N)));
    return r;
  }


  RingElem& operator-=(RingElem& r, const BigInt& N)
  {
    owner(r)->mySub(raw(r), raw(r), raw(RingElem(owner(r), N)));
    return r;
  }


  RingElem& operator*=(RingElem& r, const BigInt& N)
  {
    owner(r)->myMul(raw(r), raw(r), raw(RingElem(owner(r), N)));
    return r;
  }


  RingElem& operator/=(RingElem& r, const BigInt& N)
  {
    RingElem den(owner(r), N);
    if (IsZeroDivisor(den))  CoCoA_THROW_ERROR1(ERR::DivByZero);
    if (!owner(r)->myIsDivisible(raw(r), raw(r), raw(den)))
      CoCoA_THROW_ERROR1(ERR::BadQuot);
    return r;
  }


  std::ostream& operator<<(std::ostream& out, const ring& R)
  {
    if (!out) return out;  // short-cut for bad ostreams
    R->myOutputSelf(out);
    return out;
  }


  OpenMathOutput& operator<<(OpenMathOutput& OMOut, const ring& R)
  {
    R->myOutputSelf_OM(OMOut);
    return OMOut;
  }


  // Comparisons for arithmetically ordered rings

  int sign(ConstRefRingElem x)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(x)->mySign(raw(x));
  }


  RingElem abs(ConstRefRingElem x)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    if (x < 0)  return -x;
    return x;
  }


  BigInt floor(ConstRefRingElem x)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(x)->myFloor(raw(x));
  }


  BigInt ceil(ConstRefRingElem x)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(x)->myCeil(raw(x));
  }


  BigInt NearestInt(ConstRefRingElem x)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(x)->myNearestInt(raw(x));
  }


  int CmpDouble(ConstRefRingElem x, double z)
  {
    // This impl is simple rather than efficient.
    const BigRat Q = ConvertTo<BigRat>(z);
    return cmp(den(Q)*x, num(Q));
  }


  int cmp(ConstRefRingElem x, ConstRefRingElem y)
  {
    if (owner(x) != owner(y))
      CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(x)->myCmp(raw(x), raw(y));
  }


  int CmpAbs(ConstRefRingElem x, ConstRefRingElem y)
  {
    if (owner(x) != owner(y))
      CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(x)->myCmpAbs(raw(x), raw(y));
  }


  bool operator<(ConstRefRingElem x, ConstRefRingElem y)
  {
    if (owner(x) != owner(y))
      CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(x)->myCmp(raw(x), raw(y)) < 0;
  }


  bool operator<=(ConstRefRingElem x, ConstRefRingElem y)
  {
    if (owner(x) != owner(y))
      CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(x)->myCmp(raw(x), raw(y)) <= 0;
  }


  bool operator>(ConstRefRingElem x, ConstRefRingElem y)
  {
    if (owner(x) != owner(y))
      CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(x)->myCmp(raw(x), raw(y)) > 0;
  }


  bool operator>=(ConstRefRingElem x, ConstRefRingElem y)
  {
    if (owner(x) != owner(y))
      CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(x)->myCmp(raw(x), raw(y)) >= 0;
  }


  /////////////////////////////////////////////////////////////////////////////
  // First arg is a MachineInt, second is a RingElem

  int cmp(const MachineInt& x, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    if (IsZero(x)) return -sign(y);
    return owner(y)->myCmp(raw(RingElem(owner(y),x)), raw(y));
  }


  bool operator<(const MachineInt& x, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    if (IsZero(x)) return sign(y) > 0;
    return owner(y)->myCmp(raw(RingElem(owner(y),x)), raw(y)) < 0;
  }


  bool operator<=(const MachineInt& x, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    if (IsZero(x)) return sign(y) >= 0;
    return owner(y)->myCmp(raw(RingElem(owner(y),x)), raw(y)) <= 0;
  }


  bool operator>(const MachineInt& x, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    if (IsZero(x)) return sign(y) < 0;
    return owner(y)->myCmp(raw(RingElem(owner(y),x)), raw(y)) > 0;
  }


  bool operator>=(const MachineInt& x, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    if (IsZero(x)) return sign(y) <= 0;
    return owner(y)->myCmp(raw(RingElem(owner(y),x)), raw(y)) >= 0;
  }


  /////////////////////////////////////////////////////////////////////////////
  // First arg is a RingElem, second is a MachineInt


  int cmp(ConstRefRingElem x, const MachineInt& y)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    if (IsZero(y)) return sign(x);
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),y)));
  }


  bool operator<(ConstRefRingElem x, const MachineInt& y)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    if (IsZero(y)) return sign(x) < 0;
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),y))) < 0;
  }


  bool operator<=(ConstRefRingElem x, const MachineInt& y)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    if (IsZero(y)) return sign(x) <= 0;
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),y))) <= 0;
  }


  bool operator>(ConstRefRingElem x, const MachineInt& y)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    if (IsZero(y)) return sign(x) > 0;
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),y))) > 0;
  }


  bool operator>=(ConstRefRingElem x, const MachineInt& y)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    if (IsZero(y)) return sign(x) >= 0;
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),y))) >= 0;
  }

  ///////////////////////////////////////////////////////////////////////////
  // Comparisons between BigInt and RingElem
  int cmp(const BigInt& x, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(y)->myCmp(raw(RingElem(owner(y),x)), raw(y));
  }


  bool operator<(const BigInt& x, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(y)->myCmp(raw(RingElem(owner(y),x)), raw(y)) < 0;
  }


  bool operator<=(const BigInt& x, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(y)->myCmp(raw(RingElem(owner(y),x)), raw(y)) <= 0;
  }


  bool operator>(const BigInt& x, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(y)->myCmp(raw(RingElem(owner(y),x)), raw(y)) > 0;
  }


  bool operator>=(const BigInt& x, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(y)->myCmp(raw(RingElem(owner(y),x)), raw(y)) >= 0;
  }


  /////////////////////////////////////////////////////////////////////////////
  // First arg is a RingElem, second is a const BigInt&


  int cmp(ConstRefRingElem x, const BigInt& y)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),y)));
  }


  bool operator<(ConstRefRingElem x, const BigInt& y)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),y))) < 0;
  }


  bool operator<=(ConstRefRingElem x, const BigInt& y)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),y))) <= 0;
  }


  bool operator>(ConstRefRingElem x, const BigInt& y)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),y))) > 0;
  }


  bool operator>=(ConstRefRingElem x, const BigInt& y)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),y))) >= 0;
  }


  /////////////////////////////////////////////////////////////////////////////
  // More syntactic sugar: arithmetic between RingElems and BigRats

  bool operator==(ConstRefRingElem r, const BigRat& Q)
  {
    return owner(r)->myIsEqual(raw(r), raw(RingElem(owner(r), Q)));
  }


  RingElem operator+(ConstRefRingElem r, const BigRat& Q)
  {
    RingElem ans(owner(r));
    owner(r)->myAdd(raw(ans), raw(r), raw(RingElem(owner(r), Q)));
    return ans;
  }


  RingElem operator-(ConstRefRingElem r, const BigRat& Q)
  {
    RingElem ans(owner(r));
    owner(r)->mySub(raw(ans), raw(r), raw(RingElem(owner(r), Q)));
    return ans;
  }


  RingElem operator*(ConstRefRingElem r, const BigRat& Q)
  {
    RingElem ans(owner(r));
    owner(r)->myMul(raw(ans), raw(r), raw(RingElem(owner(r), Q)));
    return ans;
  }


  RingElem operator/(ConstRefRingElem r, const BigRat& Q)
  {
    RingElem ans(owner(r), Q); // also use ans to convert N to RingElem
    if (IsZeroDivisor(ans))  CoCoA_THROW_ERROR1(ERR::DivByZero);
    if (!owner(r)->myIsDivisible(raw(ans), raw(r), raw(ans)))
      CoCoA_THROW_ERROR1(ERR::BadQuot);
    return ans;
  }


  bool operator==(const BigRat& Q, ConstRefRingElem r)
  {
    return owner(r)->myIsEqual(raw(r), raw(RingElem(owner(r), Q)));
  }


  RingElem operator+(const BigRat& Q, ConstRefRingElem r)
  {
    RingElem ans(owner(r), Q);
    owner(r)->myAdd(raw(ans), raw(ans), raw(r));
    return ans;
  }


  RingElem operator-(const BigRat& Q, ConstRefRingElem r)
  {
    RingElem ans(owner(r), Q);
    owner(r)->mySub(raw(ans), raw(ans), raw(r));
    return ans;
  }


  RingElem operator*(const BigRat& Q, ConstRefRingElem r)
  {
    RingElem ans(owner(r), Q);
    owner(r)->myMul(raw(ans), raw(ans), raw(r));
    return ans;
  }


  RingElem operator/(const BigRat& Q, ConstRefRingElem r)
  {
    if (IsZeroDivisor(r)) CoCoA_THROW_ERROR1(ERR::DivByZero);
    RingElem ans(owner(r), Q);
    if (!owner(r)->myIsDivisible(raw(ans), raw(ans), raw(r)))
      CoCoA_THROW_ERROR1(ERR::BadQuot);
    return ans;
  }


  RingElem& operator+=(RingElem& r, const BigRat& Q)
  {
    owner(r)->myAdd(raw(r), raw(r), raw(RingElem(owner(r), Q)));
    return r;
  }


  RingElem& operator-=(RingElem& r, const BigRat& Q)
  {
    owner(r)->mySub(raw(r), raw(r), raw(RingElem(owner(r), Q)));
    return r;
  }


  RingElem& operator*=(RingElem& r, const BigRat& Q)
  {
    owner(r)->myMul(raw(r), raw(r), raw(RingElem(owner(r), Q)));
    return r;
  }


  RingElem& operator/=(RingElem& r, const BigRat& Q)
  {
    RingElem den(owner(r), Q);
    if (IsZeroDivisor(den))  CoCoA_THROW_ERROR1(ERR::DivByZero);
    if (!owner(r)->myIsDivisible(raw(r), raw(r), raw(den)))
      CoCoA_THROW_ERROR1(ERR::BadQuot);
    return r;
  }


  ///////////////////////////////////////////////////////////////////////////
  // Comparisons between BigRat and RingElem
  int cmp(const BigRat& Q, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(y)->myCmp(raw(RingElem(owner(y),Q)), raw(y));
  }


  bool operator<(const BigRat& Q, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(y)->myCmp(raw(RingElem(owner(y),Q)), raw(y)) < 0;
  }


  bool operator<=(const BigRat& Q, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(y)->myCmp(raw(RingElem(owner(y),Q)), raw(y)) <= 0;
  }


  bool operator>(const BigRat& Q, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(y)->myCmp(raw(RingElem(owner(y),Q)), raw(y)) > 0;
  }


  bool operator>=(const BigRat& Q, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(y)->myCmp(raw(RingElem(owner(y),Q)), raw(y)) >= 0;
  }


  /////////////////////////////////////////////////////////////////////////////
  // First arg is a RingElem, second is a BigRat


  int cmp(ConstRefRingElem x, const BigRat& Q)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),Q)));
  }


  bool operator<(ConstRefRingElem x, const BigRat& Q)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),Q))) < 0;
  }


  bool operator<=(ConstRefRingElem x, const BigRat& Q)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),Q))) <= 0;
  }


  bool operator>(ConstRefRingElem x, const BigRat& Q)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),Q))) > 0;
  }


  bool operator>=(ConstRefRingElem x, const BigRat& Q)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),Q))) >= 0;
  }




} // end of namespace CoCoA
