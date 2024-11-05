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
    CoCoA_THROW_ERROR(ERR::NYI, "IsPID for quotient ring");
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
    CoCoA_THROW_ERROR("Defined only for finite fields", "LogCardinality");
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
  { CoCoA_THROW_ERROR("Only for finite fields", "RingBase::myFieldGen"); }


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
      CoCoA_THROW_ERROR(ERR::EmbedBigRatFailed, "RingBase::myNew(BigRat) -- generic embedding function");
    }
    return release(n);
  }


  RingElemRawPtr RingBase::myNew(const symbol& s) const
  {
    vector<symbol> syms = symbols(ring(this));   // makes a copy of vec
    syms.push_back(s);
    if ( AreDistinct(syms) )  // should be IsInSymbols(...)
      CoCoA_THROW_ERROR("symbol not in ring", "myNew(symbol)");
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
      CoCoA_THROW_ERROR(ERR::EmbedBigRatFailed, "RingBase::myAssign(RawPtr,BigRat) -- generic embedding function");
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
//     if (myIsZero(rawx)) return CoCoA_THROW_ERROR(ERR::ReqNonZeroRingElem, "RingBase::myIsIrred");
//    if (IamField()) return false;  // subsumed by line below
//     if (myIsInvertible(rawx)) return false;
//     if (IamTrueGCDDomain()) return CoCoA_THROW_ERROR(ERR::NYI, "RingBase::myIsIrred(rawx) for GCDDomain");
//     CoCoA_THROW_ERROR(ERR::NotTrueGCDDomain, "RingBase::myIsIrred(rawx)");
    CoCoA_THROW_ERROR(ERR::ShouldNeverGetHere, "RingBase::myIsIrred(rawx)");
    return true; // NEVER REACHED; just to keep compiler quiet.
  }
  

  // Default defn for myGcd: throws NotTrueGCDDomain error.
  void RingBase::myGcd(RawPtr /*rawlhs*/, ConstRawPtr /*rawx*/, ConstRawPtr /*rawy*/) const
  {
    CoCoA_THROW_ERROR(ERR::NotTrueGCDDomain, "RingBase::myGcd");
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
  //   if (IamTrueGCDDomain()) CoCoA_THROW_ERROR(ERR::NYI, "RingBase::myExgcd(,,)");
  //   CoCoA_THROW_ERROR(ERR::NotTrueGCDDomain, "RingBase::myExgcd(,,)");
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
//     //    CoCoA_THROW_ERROR(ERR::ShouldNeverGetHere, "RingBase::myIsPrintedWithMinus");
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
    CoCoA_THROW_ERROR(ERR::ShouldNeverGetHere, "RingBase::myCmp should never be called");
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
    CoCoA_THROW_ERROR(ERR::ShouldNeverGetHere, "RingBase::myFloor should never be called");
    return BigInt(0); // NEVER REACHED; just to keep compiler quiet.
  }

  // This function should never be called.  It could be called if you forget
  // to define myCeil for an ordered ring.
  BigInt RingBase::myCeil(ConstRawPtr /*rawx*/) const
  {
    CoCoA_THROW_ERROR(ERR::ShouldNeverGetHere, "RingBase::myCeil should never be called");
    return BigInt(0); // NEVER REACHED; just to keep compiler quiet.
  }

  // This function should never be called.  It could be called if you forget
  // to define myNearestInt for an ordered ring.
  BigInt RingBase::myNearestInt(ConstRawPtr /*rawx*/) const
  {
    CoCoA_THROW_ERROR(ERR::ShouldNeverGetHere, "RingBase::myNearestInt should never be called");
    return BigInt(0); // NEVER REACHED; just to keep compiler quiet.
  }


  // In myPowerBigExp we are guaranteed that the exponent N is too
  // large to fit into a long (i.e. it is genuinely large).
  // We also know that x is not -1,0,1; i.e. a non-trivial case.
  // By default we simply complain that the exponent is too big.
  void RingBase::myPowerBigExp(RawPtr /*rawlhs*/, ConstRawPtr /*rawx*/, const BigInt& /*N*/) const
  {
    CoCoA_THROW_ERROR(ERR::ExpTooBig, "power(r,N)");
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
    CoCoA_THROW_ERROR(ERR::ShouldNeverGetHere, "GCD in field");
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
      CoCoA_THROW_ERROR(ERR::MixedRings, "RingElem == RingElem");
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
      if (IsZeroDivisor(y)) CoCoA_THROW_ERROR(ERR::DivByZero, "RingElem / RingElem");
      RingElem ans(Rx);
      if (!Rx->myIsDivisible(raw(ans), raw(x), raw(y)))
        CoCoA_THROW_ERROR(ERR::BadQuot, "RingElem / RingElem");
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
      CoCoA_THROW_ERROR(ERR::MixedRings, "gcd(RingElem, RingElem)");

    if (!IsTrueGCDDomain(Rx)) CoCoA_THROW_ERROR(ERR::NotTrueGCDDomain, "gcd(x,y)");
    RingElem ans(Rx);
    Rx->myGcd(raw(ans), raw(x), raw(y));
    return ans;
  }


  RingElem lcm(ConstRefRingElem x, ConstRefRingElem y)
  {
    const ring& Rx = owner(x);
    const ring& Ry = owner(y);
    if (Rx != Ry)
      CoCoA_THROW_ERROR(ERR::MixedRings, "lcm(RingElem, RingElem)");

    if (!IsTrueGCDDomain(Rx)) CoCoA_THROW_ERROR(ERR::NotTrueGCDDomain, "lcm(x,y)");
    RingElem ans(Rx);
    Rx->myLcm(raw(ans), raw(x), raw(y));
    return ans;
  }


  bool IsCoprime(ConstRefRingElem x, ConstRefRingElem y)
  {
    const ring& Rx = owner(x);
    const ring& Ry = owner(y);
    if (Rx != Ry)
      CoCoA_THROW_ERROR(ERR::MixedRings, "IsCoprime(RingElem, RingElem)");

    if (!IsTrueGCDDomain(Rx)) CoCoA_THROW_ERROR(ERR::NotTrueGCDDomain, "IsCoprime(x,y)");
    const bool RingIsGood = (IsZZ(Rx) || (IsPolyRing(Rx) && IsField(CoeffRing(Rx))));
    if (!RingIsGood) CoCoA_THROW_ERROR("Only over ZZ or k[x,y,...]", "IsCoprime(x,y)");
    RingElem ans(Rx);
    Rx->myGcd(raw(ans), raw(x), raw(y));
    return IsInvertible(ans);
  }


  // Naive impl -- taken from Fabio Rossi's old CoCoA-3 "galois.pkg"
    RingElem GcdList(const std::vector<RingElem> &L)
    {
      const size_t n = L.size();
      if (n == 0) CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);

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
  //     CoCoA_THROW_ERROR(ERR::MixedRings, FnName);
  //   if (!IsTrueGCDDomain(R))
  //     CoCoA_THROW_ERROR(ERR::NotTrueGCDDomain, FnName);
  //   RingElem gcd(R);
  //   // BUG: clobbers xcoeff and ycoeff always; use locals and swap into xcoeff, ycoeff???
  //   xcoeff = gcd; // must force into same ring
  //   ycoeff = gcd; // must force into same ring
  //   R->myExgcd(raw(gcd), raw(xcoeff), raw(ycoeff), raw(x), raw(y));
  //   return gcd;
  // }



  void GcdQuot(RingElem& gcd, RingElem& quot1, RingElem& quot2, ConstRefRingElem x, ConstRefRingElem y)
  {
    const char* const FnName = "GcdQuot(g,q1,q2,r1,r2)";
    const ring& R = owner(x);
    if (owner(y) != R)
      CoCoA_THROW_ERROR(ERR::MixedRings, FnName);
    if (!IsTrueGCDDomain(R))
      CoCoA_THROW_ERROR(ERR::NotTrueGCDDomain, FnName);
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
        CoCoA_THROW_ERROR(ERR::MixedRings, "RingElem += RingElem -- refuse to promote LHS");
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
        CoCoA_THROW_ERROR(ERR::MixedRings, "RingElem -= RingElem -- refuse to promote LHS");
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
        CoCoA_THROW_ERROR(ERR::MixedRings, "RingElem *= RingElem -- refuse to promote LHS");
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
      if (IsPolyRing(Rx) && CoeffRing(Rx) == Ry)  { if (IsZeroDivisor(y))  CoCoA_THROW_ERROR(ERR::DivByZero, "RingElem /= RingElem"); PolyRing(Rx)->myDivByCoeff(raw(x),raw(y)); return x; }
      const RingHom promote = AutomaticConversionHom(Rx,Ry,ErrMixed); // throws ErrMixed if not permitted 
      if (codomain(promote) == Ry)
        CoCoA_THROW_ERROR(ERR::MixedRings, "RingElem /= RingElem -- refuse to promote LHS");
      return x /= promote(y);
    }

    if (IsZeroDivisor(y))  CoCoA_THROW_ERROR(ERR::DivByZero, "RingElem /= RingElem");
    if (!Rx->myIsDivisible(raw(x), raw(x), raw(y)))
      CoCoA_THROW_ERROR(ERR::BadQuot, "RingElem /= RingElem");
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
    if (IsZero(x)) CoCoA_THROW_ERROR(ERR::ReqNonZeroRingElem, "IsIrred");
    if (IsField(owner(x))) CoCoA_THROW_ERROR(ERR::NotElemGCDDomain, "IsIrred in field");
    if (!IsTrueGCDDomain(owner(x))) CoCoA_THROW_ERROR(ERR::NotElemGCDDomain, "IsIrred for non TrueGCDDomain");
    if (IsInvertible(x)) CoCoA_THROW_ERROR(ERR::InvertibleRingElem, "IsIrred -- elem is invertible");
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
      if (!AllowFieldFlag && IsField(Rnum)) CoCoA_THROW_ERROR("Not allowed in a field", "IsDivisible");
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
      if (!AllowFieldFlag && IsField(Rnum)) CoCoA_THROW_ERROR("Not allowed in a field", "IsDivisible");
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
      CoCoA_THROW_ERROR(ERR::BadPwrZero, "power(RingElem, n)");
    if (IsNegative(n) && !IsInvertible(x))
      CoCoA_THROW_ERROR(ERR::NotUnit, "power(RingElem, n) and n < 0");
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
      CoCoA_THROW_ERROR(ERR::BadPwrZero, "power(RingElem, N)");
    if (N < 0 && !IsInvertible(x))
      CoCoA_THROW_ERROR(ERR::NotUnit, "power(RingElem, N) and N < 0");

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
      CoCoA_THROW_ERROR("Arg must be finite field", "FieldGen");
    return Fq->myFieldGen();
  }


  bool IsPthPower(ConstRefRingElem x)  ///< only in Fp and Fp[x,y,...]
  {
    const BigInt P = characteristic(owner(x));
    long p;
    if (!IsConvertible(p, P) || !IsPrime(p))
      CoCoA_THROW_ERROR(ERR::BadArg, "IsPthPower");
    if (IsFiniteField(owner(x))) return true;
    if (IsFractionField(owner(x))) return IsPthPower(num(x)) && IsPthPower(den(x));
    if (!IsPolyRing(owner(x)))
      CoCoA_THROW_ERROR(ERR::BadArg, "IsPthPower");
    if (!IsSparsePolyRing(owner(x))) CoCoA_THROW_ERROR(ERR::NYI, "IsPthPower -- not a sparse poly ring");
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
      CoCoA_THROW_ERROR(ERR::ArgTooBig, "PthRoot");  // BUG -- weak impl!
    if (IsFractionField(owner(x)))
    { // BUG BUG IMPERFECT IMPL -- could fail if num & den have non-trivial invertible factors without PthRoots
      const RingElem N = PthRoot(num(x));
      const RingElem D = PthRoot(den(x));
      const RingHom phi = EmbeddingHom(owner(x));
      return phi(N)/phi(D);
    }
    if (!IsPolyRing(owner(x)))
      CoCoA_THROW_ERROR(ERR::BadArg, "PthRoot");
    if (!IsSparsePolyRing(owner(x))) CoCoA_THROW_ERROR(ERR::NYI, "PthRoot -- not a sparse poly ring");
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
      CoCoA_THROW_ERROR(ERR::BadRing, "radical");
    if (IsZZ(R))
      return RingElem(R, radical(ConvertTo<BigInt>(x)));
    return RadicalOfPoly(x);
  }



  RingElem binomial(RingElem x, const MachineInt& n)
  {
    if (IsNegative(n)) CoCoA_THROW_ERROR(ERR::ReqNonNegative, "binomial(RingElem,int) <-- arg 2");
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
    if (IsZeroDivisor(ans)) CoCoA_THROW_ERROR(ERR::DivByZero, "RingElem / MachineInt");
    if (!owner(r)->myIsDivisible(raw(ans), raw(r), raw(ans)))
      CoCoA_THROW_ERROR(ERR::BadQuot, "RingElem / MachineInt");
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
    if (IsZeroDivisor(r))  CoCoA_THROW_ERROR(ERR::DivByZero, "MachineInt / RingElem");
    RingElem ans(owner(r), n); // use ans to convert n to a RingElem
    if (!owner(r)->myIsDivisible(raw(ans), raw(ans), raw(r)))
      CoCoA_THROW_ERROR(ERR::BadQuot, "MachineInt / RingElem");
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
    if (IsZeroDivisor(den)) CoCoA_THROW_ERROR(ERR::DivByZero, "RingElem /= MachineInt");
    if (!owner(r)->myIsDivisible(raw(r), raw(r), raw(den)))
      CoCoA_THROW_ERROR(ERR::BadQuot, "RingElem /= MachineInt");
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
    if (IsZeroDivisor(ans))  CoCoA_THROW_ERROR(ERR::DivByZero, "RingElem / BigInt");
    if (!owner(r)->myIsDivisible(raw(ans), raw(r), raw(ans)))
      CoCoA_THROW_ERROR(ERR::BadQuot, "RingElem / BigInt");
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
    if (IsZeroDivisor(r)) CoCoA_THROW_ERROR(ERR::DivByZero, "BigInt / RingElem");
    RingElem ans(owner(r), N);
    if (!owner(r)->myIsDivisible(raw(ans), raw(ans), raw(r)))
      CoCoA_THROW_ERROR(ERR::BadQuot, "BigInt / RingElem");
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
    if (IsZeroDivisor(den)) CoCoA_THROW_ERROR(ERR::DivByZero, "RingElem /= BigInt");
    if (!owner(r)->myIsDivisible(raw(r), raw(r), raw(den)))
      CoCoA_THROW_ERROR(ERR::BadQuot, "RingElem /= BigInt");
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
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "sign(RingElem)");
    return owner(x)->mySign(raw(x));
  }


  RingElem abs(ConstRefRingElem x)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "abs(RingElem)");
    if (x < 0) return -x;
    return x;
  }


  BigInt floor(ConstRefRingElem x)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "floor(RingElem)");
    return owner(x)->myFloor(raw(x));
  }


  BigInt ceil(ConstRefRingElem x)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "ceil(RingElem)");
    return owner(x)->myCeil(raw(x));
  }


  BigInt NearestInt(ConstRefRingElem x)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "NearestInt(RingElem)");
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
      CoCoA_THROW_ERROR(ERR::MixedRings, "cmp(RingElem, RingElem)");
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "cmp(RingElem, RingElem)");
    return owner(x)->myCmp(raw(x), raw(y));
  }


  int CmpAbs(ConstRefRingElem x, ConstRefRingElem y)
  {
    if (owner(x) != owner(y))
      CoCoA_THROW_ERROR(ERR::MixedRings, "CmpAbs(RingElem, RingElem)");
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "CmpAbs(RingElem, RingElem)");
    return owner(x)->myCmpAbs(raw(x), raw(y));
  }


  bool operator<(ConstRefRingElem x, ConstRefRingElem y)
  {
    if (owner(x) != owner(y))
      CoCoA_THROW_ERROR(ERR::MixedRings, "RingElem < RingElem");
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "RingElem < RingElem");
    return owner(x)->myCmp(raw(x), raw(y)) < 0;
  }


  bool operator<=(ConstRefRingElem x, ConstRefRingElem y)
  {
    if (owner(x) != owner(y))
      CoCoA_THROW_ERROR(ERR::MixedRings, "RingElem <= RingElem");
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "RingElem <= RingElem");
    return owner(x)->myCmp(raw(x), raw(y)) <= 0;
  }


  bool operator>(ConstRefRingElem x, ConstRefRingElem y)
  {
    if (owner(x) != owner(y))
      CoCoA_THROW_ERROR(ERR::MixedRings, "RingElem > RingElem");
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "RingElem > RingElem");
    return owner(x)->myCmp(raw(x), raw(y)) > 0;
  }


  bool operator>=(ConstRefRingElem x, ConstRefRingElem y)
  {
    if (owner(x) != owner(y))
      CoCoA_THROW_ERROR(ERR::MixedRings, "RingElem >= RingElem");
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "RingElem >= RingElem");
    return owner(x)->myCmp(raw(x), raw(y)) >= 0;
  }


  /////////////////////////////////////////////////////////////////////////////
  // First arg is a MachineInt, second is a RingElem

  int cmp(const MachineInt& x, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "cmp(MachineInt, RingElem)");
    if (IsZero(x)) return -sign(y);
    return owner(y)->myCmp(raw(RingElem(owner(y),x)), raw(y));
  }


  bool operator<(const MachineInt& x, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "MachineInt < RingElem");
    if (IsZero(x)) return sign(y) > 0;
    return owner(y)->myCmp(raw(RingElem(owner(y),x)), raw(y)) < 0;
  }


  bool operator<=(const MachineInt& x, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "MachineInt <= RingElem");
    if (IsZero(x)) return sign(y) >= 0;
    return owner(y)->myCmp(raw(RingElem(owner(y),x)), raw(y)) <= 0;
  }


  bool operator>(const MachineInt& x, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "MachineInt > RingElem");
    if (IsZero(x)) return sign(y) < 0;
    return owner(y)->myCmp(raw(RingElem(owner(y),x)), raw(y)) > 0;
  }


  bool operator>=(const MachineInt& x, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "MachineInt >= RingElem");
    if (IsZero(x)) return sign(y) <= 0;
    return owner(y)->myCmp(raw(RingElem(owner(y),x)), raw(y)) >= 0;
  }


  /////////////////////////////////////////////////////////////////////////////
  // First arg is a RingElem, second is a MachineInt


  int cmp(ConstRefRingElem x, const MachineInt& y)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "cmp(RingElem, MachineInt)");
    if (IsZero(y)) return sign(x);
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),y)));
  }


  bool operator<(ConstRefRingElem x, const MachineInt& y)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "RingElem < MachineInt");
    if (IsZero(y)) return sign(x) < 0;
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),y))) < 0;
  }


  bool operator<=(ConstRefRingElem x, const MachineInt& y)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "RingElem <= MachineInt");
    if (IsZero(y)) return sign(x) <= 0;
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),y))) <= 0;
  }


  bool operator>(ConstRefRingElem x, const MachineInt& y)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "RingElem > MachineInt");
    if (IsZero(y)) return sign(x) > 0;
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),y))) > 0;
  }


  bool operator>=(ConstRefRingElem x, const MachineInt& y)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "RingElem >= MachineInt");
    if (IsZero(y)) return sign(x) >= 0;
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),y))) >= 0;
  }

  ///////////////////////////////////////////////////////////////////////////
  // Comparisons between BigInt and RingElem
  int cmp(const BigInt& x, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "cmp(BigInt, RingElem)");
    return owner(y)->myCmp(raw(RingElem(owner(y),x)), raw(y));
  }


  bool operator<(const BigInt& x, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "BigInt < RingElem");
    return owner(y)->myCmp(raw(RingElem(owner(y),x)), raw(y)) < 0;
  }


  bool operator<=(const BigInt& x, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "BigInt <= RingElem");
    return owner(y)->myCmp(raw(RingElem(owner(y),x)), raw(y)) <= 0;
  }


  bool operator>(const BigInt& x, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "BigInt > RingElem");
    return owner(y)->myCmp(raw(RingElem(owner(y),x)), raw(y)) > 0;
  }


  bool operator>=(const BigInt& x, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "BigInt >= RingElem");
    return owner(y)->myCmp(raw(RingElem(owner(y),x)), raw(y)) >= 0;
  }


  /////////////////////////////////////////////////////////////////////////////
  // First arg is a RingElem, second is a const BigInt&


  int cmp(ConstRefRingElem x, const BigInt& y)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "cmp(RingElem, BigInt)");
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),y)));
  }


  bool operator<(ConstRefRingElem x, const BigInt& y)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "RingElem < BigInt");
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),y))) < 0;
  }


  bool operator<=(ConstRefRingElem x, const BigInt& y)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "RingElem <= BigInt");
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),y))) <= 0;
  }


  bool operator>(ConstRefRingElem x, const BigInt& y)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "RingElem > BigInt");
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),y))) > 0;
  }


  bool operator>=(ConstRefRingElem x, const BigInt& y)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "RingElem >= BigInt");
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
    if (IsZeroDivisor(ans)) CoCoA_THROW_ERROR(ERR::DivByZero, "RingElem / BigRat");
    if (!owner(r)->myIsDivisible(raw(ans), raw(r), raw(ans)))
      CoCoA_THROW_ERROR(ERR::BadQuot, "RingElem / BigRat");
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
    if (IsZeroDivisor(r)) CoCoA_THROW_ERROR(ERR::DivByZero, "BigRat / RingElem");
    RingElem ans(owner(r), Q);
    if (!owner(r)->myIsDivisible(raw(ans), raw(ans), raw(r)))
      CoCoA_THROW_ERROR(ERR::BadQuot, "BigRat / RingElem");
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
    if (IsZeroDivisor(den)) CoCoA_THROW_ERROR(ERR::DivByZero, "RingElem /= BigRat");
    if (!owner(r)->myIsDivisible(raw(r), raw(r), raw(den)))
      CoCoA_THROW_ERROR(ERR::BadQuot, "RingElem /= BigRat");
    return r;
  }


  ///////////////////////////////////////////////////////////////////////////
  // Comparisons between BigRat and RingElem
  int cmp(const BigRat& Q, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "cmp(BigRat, RingElem)");
    return owner(y)->myCmp(raw(RingElem(owner(y),Q)), raw(y));
  }


  bool operator<(const BigRat& Q, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "BigRat < RingElem");
    return owner(y)->myCmp(raw(RingElem(owner(y),Q)), raw(y)) < 0;
  }


  bool operator<=(const BigRat& Q, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "BigRat <= RingElem");
    return owner(y)->myCmp(raw(RingElem(owner(y),Q)), raw(y)) <= 0;
  }


  bool operator>(const BigRat& Q, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "BigRat > RingElem");
    return owner(y)->myCmp(raw(RingElem(owner(y),Q)), raw(y)) > 0;
  }


  bool operator>=(const BigRat& Q, ConstRefRingElem y)
  {
    if (!IsOrderedDomain(owner(y)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "BigRat >= RingElem");
    return owner(y)->myCmp(raw(RingElem(owner(y),Q)), raw(y)) >= 0;
  }


  /////////////////////////////////////////////////////////////////////////////
  // First arg is a RingElem, second is a BigRat


  int cmp(ConstRefRingElem x, const BigRat& Q)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "cmp(RingElem, BigRat)");
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),Q)));
  }


  bool operator<(ConstRefRingElem x, const BigRat& Q)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "RingElem < BigRat");
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),Q))) < 0;
  }


  bool operator<=(ConstRefRingElem x, const BigRat& Q)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "RingElem <= BigRat");
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),Q))) <= 0;
  }


  bool operator>(ConstRefRingElem x, const BigRat& Q)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "RingElem > BigRat");
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),Q))) > 0;
  }


  bool operator>=(ConstRefRingElem x, const BigRat& Q)
  {
    if (!IsOrderedDomain(owner(x)))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "RingElem >= BigRat");
    return owner(x)->myCmp(raw(x), raw(RingElem(owner(x),Q))) >= 0;
  }




} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/ring.C,v 1.111 2024/07/15 17:24:45 bigatti Exp $
// $Log: ring.C,v $
// Revision 1.111  2024/07/15 17:24:45  bigatti
// Summary: Error codes:  added CoCoA_THROW_ERROR1(ErrID);
// Removed BadRowIndex, BadColIndex;  Renamed  ERR::Empty --> ERR::ReqNonEmpty
//
// Revision 1.110  2024/07/02 20:42:07  abbott
// Summary: Renaming of errors (redmine 308)
//
// Revision 1.109  2024/07/02 15:37:50  bigatti
// Summary: Changed error codes: LogZero into ReqNonZero
// and Not... into ReqNonNegative, ReqNonNegativeGrading, ReqPositive, ReqPositiveGrading
//
// Revision 1.108  2024/07/02 12:22:14  bigatti
// Summary: error codes, first batch:
// ReqUnivariate, ReqNonZero, ReqNonZeroGradingDim, ReqNonZeroModulus,  ReqNonZeroRingElem
//
// Revision 1.107  2024/06/27 19:39:20  abbott
// Summary: Added short-cuts in op*= and op/= for mult/div poly by a coeff
//
// Revision 1.106  2024/03/18 08:48:54  abbott
// Summary: Removed cruft; improved indentation
//
// Revision 1.105  2023/03/09 22:35:16  abbott
// Summary: Corrected some err mesgs
//
// Revision 1.104  2023/03/08 19:43:57  abbott
// Summary: Made myCharacteristic a function (instead of proc with OUT param) redmine 536
//
// Revision 1.103  2023/03/06 22:08:08  abbott
// Summary: Added variants of IsCoprime (mixed RingElem, integer)
//
// Revision 1.102  2023/01/31 12:30:18  abbott
// Summary: Moved GcdList to here (just not to lose it; currently never called)
//
// Revision 1.101  2022/11/30 20:37:50  abbott
// Summary: Restricted IsCoprime to ZZ or k[x,y,...] (redmine 1710)
//
// Revision 1.100  2022/10/25 17:29:18  abbott
// Summary: Added FieldGen (finite field) and mem fn myFieldGen (default defn gives error)
//
// Revision 1.99  2022/10/10 10:00:16  abbott
// Summary: Added move-assign op
//
// Revision 1.98  2022/08/05 13:13:28  abbott
// Summary: Auto promotion (only on RHS) for +=, -=, *=, /=
//
// Revision 1.97  2022/06/08 15:42:20  abbott
// Summary: Swap fn made inline in header file
//
// Revision 1.96  2022/02/18 14:12:02  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.95  2022/02/08 20:18:55  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.94  2021/08/03 12:24:45  abbott
// Summary: Better err mesg
//
// Revision 1.93  2021/02/17 17:51:40  abbott
// Summary: myLogCardinality now throws by default (i.e. if not a finite field); redmine 898
//
// Revision 1.92  2021/01/07 15:23:38  abbott
// Summary: Corrected copyright
//
// Revision 1.91  2020/10/27 20:02:00  abbott
// Summary: GcdQuot now sets correct ring for OUT params (redmine 1085)
//
// Revision 1.90  2020/10/14 20:03:07  abbott
// Summary: Implemented IsDivisible_AllowFields (redmine 1500)
//
// Revision 1.89  2020/10/05 19:28:44  abbott
// Summary: Added comment in power (maybe handle simple cases, i.e. power of 1,0 or -1)
//
// Revision 1.88  2020/06/22 15:43:37  abbott
// Summary: Use new CoCoA_STATIC_ERROR_MESG macro
//
// Revision 1.87  2020/06/20 19:12:17  abbott
// Summary: AutomaticConversionHom now requires 3rd arg (FnName of caller)
//
// Revision 1.86  2020/06/19 15:06:57  abbott
// Summary: Added code to allow dyadic arith ops in mixed rings
//
// Revision 1.85  2020/06/17 15:49:30  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.84  2020/06/17 15:21:51  bigatti
// -- added HMul for multiplication across rings
//
// Revision 1.83  2020/02/11 16:56:42  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.82  2020/02/11 16:12:20  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.81  2019/03/27 13:51:07  bigatti
// (abbott) removed ThreadsafeCounter --> using atomic
//
// Revision 1.80  2019/03/18 13:16:26  abbott
// Summary: Added new include (after splitting NumTheory); commented out exgcd.
//
// Revision 1.79  2018/05/18 16:39:39  bigatti
// -- added include SparsePolyOps-RingElem.H
//
// Revision 1.78  2018/05/18 12:25:54  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.77  2018/05/17 15:58:24  bigatti
// -- added include SparsePolyIter
//
// Revision 1.76  2018/04/20 18:51:26  abbott
// Summary: Changed ctors for BigInt/BigRat from string or from MPZ/MPQ
//
// Revision 1.75  2018/02/27 17:30:23  abbott
// Summary: Renamed NumTheory_prime to NumTheory-prime; changed includes
//
// Revision 1.74  2018/02/27 10:56:07  abbott
// Summary: Added include NumTheory_prime
//
// Revision 1.73  2017/11/14 15:14:18  abbott
// Summary: Changed defn of IsIrred: now gives error when it should
//
// Revision 1.72  2017/10/17 15:44:26  abbott
// Summary: Added new fn IsCoprime
//
// Revision 1.71  2017/09/06 11:56:30  abbott
// Summary: Changed ERR::SERIOUS into ERR::ShouldNeverGetHere
//
// Revision 1.70  2017/04/27 14:03:50  bigatti
// -- added RingElem(ring, string)
//
// Revision 1.69  2016/11/11 14:15:34  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.68  2016/10/20 18:04:44  abbott
// Summary: Added radical
//
// Revision 1.67  2016/03/25 20:13:09  abbott
// Summary: Introduced new mem fns myCeil, myNearestInt; renamed NearestInteger to NearestInt.  New impls for floor, ceil, NearestInt
//
// Revision 1.66  2015/06/26 15:02:38  abbott
// Summary: Added context arg to call to CheckForInterrupt
// Author: JAA
//
// Revision 1.65  2015/04/30 08:47:31  bigatti
// -- added myIsZeroDivisor (and changed IsZeroDivisor accordingly)
//
// Revision 1.64  2014/07/11 15:38:09  bigatti
// -- added myOutputSelfShort and myOutputSelfLong
//
// Revision 1.63  2014/07/08 08:40:36  abbott
// Summary: Removed AsFractionField
// Author: JAA
//
// Revision 1.62  2014/07/07 13:28:41  abbott
// Summary: Removed AsPolyRing; Removed AsSparsePolyRing
// Author: JAA
//
// Revision 1.61  2014/06/14 19:45:08  abbott
// Summary: Added new fn CmpAbs for RingElem (via memfn myCmpAbs)
// Author: JAA
//
// Revision 1.60  2014/05/16 12:27:01  abbott
// Summary: Changed defn of NearestInteger to be compatible with round/RoundDiv
// Author: JAA
//
// Revision 1.59  2014/04/30 16:28:53  abbott
// Summary: Replaced size_t by long
// Author: JAA
//
// Revision 1.58  2014/04/08 15:41:28  abbott
// Summary: Replaced test for IsZero by IsZeroDivisor (in divisions)
// Author: JAA
//
// Revision 1.57  2014/03/27 17:18:43  abbott
// Summary: Improved IsZeroDivisor (now uses IsIntegralDomain3)
// Author: JAA
//
// Revision 1.56  2014/03/26 17:05:59  abbott
// Summary: binomial(RingElem,int) now delegates to binomial(int,int) when it can
// Author: JAA
//
// Revision 1.55  2014/03/26 16:29:27  abbott
// Summary: Added binomial(RingElem,int)
// Author: JAA
//
// Revision 1.54  2014/03/06 17:02:57  abbott
// Summary: Removed some incorrect assertions (which complained about 0^0)
// Author: JAA
//
// Revision 1.53  2014/03/06 15:50:19  abbott
// Summary: Zero to power zero now gives 1 (previously it was inconsistent)
// Author: JAA
//
// Revision 1.52  2014/01/29 13:03:30  abbott
// Summary: Added ctors for RingElem from mpz_t and mpq_t.  Improved impls for IsPthPower and PthRoot
// Author: John Abbott
//
// Revision 1.51  2014/01/28 11:11:09  bigatti
// -- added  IsDivisible(RingElem& lhs, A, B);  with all(?) variants
//
// Revision 1.50  2013/10/15 16:24:45  abbott
// Improved impl of PthRoot so that it works for large finite fields.
//
// Revision 1.49  2013/06/28 17:03:51  abbott
// Modified semantics of IdealBase::myDivMod;
// it now returns a boolean.
// Several consequential changes.
//
// Revision 1.48  2013/05/29 17:03:03  bigatti
// -- added IsZeroDivisor (to be improved)
//
// Revision 1.47  2013/05/08 12:14:27  abbott
// Removed include directive which is no longer needed.
//
// Revision 1.46  2013/03/26 14:57:03  abbott
// Removed use of obsolete "convert" procedure.
//
// Revision 1.45  2013/03/15 14:58:01  bigatti
// -- added   RingElem::RingElem(const ring& R, ConstRefRingElem rhs):
//
// Revision 1.44  2013/02/11 11:45:24  abbott
// Fixed embarrassing bug: IsNegative(n) < 0
//
// Revision 1.43  2012/10/24 13:02:34  abbott
// Renamed class ConstRefRingElem to RingElemAlias.
// ConstRefRingElem is now a typedef.  Numerous consequential changes.
//
// Revision 1.42  2012/10/17 12:23:29  abbott
// Removed class RefRingElem;  replaced  RefRingElem  by RingElem&.
//
// Revision 1.41  2012/10/03 15:23:25  abbott
// Added default ctor for ring (produces RingZZ).
// NB assignment of rings was already working, apparently!
// Added default ctor for RingElem (produces 0 in RingZZ).
// Added "ring changing" assignment for RingElem.
// Added "ring changing" swap for RingElem.
//
// Revision 1.40  2012/06/14 14:42:24  abbott
// Added lcm for RingElem and MachineInt/BigInt --  were missing.
//
// Revision 1.39  2012/05/30 16:09:11  bigatti
// -- applied "3" convention on bool3
//
// Revision 1.38  2012/05/29 14:56:15  abbott
// Separated ThreadsafeCounter from symbol; also employed it in ring.C.
//
// Revision 1.37  2012/05/29 07:45:23  abbott
// Implemented simplification change to bool3:
//  changed names of the constants,
//  changed names of the testing fns.
//
// Revision 1.36  2012/05/28 10:35:32  abbott
// Changed default defn of IsTrueGCDDomain (makes RingQQImpl a bit simpler).
//
// Revision 1.35  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.34  2012/05/24 13:05:43  abbott
// RingBase::myGcd (default defn) now throws NotTrueGCDDomain error.
//
// Revision 1.33  2012/05/22 10:02:37  abbott
// Removed IsGCDDomain; substituted by IsTrueGCDDomain.
// Added IsFractionFieldOfGCDDomain.
//
// Revision 1.32  2012/04/27 15:01:20  abbott
// Added several fns:
//   IsFiniteField, LogCardinality, IsPthPower, PthRoot
//   IsPID, IsPID3
//
// Revision 1.31  2012/02/10 10:34:09  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.30  2012/02/08 13:48:39  bigatti
// -- changed Z,Q --> ZZ,QQ
//
// Revision 1.29  2011/11/09 14:29:37  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.28  2011/08/24 10:32:04  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.27  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.26  2011/08/12 15:26:54  abbott
// Improved pseudo-ctor for RingElem from a BigRat: better exception safety.
// Modified RingBase::myIsZeroAddMul -- would be much cleaner if args
// were RingElems rather than raw ptrs.
//
// Revision 1.25  2011/06/23 16:04:46  abbott
// Added IamExact mem fn for rings.
// Added myRecvTwinFloat mem fn for rings.
// Added first imple of RingHom from RingTwinFloat to other rings.
//
// Revision 1.24  2011/05/19 14:50:13  abbott
// Removed defns of old form of (my)IsRational.
// Added default defn of myIsDouble.
// Added defn of IsDouble.
// Removed old commented out impl of printing of std::vector<RingElem>.
//
// Revision 1.23  2011/03/30 09:13:04  bigatti
// -- added myNew(const symbol& s)  and fixed  RingElem ctor
//
// Revision 1.22  2011/03/16 15:17:20  bigatti
// -- fixed a CoCoA_ASSERT
//
// Revision 1.21  2011/03/11 21:50:59  abbott
// Fixed an off-by-one error in myBinaryPower.
//
// Revision 1.20  2011/03/10 16:39:33  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.19  2011/02/28 14:20:44  bigatti
// -- just some comment
//
// Revision 1.18  2011/01/28 11:13:56  bigatti
// -- moved default impl for IsPrintedWithMinus to ring.H
//
// Revision 1.17  2011/01/19 16:34:30  bigatti
// -- added lcm/myLcm
//
// Revision 1.16  2010/12/20 15:19:29  bigatti
// -- modified IsZeroAddMul with temporary variable (slug found with cyclotomic)
//
// Revision 1.15  2010/10/01 15:45:39  bigatti
// -- added mySymbolValue
// -- added RingElem(R, sym)
//
// Revision 1.14  2009/12/11 11:46:32  abbott
// Changed fn  convert  into  IsConvertible.
// Added template procedure  convert.
// New version because change is not backward compatible.
//
// Revision 1.13  2009/12/03 17:37:22  abbott
// Moved fns characteristic & swap into .C file, so that .H does not
// need to include ZZ.H and error.H.
//
// Revision 1.12  2009/10/29 18:30:09  abbott
// Changed the ring ID value to be a long (previously was size_t).
// Some minor cleaning to include directives.
//
// Revision 1.11  2009/07/02 16:32:10  abbott
// Consequential changes stemming from new class BigRat, and modified interface to the member
// function RingBase::myIsRational.  Also some new conversion functions.
//
// Revision 1.10  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.9  2008/11/18 10:24:48  abbott
// Added floor, ceil, NearestInteger and CmpDouble functions.
// Also added myFloor member fn.  Plus some very minor code tidying.
//
// Revision 1.8  2008/04/22 14:44:26  abbott
// Removed top level function square.
// Internal powering routine myBinaryPower now uses mySquare.
//
// Revision 1.7  2008/04/22 13:09:15  abbott
// Removed IsPositive and IsNegative functions for ZZs.
// Comparison between RingElem and 0 is now handled specially (specially fast?).
//
// Revision 1.6  2008/04/15 14:58:06  bigatti
// -- added mySquare, square
//
// Revision 1.5  2008/03/12 16:29:37  bigatti
// -- added: IsIrred
//
// Revision 1.4  2007/10/30 17:14:06  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.3  2007/05/22 22:45:14  abbott
// Changed fn name IsUnit to IsInvertible.
//
// Revision 1.2  2007/03/23 18:38:42  abbott
// Separated the "convert" functions (and CheckedCast) into their own files.
// Many consequential changes.  Also corrected conversion to doubles.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.9  2007/03/08 18:22:28  cocoa
// Just whitespace cleaning.
//
// Revision 1.8  2006/12/07 17:23:46  cocoa
// -- for compilation with _Wextra: commented out names of unused arguments
//
// Revision 1.7  2006/12/06 17:40:10  cocoa
// -- style: RawPtr args are now called "raw.."
//
// Revision 1.6  2006/11/20 15:55:02  cocoa
// ring is now a class again.  Improved definitions of operator-> in derived classes.
//
// Revision 1.5  2006/11/02 13:25:43  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.4  2006/10/27 19:09:45  cocoa
// Replaced some member functions of CoCoA::symbol by friend functions.
// Removed some include dependency on symbol.H
//
// Revision 1.3  2006/10/16 23:18:58  cocoa
// Corrected use of std::swap and various special swap functions.
// Improved myApply memfn for homs of RingDistrMPolyInlPP.
//
// Revision 1.2  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.11  2006/05/29 16:22:37  cocoa
// Third time lucky???
// Added myIsInteger member function to all rings (NYI for RingFloat).
//
// Revision 1.10  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.9  2006/04/27 13:43:35  cocoa
// Added abs function (only for arithmetically ordered rings).
//
// Revision 1.8  2006/04/21 15:01:36  cocoa
// Changed default implementation of RingBase::myGcd -- it now gives a SERIOUS
// error.  All fields must now handle a call to gcd explicitly: they can use
// the new myGcdInField function.  It's now cleaner than it was.
//
// Revision 1.7  2006/04/04 17:20:53  cocoa
// -- fixed: myPower trivial case n==0
//
// Revision 1.6  2006/03/27 12:21:25  cocoa
// Minor silly changes to reduce number of complaints from some compiler or other.
//
// Revision 1.5  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.4  2006/03/14 15:01:49  cocoa
// Improved the implementation of ring member fns for computing powers.
// Should keep Intel C++ compiler quieter too.
//
// Revision 1.3  2006/03/12 21:28:33  cocoa
// Major check in after many changes
//
// Revision 1.2  2006/02/14 16:22:20  cocoa
// -- defined "operator<<" for vector<RingElem>&  in ring.H/C
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.4  2005/09/22 18:04:17  cocoa
// It compiles; the tests run OK.  The examples compile.
// No documentation -- the mindless eurocrats have rendered
// me mindless too.
//
// Revision 1.3  2005/07/08 15:09:28  cocoa
// Added new symbol class (to represent names of indets).
// Integrated the new class into concrete polynomial rings
// and PPMonoid -- many consequential changes.
// Change ctors for the "inline" sparse poly rings: they no
// longer expect a PPMonoid, but build their own instead
// (has to be a PPMonoidOv).
//
// Revision 1.2  2005/06/30 16:02:18  cocoa
// -- myIsPrintedWithMinus now returns false (used to throw)
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.6  2005/04/27 16:14:56  cocoa
// Cleaned up example programs -- added "free use" permit.
// Changed a couple of ErrorInfo object names, and added
// ERR::NotTrueGCDDomain.
//
// Revision 1.5  2005/04/20 15:40:47  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.4  2005/04/19 14:06:03  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.3  2005/02/11 16:45:24  cocoa
// Removed the useless and misleading functions myInit and myKill
// from the SmallFp*Impl classes; various consequential changes.
//
// Revision 1.2  2005/02/11 14:15:20  cocoa
// New style ring elements and references to ring elements;
// I hope I have finally got it right!
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.17  2004/11/18 18:33:40  cocoa
// Now every ring know its own "one" element (as well as "zero").
// Several consequential changes.
//
// Revision 1.16  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.15  2004/11/11 14:05:02  cocoa
// -- minor changes for doxygen
//
// Revision 1.14  2004/11/11 11:56:09  cocoa
// (1) Tidied makefiles, and introduced common.mki
// (2) Improved several tests, and cleaned them so that they
//     handle sanely any otherwise unhandled exceptions.
//
// Revision 1.13  2004/11/04 18:47:42  cocoa
// (1) Ring member functions which previously expected mpz_t args
//     now expect ZZ args.  Numerous minor consequential changes.
// (2) Renamed function which gives access to the mpz_t value inside
//     a ZZ object: previously was raw(...), now is mpzref(...).
//     Plenty of calls had to be altered.
//
// Revision 1.12  2004/07/27 16:03:38  cocoa
// Added IsCommutative test and IamCommutative member function
// to all rings.  Tidied geobuckets a little.
//
// Revision 1.11  2004/07/20 15:04:06  cocoa
// The next step in the new "ring element" conversion process:
// handling the case of creating a "const RefRingElem" object
// (since C++ refuses to do this properly itself).
//
// Revision 1.10  2004/07/16 15:45:12  cocoa
// First stage of new RingElem implementation completed.
//
// Revision 1.9  2004/05/27 16:14:02  cocoa
// Minor revision for new coding conventions.
//
// Revision 1.8  2004/05/24 15:52:13  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.7  2004/04/08 15:33:34  cocoa
// Added function IsInteger, and the related RingBase::myIsInteger
// virtual function, plus all necessary implementations.
//
// Revision 1.6  2004/03/20 17:46:10  cocoa
// Check in prior to departure to RWCA
//
// Revision 1.5  2004/02/03 16:16:20  cocoa
// Removed pointless IamGCDDomain functions from several concrete rings.
// Added IamOrderedDomain functions where appropriate.
// Tidied ctors for the small finite fields.
//
// Revision 1.4  2003/11/14 13:06:05  cocoa
// -- New function "myIsPrintAtom" for printing polynomials and fractions
//
// Revision 1.3  2003/10/17 10:51:06  cocoa
// Major cleaning, and new naming convention.
//
// Revision 1.2  2003/10/09 12:13:44  cocoa
// New coding convention for rings.
//
// Revision 1.18  2003/06/23 16:36:46  abbott
// Prior to public release:
//   fixed wrongly placed zero test in operator/(ConstRefRingElem,long)
//
// Revision 1.17  2003/04/16 14:16:46  abbott
// Added code to handle some trivial cases in SequentialPower and
// BinaryPower.  Added definition of new dispatch function IsDivisible.
//
// Revision 1.16  2003/01/21 15:38:48  abbott
// Many minor syntactic modifications consequent from the new
// names and structure in ring.H.
//
// Revision 1.15  2002/12/13 15:16:03  abbott
// Added functions IsOne and IsMinusOne (both at global level and
// as member functions of an AbstractRing).  Commented out the
// member function for testing equality between a RingElem and a
// machine integer.  Added member functions for initializing or
// assigning from a big integer (of type ZZ).  Added syntactic
// sugar operators for arithmetic between RingElems and ZZs.
//
// Revision 1.14  2002/12/11 11:49:18  abbott
// Checking in for safety -- this code does at least compile and work.
// Exponents in power functions must now be non-negative (unsigned int).
// Added functions between RingElems and ZZs (big integers).
// Added two member functions for initialising/assigning from big
// integers.
//
// Made op == and != inline because Anna wanted them so; probably revoke
// this after IsOne has been invented.
//
// Revision 1.13  2002/11/14 16:25:32  abbott
// Extensive consequential changes due to renaming in ring.H.
// ring renamed to AbstractRing.
// New typedef RingElem for AbstractRing::elem.
// New typedef ConstRawValue to be used as arg type in ring member functions.
// Old alias subclass has been largely replaced by new ConstRefElem subclass.
// New typedef ConstRefRingElem to be used as arg type for read only RingElems :-(
//
// Revision 1.12  2002/07/05 15:18:17  abbott
// Added definition of member fn SequentialPower (previously simply called power).
// Redefined op/ so that it calls IsDivisible, and gives a useful error mesg
// if the quotient does not exist in the ring.
//
// Revision 1.11  2002/06/27 16:06:34  abbott
// Added BinaryPower and BinaryPowerLoop default definitions.
//
// Revision 1.10  2002/06/21 14:40:16  abbott
// Rewrote the default implementation of IsZeroAddMul; now it is
// cleaner, and possibly more exception safe.
// Added definitions of the "syntactic sugar" arithmetic operations
// (between ring::elems and longs).
//
// Revision 1.9  2002/05/30 13:34:14  abbott
// Moved equals and not-equals operators to ring.H so they can be inline.
//
// Revision 1.8  2002/05/15 14:55:26  abbott
// Added power function default definition.
// Added code to handle alias ring::elem objects.
// General cleaning, and reorganization of indentation.
//
// Revision 1.7  2002/02/15 11:55:11  bigatti
// - added IsGCDDomain, gcd, "=", "==", "!="
//
// Revision 1.6  2002/01/30 15:00:14  abbott
// Tidied the "using" declarations.
// Tidied the "unimplemented: mixed ring arith" error messages.
// Removed the definitions of the two "raw" functions to ring.H.
//
// Revision 1.5  2001/12/07 18:20:40  abbott
// Changed names in accordance with new coding conventions
// (in particular ring::raw_elem became ring::RawValue).
//
// Revision 1.4  2001/11/16 18:45:06  bigatti
// added:   using namespace std;
// for compatibility with gcc-3
//
// Revision 1.3  2001/11/07 20:52:53  abbott
// Added implementations of the functions elem::raw(...).
// These could be inline, but are not yet.  I'll wait until I'm
// sure they're needed before thinking too much about it.
//
// Revision 1.2  2001/10/29 20:34:38  abbott
// Several changes: numerous minor fixes, and the more important correction
// of the standard binary arithmetic operators so that they no longer leak
// heap copies of every value computed (red face time).  The simplistic
// fixes rather imply a strong need for reference counting (or some other
// way of achieving cheap copying) of elem values -- to be discussed!
//
// Revision 1.1  2001/10/05 12:44:31  abbott
// Initial revision
//
