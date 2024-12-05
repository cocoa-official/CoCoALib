//   Copyright (c)  2005-2018  John Abbott and Anna M. Bigatti

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


// Source code for classes ideal and IdealBase

#include "CoCoA/ideal.H"

#include "CoCoA/OpenMath.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/VectorOps.H"  // for HasUniqueOwner
#include "CoCoA/utils.H"  // for len

#include <algorithm>
using std::copy;
#include <iostream>
using std::ostream;
//#include <vector>
using std::vector;


namespace CoCoA
{

  // C++ needs this function to be defined
  IdealBase::~IdealBase()
  {}


  //---------------------------------------------------------------------------

  ideal::ideal(IdealBase* IPtr):
      myPtr(IPtr)
  {
    IPtr->myRefCountInc();
  }


  ideal::ideal(ConstRefRingElem r)
  {
    vector<RingElem> gens;
    gens.push_back(RingElem(r));
    ideal tmp = owner(r)->myIdealCtor(gens);
    myPtr = tmp.myPtr;
    myPtr->myRefCountInc();
  }


  ideal::ideal(ConstRefRingElem r1, ConstRefRingElem r2)
  {
    vector<RingElem> gens;
    gens.push_back(RingElem(r1));
    gens.push_back(RingElem(r2));
    if (!HasUniqueOwner(gens))
      CoCoA_THROW_ERROR(ERR::MixedRings, "ideal(r1, r2)");
    ideal tmp = owner(r1)->myIdealCtor(gens);
    myPtr = tmp.myPtr;
    myPtr->myRefCountInc();
  }


  ideal::ideal(ConstRefRingElem r1, ConstRefRingElem r2, ConstRefRingElem r3)
  {
    vector<RingElem> gens;
    gens.push_back(RingElem(r1));
    gens.push_back(RingElem(r2));
    gens.push_back(RingElem(r3));
    if (!HasUniqueOwner(gens))
      CoCoA_THROW_ERROR(ERR::MixedRings, "ideal(r1, r2, r3)");
    ideal tmp = owner(r1)->myIdealCtor(gens);
    myPtr = tmp.myPtr;
    myPtr->myRefCountInc();
  }


  ideal::ideal(ConstRefRingElem r1, ConstRefRingElem r2, ConstRefRingElem r3, ConstRefRingElem r4)
  {
    vector<RingElem> gens;
    gens.push_back(RingElem(r1));
    gens.push_back(RingElem(r2));
    gens.push_back(RingElem(r3));
    gens.push_back(RingElem(r4));
    if (!HasUniqueOwner(gens))
      CoCoA_THROW_ERROR(ERR::MixedRings, "ideal(r1, r2, r3, r4)");
    ideal tmp = owner(r1)->myIdealCtor(gens);
    myPtr = tmp.myPtr;
    myPtr->myRefCountInc();
  }


  ideal::ideal(const std::vector<RingElem>& gens)
  {
    if (gens.empty()) CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
    if (!HasUniqueOwner(gens)) CoCoA_THROW_ERROR(ERR::MixedRings, "ideal(gens)");
    ideal tmp = owner(gens[0])->myIdealCtor(gens);
    myPtr = tmp.myPtr;
    myPtr->myRefCountInc();
  }


  ideal::ideal(const ring& R, const std::vector<RingElem>& gens)
  {
    if (!gens.empty())
      if (owner(gens[0]) != R || !HasUniqueOwner(gens))
        CoCoA_THROW_ERROR(ERR::MixedRings, "ideal(R, gens)");
    ideal tmp = R->myIdealCtor(gens);
    myPtr = tmp.myPtr;
    myPtr->myRefCountInc();
  }


  ideal::ideal(const ideal& copy)
  {
    copy->myRefCountInc();
    myPtr = copy.myPtr;
  }


  ideal& ideal::operator=(const ideal& rhs)
  {
    // This impl is valid even if lhs and rhs belong to different rings
    rhs->myRefCountInc();
    myPtr->myRefCountDec();
    myPtr = rhs.myPtr;
    return *this;
  }


  ideal::~ideal()
  {
    myPtr->myRefCountDec();
  }


  IdealBase* MakeUnique(ideal& I)
  {
    if (I->myRefCountIsOne()) return I.myPtr;
    IdealBase* NewPtr(I->myClone());
    I->myRefCountDec();  // after myClone for exc. safety
    I.myPtr = NewPtr;
    I.myPtr->myRefCount = 1;
    return I.myPtr;
  }


  //---------------------------------------------------------------------------
  // Syntactic sugar functions

  RingElem operator%(const MachineInt& n, const ideal& I)
  {
    RingElem ans(RingOf(I), n);
    I->myReduceMod(raw(ans));
    return ans;
  }


  RingElem operator%(const BigInt& N, const ideal& I)
  {
    RingElem ans(RingOf(I), N);
    I->myReduceMod(raw(ans));
    return ans;
  }


  RingElem operator%(const BigRat& Q, const ideal& I)
  {
    RingElem ans(RingOf(I), Q);
    I->myReduceMod(raw(ans));
    return ans;
  }


  RingElem operator%(ConstRefRingElem r, const ideal& I)
  {
    if (owner(r) != RingOf(I))
      CoCoA_THROW_ERROR(ERR::MixedRings, "r%I  -- reduction of RingElem modulo an ideal");
    RingElem ans(r);
    I->myReduceMod(raw(ans));
    return ans;
  }


  // Two separate impls in case ring is not commutative
  ideal operator*(ConstRefRingElem r, const ideal& I)
  {
    if (owner(r) != RingOf(I))
      CoCoA_THROW_ERROR(ERR::MixedRings, "r*I, product of RingElem and ideal");
    if (IsZero(r)) return ideal(r);
    if (IsInvertible(r)) return I;
    vector<RingElem> g = gens(I);
    const long n = len(g);
    for (long i=0; i < n; ++i)
      g[i] = r*g[i];
    return ideal(RingOf(I), g);
  }

  ideal operator*(const ideal& I, ConstRefRingElem r)
  {
    if (owner(r) != RingOf(I))
      CoCoA_THROW_ERROR(ERR::MixedRings, "I*r, product of ideal and RingElem");
    if (IsZero(r)) return ideal(r);
    if (IsInvertible(r)) return I;
    vector<RingElem> g = gens(I);
    const long n = len(g);
    for (long i=0; i < n; ++i)
      g[i] = g[i]*r;
    return ideal(RingOf(I), g);
  }


  ideal operator+(const ideal& I, const ideal& J)
  {
    if (RingOf(I) != RingOf(J))
      CoCoA_THROW_ERROR(ERR::MixedRings, "ideal+ideal");
    ideal ans(I);
    MakeUnique(ans)->myAdd(J);
    return ans;
  }


  ideal& operator+=(ideal& I, const ideal& J)
  {
    if (RingOf(I) != RingOf(J))
      CoCoA_THROW_ERROR(ERR::MixedRings, "ideal+=ideal");
    MakeUnique(I)->myAdd(J);
    return I;
  }


  ideal operator*(const ideal& I, const ideal& J)
  {
    if (RingOf(I) != RingOf(J))
      CoCoA_THROW_ERROR(ERR::MixedRings, "ideal*ideal");
    ideal ans(I);
    MakeUnique(ans)->myMul(J);
    return ans;
  }


  ideal& operator*=(ideal& I, const ideal& J)
  {
    if (RingOf(I) != RingOf(J))
      CoCoA_THROW_ERROR(ERR::MixedRings, "ideal*=ideal");
    MakeUnique(I)->myMul(J);
    return I;
  }


  ideal power(const ideal& I, const MachineInt& n)
  {
    const ring R = RingOf(I);
    if (IsNegative(n)) CoCoA_THROW_ERROR(ERR::ReqNonNegative, "power(I,n)");
    unsigned long N = AsUnsignedLong(n);
    if (N == 0) return ideal(one(R));

    // An iterative implementation of binary powering.
    unsigned long bit = 1; while (bit <= N/2) bit <<= 1;
    ideal ans = I;
    while (bit != 1)
    {
      ans *= ans;
      bit >>= 1;
      if (N&bit) ans *= I;
    }
    return ans;
  }


  ideal power(const ideal& I, const BigInt& N)
  {
    long n;
    if (!IsConvertible(n, N)) CoCoA_THROW_ERROR(ERR::ExpTooBig, "power(I,N)");
    return power(I, n);
  }


  ///// ideal intersection(const ideal& I, const ideal& J)
 ideal intersect(const ideal& I, const ideal& J)
  {
    if (RingOf(I) != RingOf(J))
      CoCoA_THROW_ERROR(ERR::MixedRings, "intersect(ideal,ideal)");
    if (IsZero(I) || IsZero(J))
      return ideal(RingOf(I), vector<RingElem>(0));
    // case IsOne in general can be very expensive
    ideal ans(I);
    MakeUnique(ans)->myIntersect(J);
    return ans;
  }


  ideal colon(const ideal& I, const ideal& J)
  {
    if (RingOf(I) != RingOf(J))
      CoCoA_THROW_ERROR(ERR::MixedRings, "colon(ideal,ideal)");
    ideal ans(I);
    MakeUnique(ans)->myColon(J);
    return ans;
  }


  ///  ideal saturation(const ideal& I, const ideal& J)
  ideal saturate(const ideal& I, const ideal& J)
  {
    if (RingOf(I) != RingOf(J))
      CoCoA_THROW_ERROR(ERR::MixedRings, "saturate(ideal,ideal)");
    ideal ans(I);
    MakeUnique(ans)->mySaturate(J);
    return ans;
  }


  bool IsContained(const ideal& I, const ideal& J)
  {
    const vector<RingElem>& gensI = gens(I);
    const long NumGensI = len(gensI);
    for (long i=0; i < NumGensI; ++i)
      if (!IsElem(gensI[i], J)) return false;
    return true;
  }


  bool operator==(const ideal& I, const ideal& J)
  {
    //??? check first whether the myIdealPtrs are equal???
    return IsContained(I, J) && IsContained(J, I);
  }


  bool IsElem(ConstRefRingElem r, const ideal& I)
  {
    if (owner(r) != RingOf(I)) CoCoA_THROW_ERROR(ERR::MixedRings, "IsElem(r, I)");
    return I->IhaveElem(raw(r));
  }


  bool IsElem(const BigInt& r, const ideal& I)
  { return IsElem(RingElem(RingOf(I),r), I); }

  bool IsElem(const BigRat& r, const ideal& I)
  { return IsElem(RingElem(RingOf(I),r), I); }

  bool IsElem(const MachineInt& r, const ideal& I)
  { return IsElem(RingElem(RingOf(I),r), I); }


  std::ostream& operator<<(std::ostream& out, const ideal& I)
  {
    if (!out) return out;  // short-cut for bad ostreams
    I->myOutputSelf(out);
    return out;
  }


  OpenMathOutput& operator<<(OpenMathOutput& OMOut, const ideal& I)
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("cocoa", "ideal");
    OMOut << RingOf(I);
    I->myOutputSelf_OM(OMOut);
    OMOut->mySendApplyEnd();
    return OMOut;
  }


  //---------------------------------------------------------------------------
  // Functions to do with IdealBase
  // The functions for user assertions must check consistency of the assertions:


  bool IdealBase::IamOne() const
  { return IhaveElem(raw(one(myRing()))); }


  bool IdealBase::IamMaximal() const
  {
    if (IamOne()) CoCoA_THROW_ERROR("Must be a proper ideal", "IamMaximal");
    if (IsUncertain3(IamMaximal3Flag))  myTestIsMaximal();
    CoCoA_ASSERT_ALWAYS(!IsUncertain3(IamMaximal3Flag)); // paranoia2020
    return IsTrue3(IamMaximal3Flag);
  }


  bool IdealBase::IamPrimary() const
  {
    if (IamOne()) CoCoA_THROW_ERROR("Must be a proper ideal", "IamPrimary");
    if (IsUncertain3(IamPrimary3Flag))  myTestIsPrimary();
    CoCoA_ASSERT_ALWAYS(!IsUncertain3(IamPrimary3Flag)); // paranoia2020
    return IsTrue3(IamPrimary3Flag);
  }


  bool IdealBase::IamPrime() const
  {
    if (IamOne()) CoCoA_THROW_ERROR("Must be a proper ideal", "IamPrime");
    if (IsUncertain3(IamPrime3Flag))  myTestIsPrime();
    CoCoA_ASSERT_ALWAYS(!IsUncertain3(IamPrime3Flag)); // paranoia2020
    return IsTrue3(IamPrime3Flag);
  }


  bool IdealBase::IamRadical() const
  {
    if (IamOne()) CoCoA_THROW_ERROR("Must be a proper ideal", "IamRadical");
    if (IsUncertain3(IamRadical3Flag))  myTestIsRadical();
    CoCoA_ASSERT_ALWAYS(!IsUncertain3(IamRadical3Flag)); // paranoia2020
    return IsTrue3(IamRadical3Flag);
  }


  //---  Flags assignments ----------------------------------
  // (maximality implies primality, etc)
  bool IdealBase::myAssignMaximalFlag(bool b) const
  {
    if (!IsUncertain3(IamMaximal3Flag))
    {
      if ( IsTrue3(IamMaximal3Flag) != b )
        CoCoA_THROW_ERROR("Contradictory user assertion", "myAssignMaximalFlag");
      return b;
    }
    if (b)
    { // user asserts  maximal true
      IamMaximal3Flag = true3;
      myAssignPrimeFlag(true);
    }
    else
    { // user asserts  maximal false
      IamMaximal3Flag = false3;
      if (IsTrue3(IsPID3(myRing())))  IamPrime3Flag = false3; // no recursion
    }
    return b;
  }


  bool IdealBase::myAssignPrimaryFlag(bool b) const
  {
    if (!IsUncertain3(IamPrimary3Flag))
    {
      if ( IsTrue3(IamPrimary3Flag) != b )
        CoCoA_THROW_ERROR("Contradictory user assertion", "myAssignPrimaryFlag");
      return b;
    }
    if (b) // user asserts  primary true
      IamPrimary3Flag = true3;
    else   // user asserts  primary false
    {
      IamPrimary3Flag = false3;
      myAssignPrimeFlag(false);
    }
    return b;
  }


  bool IdealBase::myAssignPrimeFlag(bool b) const
  {
    if (!IsUncertain3(IamPrime3Flag))
    {
      if ( IsTrue3(IamPrime3Flag) != b )
        CoCoA_THROW_ERROR("Contradictory user assertion", "myAssignPrimeFlag");
      return b;
    }
    if (b) // user asserts  prime true
    {
      IamPrime3Flag = true3;
      if (IsTrue3(IsPID3(myRing())))  IamMaximal3Flag = true3; // no recursion
      myAssignPrimaryFlag(true);
      myAssignRadicalFlag(true);
    }
    else  // user asserts  prime false
    {
      IamPrime3Flag   = false3;
      myAssignMaximalFlag(false);
    }
    return b;
  }


  bool IdealBase::myAssignRadicalFlag(bool b) const
  {
    if (!IsUncertain3(IamRadical3Flag))
    {
      if ( IsTrue3(IamRadical3Flag) != b )
        CoCoA_THROW_ERROR("Contradictory user assertion", "myAssignRadicalFlag");
      return b;
    }
    if (b) // user asserts  radical true
      IamRadical3Flag = true3;
    else   // user asserts  radical false
    {
      IamRadical3Flag = false3;
      myAssignPrimeFlag(b);
    }
    return b;
  }


  //--- printing ----------------------------------------------
  // Simplistic default definition
  void IdealBase::myOutputSelf(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "ideal(";
    const vector<RingElem>& g = myGens();
    const long NumGens = len(g);
    for (long i=0; i < NumGens; ++i)
    {
      out << g[i];
      if (i != NumGens-1) out << ",  ";
    }
    out << ")";
  }


  void IdealBase::myOutputSelf_OM(OpenMathOutput& OMOut) const
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("cocoa", "ideal");
    const vector<RingElem>& G = myGens();
    const long NumGens = len(G);
    OMOut << NumGens;                // Number of gens, should be an attribute???
    for (long i=0; i < NumGens; ++i) // To be reconsidered ???
      OMOut << G[i];                 // ???
    OMOut->mySendApplyEnd();
  }


  //--- member functions only for SparsePolyRing ---------------------------
  // ideal IdealBase::myElim(const std::vector<RingElem>& /*v*/)
  // { // default implementation
  //   CoCoA_THROW_ERROR(ERR::NYI, "myElim (only for SparsePolyRing)");
  // }

  void IdealBase::myAssignElim(const ideal& /*I*/, const std::vector<RingElem>& /*v*/)
  { // default implementation
    CoCoA_THROW_ERROR(ERR::NYI, "myAssignElim (only for SparsePolyRing)");
  }

  void IdealBase::mySaturate(const ideal& /*v*/)
  { // default implementation
    CoCoA_THROW_ERROR(ERR::NYI, "myElim (only for SparsePolyRing)");
  }

  void IdealBase::myMinimalize()
  { // default implementation
    CoCoA_THROW_ERROR(ERR::NYI, "myMinimalize (only for SparsePolyRing)");
  }

  std::vector<ideal> IdealBase::myPrimaryDecomposition() const
  { // default implementation
    CoCoA_THROW_ERROR(ERR::NYI, "myPrimaryDecomposition (only for SparsePolyRing)");
    return std::vector<ideal>(); // just to keep the compiler quiet
  }



} // end of namespace CoCoA
