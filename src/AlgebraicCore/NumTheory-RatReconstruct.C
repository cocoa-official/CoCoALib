//   Copyright (c)  1999,2009-2011  John Abbott and Anna M. Bigatti

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


#include "CoCoA/NumTheory-RatReconstruct.H"

#include "CoCoA/NumTheory-gcd.H"
#include "CoCoA/NumTheory-ContFrac.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRatOps.H"
#include "CoCoA/assert.H"
#include "CoCoA/error.H"
#include "CoCoA/utils.H"

#include <algorithm>
// using std::min;
using std::max;
// using std::swap;
#include <iostream>
using std::ostream;
#include <vector>
using std::vector;

namespace CoCoA
{

  //////////////////////////////////////////////////////////////////
  // Heuristic fault-tolerant rational reconstruction (continued fraction method)

  RatReconstructByContFrac::RatReconstructByContFrac(MachineInt LogEps):
      myCRT(),
      myLogEps(AsSignedLong(LogEps)),
      myResultIsUpToDate(true),
      myResultIsConvincing(false),
      myResult(0,1),
      myBadFactor(0) // obviously "wrong" initial value
  {
    if (myLogEps < 3) CoCoA_THROW_ERROR(ERR::BadArg, "RatReconstructByContFrac ctor");
    if (myLogEps > 10000) CoCoA_THROW_ERROR(ERR::ArgTooBig, "RatReconstructByContFrac ctor");
  }



  void RatReconstructByContFrac::myAddInfo(const MachineInt& r, const MachineInt& m)
  {
    myCRT.myAddInfo(r,m);
    myResultIsUpToDate = false;
    myResultIsConvincing = false;
  }

  void RatReconstructByContFrac::myAddInfo(const BigInt& R, const BigInt& M)
  {
    myCRT.myAddInfo(R,M);
    myResultIsUpToDate = false;
    myResultIsConvincing = false;
  }


  void RatReconstructByContFrac::myUpdateResult() const
  {
    if (myResultIsUpToDate) return;
    CoCoA_ASSERT(!myResultIsConvincing);
//    VerboseLog VERBOSE("RatRecon");
    myResultIsUpToDate = true;
    const BigInt& X = CombinedResidue(myCRT);
    const BigInt& M = CombinedModulus(myCRT);
//    VERBOSE(50) << "X=" << FloatStr(X) << "  M=" << FloatStr(M) << std::endl;

    const long LogM = FloorLog2(M);
    const long LogExtraFactor = myLogEps + FloorLog2(LogM*LogM);
    // Check for zero (allowing some faulty residues)
    if (2*FloorLog2(gcd(X,M)) > LogExtraFactor+LogM)
    {
//      VERBOSE(50) << "ZERO" << std::endl;
      myResult = 0;
      myResultIsConvincing = true;
      return;
    }

    // Main algm
    const BigRat q(X, M);
    BigInt MaxQuot(1);
    for (ContFracIter it(q); !IsEnded(it); ++it)
    {
      MaxQuot = max(MaxQuot, quot(it));
//      if (quot(it) > MaxQuot) MaxQuot = quot(it);
    }

    // Result is "up-to-date".
    // Do first check whether it is "convincing"; if not, return.
    if (FloorLog2(MaxQuot) < myLogEps)
    {
//      VERBOSE(50) << "FAIL-1:  MaxQuot=" << MaxQuot << "  LogEps=" << myLogEps << std::endl;
      return;
    }

    ContFracApproximant CFA;
    for (ContFracIter it(q); quot(it) != MaxQuot; ++it)
      CFA.myAppendQuot(quot(it));

    const BigRat& RS = CFA.myRational();
    myResult = X - M*RS;
    myBadFactor = gcd(M, den(RS));
    // Second convincingness check:
    if (FloorLog2(num(myResult))+FloorLog2(den(myResult))+2*FloorLog2(myBadFactor)+LogExtraFactor > LogM)
    {
//      VERBOSE(50) << "FAIL-2: LogNum=" << FloorLog2(num(myResult)) << "  LogDen=" << FloorLog2(den(myResult)) << "   2*LogBad=" << 2*FloorLog2(myBadFactor) << "  LogExtra=" << LogExtraFactor << "   LogM=" << LogM << std:: endl;
      return;
    }
//    VERBOSE(50) << "SUCCESS  myResult=" << myResult << std::endl;
    myResultIsConvincing = true;
  }


  const BigRat& ReconstructedRat(const RatReconstructByContFrac& reconstructor)
  {
    if (!IsConvincing(reconstructor)) // automatically updates
      CoCoA_THROW_ERROR("Result is not convincing","ReconstructedRat(RatReconstructByContFrac)");
    return reconstructor.myResult;
  }


  bool IsConvincing(const RatReconstructByContFrac& reconstructor)
  {
    reconstructor.myUpdateResult();
    return reconstructor.myResultIsConvincing;
  }


  const BigInt& BadMFactor(const RatReconstructByContFrac& reconstructor)
  {
    if (!IsConvincing(reconstructor)) // automatically updates
      CoCoA_THROW_ERROR("Result is not convincing","BadMFactor(RatReconstructByContFrac)");
    return reconstructor.myBadFactor;
  }


  std::ostream& operator<<(std::ostream& out, const RatReconstructByContFrac& reconstructor)
  {
    if (!out) return out;  // short-cut for bad ostreams
    out << "RatReconstructByContFrac(CRT=" << reconstructor.myCRT
//        << ", threshold=" << reconstructor.myThresholdValue
        << ", LogEps=" << reconstructor.myLogEps
        << ", IsConvincing=" << reconstructor.myResultIsConvincing;
    if (reconstructor.myResultIsConvincing)
      out << ", result=" << reconstructor.myResult;
    out << ")";
    return out;
  }


  //////////////////////////////////////////////////////////////////
  // Heuristic fault-tolerant rational reconstruction (2x2 lattice method)


  const long RatReconstructByLattice::ourDefaultSafetyFactor = 4096;

  BigInt RatReconstructByLattice::myCheckSafetyFactor(const BigInt& SafetyFactor)
  {
    // SafetyFactor == 0 --> use default value (see static data mem ourDefaultSafetyFactor)
    if (SafetyFactor < 0) CoCoA_THROW_ERROR(ERR::ReqNonNegative, "RatReconstructByLattice ctor");
    if (SafetyFactor > 0) return SafetyFactor;
    return BigInt(ourDefaultSafetyFactor);
  }

  RatReconstructByLattice::RatReconstructByLattice(const MachineInt& SafetyFactor):
      myCRT(),
      mySafetyFactor(myCheckSafetyFactor(BigInt(SafetyFactor))),
      myResultIsUpToDate(true),
      myResultIsConvincing(true),
      myResult(0,1),
      myBadFactor(0) // obviously "wrong" initial value
  {}

  RatReconstructByLattice::RatReconstructByLattice(const BigInt& SafetyFactor):
      myCRT(),
      mySafetyFactor(myCheckSafetyFactor(SafetyFactor)),
      myResultIsUpToDate(true),
      myResultIsConvincing(true),
      myResult(0,1),
      myBadFactor(0) // obviously "wrong" initial value
  {}


  void RatReconstructByLattice::myAddInfo(const MachineInt& r, const MachineInt& m)
  {
    myCRT.myAddInfo(r,m);
    myResultIsUpToDate = false;
    myResultIsConvincing = false;
  }

  void RatReconstructByLattice::myAddInfo(const BigInt& R, const BigInt& M)
  {
    myCRT.myAddInfo(R,M);
    myResultIsUpToDate = false;
    myResultIsConvincing = false;
  }


  void RatReconstructByLattice::myUpdateResult() const
  {
    if (myResultIsUpToDate) return;
    myResultIsUpToDate = true;
    const BigInt R = CombinedResidue(myCRT);
    const BigInt M = CombinedModulus(myCRT);

    BigInt a0 = M;
    BigInt b0;
    BigInt a1 = R;
    BigInt b1(1);

    do
    {
      const BigInt q = round(BigRat(a0*a1+b0*b1, a1*a1+b1*b1));
      a0 -= q*a1;
      b0 -= q*b1;
      swap(a0,a1);
      swap(b0,b1);
    } while (a1*a1+b1*b1 < a0*a0+b0*b0);
    if (mySafetyFactor*(a0*a0+b0*b0) >= M) return; // result is "up-to-date" but NOT "convincing"
    myResult = BigRat(a0,b0);
    myBadFactor = gcd(gcd(a0,b0),M);
    myResultIsConvincing = true;
  }


  const BigRat& ReconstructedRat(const RatReconstructByLattice& reconstructor)
  {
    if (!IsConvincing(reconstructor)) // automatically updates
      CoCoA_THROW_ERROR("Result is not convincing","ReconstructedRat(RatReconstructByLattice)");
    return reconstructor.myResult;
  }

  bool IsConvincing(const RatReconstructByLattice& reconstructor)
  {
    reconstructor.myUpdateResult();
    return reconstructor.myResultIsConvincing;
  }


  const BigInt& BadMFactor(const RatReconstructByLattice& reconstructor)
  {
    if (!IsConvincing(reconstructor)) // automatically updates
      CoCoA_THROW_ERROR("Result is not convincing","BadMFactor(RatReconstructByLattice)");
    return reconstructor.myBadFactor;
  }

  std::ostream& operator<<(std::ostream& out, const RatReconstructByLattice& reconstructor)
  {
    if (!out) return out;  // short-cut for bad ostreams
    out << "RatReconstructByLattice(CRT=" << reconstructor.myCRT
        << ", SafetyFactor=" << reconstructor.mySafetyFactor
        << ", IsConvincing=" << reconstructor.myResultIsConvincing;
    if (reconstructor.myResultIsConvincing)
      out << ", result=" << reconstructor.myResult;
    out << ")";
    return out;
  }

  //////////////////////////////////////////////////////////////////
  // FTRR

  BigInt ComputeMmax(long e, vector<long> mod)
  {
    sort(mod.begin(), mod.end());
    const long s = len(mod);
    BigInt ans(1);
    for (long i=s-e; i < s; ++i)
      ans *= mod[i];
    return ans;
  }

  // BUG BUG BUG  STOPGAP impl; <---> generic impl in tmp.H
  BigInt product(const vector<long>& mod)
  {
    BigInt ans(1);
    const long s = len(mod);
    for (long i=0; i < s; ++i)
      ans *= mod[i];
    return ans;
  }

  BigRat RatReconstructWithBounds(long e, const BigInt& P, const BigInt& Q, const std::vector<long>& res, const std::vector<long>& mod)
  {
    const BigRat FAILURE(BigRat::OneOverZero); // "impossible rational" used to indicate failure
    const long s = len(res);
    if (len(mod) != s) CoCoA_THROW_ERROR(ERR::BadArg, "FTRR");
    if (e < 0 || 2*e >= s || P < 1 || Q < 1) CoCoA_THROW_ERROR(ERR::BadArg, "FTRR");
    const BigInt Mmax = ComputeMmax(e, mod);
    if (2*P*Q*power(Mmax,2) >= product(mod)) CoCoA_THROW_ERROR(ERR::BadArg, "FTRR");

    {
      long CountZeroes = 0;
      for (long i=0; i < s; ++i)
        if (res[i]%mod[i] == 0) ++CountZeroes;
      if (CountZeroes >= s-e) return BigRat(0);
    }

    CRTMill CRT;
    for (long i=0; i < s; ++i)
      CRT.myAddInfo(res[i], mod[i]);
    const BigInt X = CombinedResidue(CRT);
    const BigInt M = CombinedModulus(CRT);

    if (gcd(X,M) > P*Mmax) return FAILURE;

    BigInt u1(1); BigInt u2;    BigInt u3(M);
    BigInt v1;    BigInt v2(1); BigInt v3(X);
    while (abs(v2) <= Q*Mmax)
    {
      const BigInt q = u3/v3; // floor division!
      u1 = u1 - q*v1;  swap(u1, v1);
      u2 = u2 - q*v2;  swap(u2, v2);
      u3 = u3 - q*v3;  swap(u3, v3);
    }
    const BigRat r = X + M*BigRat(u1, u2);
    if (abs(num(r)) > P || den(r) > Q) return FAILURE;
    long CountBadModuli=0;
    for (long i=0; i < s; ++i)
      if (gcd(u2, mod[i]) > 1) ++CountBadModuli;
    if (CountBadModuli > e) return FAILURE;
    return r;
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/NumTheory-RatReconstruct.C,v 1.8 2024/07/02 15:32:23 bigatti Exp $
// $Log: NumTheory-RatReconstruct.C,v $
// Revision 1.8  2024/07/02 15:32:23  bigatti
// Summary: Changed error codes: LogZero into ReqNonZero
//
// and Not... into ReqNonNegative, ReqNonNegativeGrading, ReqPositive, ReqPositiveGrading
//
// Revision 1.7  2022/02/18 14:11:55  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.6  2021/01/15 16:59:33  abbott
// Summary: Redmine 1563: new ctor for BigRat direct from integer
//
// Revision 1.5  2021/01/07 15:16:51  abbott
// Summary: Corrected copyright
//
// Revision 1.4  2020/06/17 15:49:24  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.3  2020/01/26 14:41:59  abbott
// Summary: Revised includes after splitting NumTheory (redmine 1161)
//
// Revision 1.2  2019/03/18 11:57:39  abbott
// Summary: Added missing includes
//
// Revision 1.1  2019/03/18 11:24:19  abbott
// Summary: Split NumTheory into several smaller files
//
//
//
