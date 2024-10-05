//   Copyright (c)  2022  John Abbott and Anna M. Bigatti

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


#include "CoCoA/NumTheory-root.H"
#include "CoCoA/NumTheory-gcd.H"
#include "CoCoA/NumTheory-factor.H"
#include "CoCoA/NumTheory-modular.H"

#include "CoCoA/BigIntOps.H"
// #include "CoCoA/assert.H"
#include "CoCoA/bool3.H"
// #include "CoCoA/config.H"
// #include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/interrupt.H"
#include "CoCoA/ProgressReporter.H"
// #include "CoCoA/utils.H"

// #include <algorithm>
// using std::min;
// // using std::max;
// // using std::swap;
// #include <cmath>
// // myThreshold uses floor,pow,sqrt
// // #include <cstdlib>
// // using std::ldiv;
// #include <limits>
// using std::numeric_limits;
// #include <vector>
// using std::vector;

namespace CoCoA
{

  namespace // anonymous
  {

    // assumes k prime
    // Quick Check if N is k-th power
    // RTN: false3 means "definitely no", o/w uncertain3 means "uncertain"
    bool3 IsPowerQuick(const BigInt& N, long k)
    {
      if (IsOdd(k)) k *= 2; // step size between number to check for primality, always even
      long p = 1+(1000000000/k)*k; // 1000000000 is "magic number"
      int counter = 0;
      while (counter < 5)
      {
        do { p += k; } while (!IsPrime(p));
        ++counter;
        const long r = (k==2)?(p-1)/k:(p-1)/(k/2); // exact division!
        if (PowerMod(N, r, p) != 1)
          return false3;
      }
      return uncertain3;
    }


    // Return smallest r such that there is k with N == r^k
    // if UPBexp != 0 then restrict to k <= UPBexp
    BigInt StarRoot_naive(BigInt N, long UPBexp)
    {
      if (N < 0) CoCoA_THROW_ERROR("arg1: N must be non-negative", "StarRoot_naive");
      if (N < 2) return N;  // or ERROR???

      long logN = FloorLog2(N);
      if (UPBexp == 0) UPBexp = logN;
      double lnN = log(N);
      long MaxK = 1+logN;
      long multiplicity = 0;
      PrimeSeq pseq;
///      ProgressReporter report(10);
      for (PrimeSeq kseq; *kseq <= MaxK  && *kseq <= UPBexp; ++kseq)
      {
        long k = *kseq;
///        report(k);
        CheckForInterrupt("StarRoot_naive");

        for (int i=0; i < 1; ++i)
        {
          long p = *pseq;
          ++pseq;
          if(N%p!=0) continue;
          const long m = FactorMultiplicity(p,N);
          if (m == 0) continue;
          multiplicity = gcd(m, multiplicity);
          if (multiplicity == 1) return N;
        }

        MaxK  = static_cast<long>(std::ceil(lnN/std::log(static_cast<double>(*pseq)))); // floor(0.01+lnN/log(*pseq)) ???
        if (k > MaxK) break;
        if (gcd(k,multiplicity) == 1) continue;
        if (IsFalse3(IsPowerQuick(N,k))) continue;
        // Quick check permits N to be k-th power, so compute k-th root and check
        BigInt r = FloorRoot(N,k);
        if (N != power(r,k)) continue;
        // N was perfect k-th power; now check if k^2-th power or k^3-th power etc.
        do
        {
          N = r;
          logN /= k;
          lnN /= k;
          multiplicity /= k;
          if (IsFalse3(IsPowerQuick(N,k))) break;
          r = FloorRoot(N,k);
        }
        while (N == power(r,k));
      }
      if (multiplicity > 1) N = FloorRoot(N, multiplicity);
      return N;
    }

  } // end of namespace anonymous

  //------------------------------------------------------------------

  // NOTE: should also check that multiplicity is a mult of k

  // Certify not a k-th power via Fermat Little Thm
  // Returns (smallest) prime p st N is not k-th power mod p.
  // Grunwald-Wang thm says this always terminates.
  long CertifyNotPower(const BigInt& N, long k)
  {
    const char* const FnName = "CertifyNotPower";
    if (N < 2) CoCoA_THROW_ERROR("arg 1: N must be >= 2", FnName);
    if (k < 2 || !IsPrime(k)) CoCoA_THROW_ERROR("arg 2: k must be prime", FnName);

    const long step = (k==2)?k:2*k;
    long p = 1;
    while (true)
    {
      p += step;
      if (!IsPrime(p)) continue;
      if (IsDivisible(N,p)) continue; // could be clever: check multiplicity (& if multiplicity is mult of k then check if what's left is a k-th power)
      CheckForInterrupt("CertifyNotPower");
      const long r = (p-1)/k;
      if (PowerMod(N, r, p) != 1)
        return p;
    }
  }


  // Similar to above, but increase prime quickly until a certain threshold
  // (works better if input is a factorial or primorial).
  long CertifyNotPower_(const BigInt& N, long k)
  {
    const char* const FnName = "CertifyNotPower";
    if (N < 2) CoCoA_THROW_ERROR("arg 1: N must be >= 2", FnName);
    if (k < 2 || !IsPrime(k)) CoCoA_THROW_ERROR("arg 2: k must be prime", FnName);
    if (IsOdd(k)) k *= 2;
    constexpr long thresh = 10000000;
    long p = 1;
    while (true)
    {
      p += k;
      if (!IsPrime(p)) continue;
      if (IsDivisible(N,p)) continue; // be clever?  (see above)
      const long r = (p-1)/k;
      if (PowerMod(N, r, p) != 1)
        return p;
      if (p < thresh) p = 2*p-1;
    }
  }




  BigInt StarRoot(BigInt N, long UPBexp)
  {
    return StarRoot_naive(N, UPBexp);
  }

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/NumTheory-root.C,v 1.6 2023/02/16 20:47:39 abbott Exp $
// $Log: NumTheory-root.C,v $
// Revision 1.6  2023/02/16 20:47:39  abbott
// Summary: Cleaned up impl of StarRoot
//
// Revision 1.5  2023/02/14 19:59:20  abbott
// Summary: Partially improved & cleaned (stil needs work!)
//
// Revision 1.4  2022/08/05 11:37:29  abbott
// Summary: Added comments
//
// Revision 1.3  2022/02/18 14:11:56  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.2  2022/02/02 09:24:46  abbott
// Summary: Minor improvements (redmine 1657)
//
// Revision 1.1  2022/01/20 19:15:57  abbott
// Summary: Added new fns StarRoot and CertifyNotPower
//
//
//
