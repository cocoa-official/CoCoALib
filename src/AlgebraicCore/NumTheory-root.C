//   Copyright (c)  2022,2026  John Abbott and Anna M. Bigatti

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
#include "CoCoA/verbose.H"
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
      if (IsOdd(k))  k *= 2; // step size between number to check for primality, always even
      long p = 1+(1000000000/k)*k; // 1000000000 is "magic number"
      int counter = 0;
      int CountUseful = 0;
      while (counter < 5)
      {
        do { p += k; } while (!IsPrime(p));
        // At this point p is a prime which is 1 modulo k
        ++counter;
        const long r = N%p;
        if (r == 0 || r == 1)  continue;
        ++CountUseful;
        const long e = (k==2)?(p-1)/k:(p-1)/(k/2); // exact division!
        if (PowerMod(r, e, p) != 1)
          return false3;
      }
      // if (CountUseful < 2)  DO_A_BIGGER_PRIME (or even 2)
      return uncertain3;
    }


/****************************************************************************/
// Delete StarRoot_BUGGY !!!
/****************************************************************************/

    // Return smallest r such that there is k with N == r^k
    // if UPBexp != 0 then restrict to k <= UPBexp
    BigInt StarRoot_BUGGY(BigInt N, long UPBexp)
    {
      if (N < 0)  CoCoA_THROW_ERROR1(ERR::ReqNonNegative);
      if (N < 2)  return N;  // or ERROR???

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

    // Return smallest r such that there is s with N == r^s
    // if UPBexp != 0 then restrict to k <= UPBexp
    //
    // This function needs some explanation: it's not really
    // advanced, but becomes involved when sorting out the details.
    // The idea is to use two distinct approaches "in parallel":
    // (1) if p|N then s | FactorMultiplicity(p,N)
    // (2) if FloorRoot(N,k) is exact then k|s -- in fact,
    //     we replace N by FloorRoot(N,k) then continue.
    //
    // The invariant for the main loop is AND of:
    // (*) N is not divisible by any prime < *pseq;
    // (*) N is not j-th power for any j < *kseq,
    //     equiv s is not divisible by any j < *kseq
    // (*) s|multiplicity
    BigInt StarRoot_improved(BigInt N, long UPBexp)
    {
      if (N < 0)  CoCoA_THROW_ERROR1(ERR::ReqNonNegative);
      if (N < 2)  return N;

      VerboseLog VERBOSE("StarRoot_improved");
      long logN = FloorLog2(N);
      if (UPBexp == 0)  UPBexp = logN;
      double lnN = log(N)+1.0/65536.0; // better to mult by (1+epsilon)  ???
      long MaxK = 1+logN;
      long multiplicity = 0; // if N = r^s then s divides multiplicity
      double lnKnownFactor = 0.0;  // may help exit loop sooner
      PrimeSeq pseq;  // primes for trial division (see inner loop below)
      long p = 1; // init must be 1; "follows" pseq in inner loop below
      PrimeSeq kseq;
      while (*kseq <= MaxK  && *kseq <= UPBexp)
      {
        CheckForInterrupt("StarRoot_improved");
        /*const*/ long k = *kseq;
        VERBOSE(200) << "OUTER LOOP (prime factor of base) p=" << *pseq << std::endl;
        VERBOSE(200) << "OUTER LOOP (prime factor of exponent) k=" << k << std::endl;
        VERBOSE(200) << "OUTER LOOP: lnKnownFactor+k*ln(p)="<<lnKnownFactor+k*std::log(p)<<"   lnN="<<lnN<<std::endl;
        // Try some prime divisors: if we find a divisor, we can update multiplicity:
        for (int i=0; i < 5; ++i, ++pseq)
        {
          p = *pseq;
          /*const*/ long m = FactorMultiplicity(p,N);
          if (m == 0)  continue;
          VERBOSE(200) << "Prime factor p = " << p << std::endl;
          VERBOSE(200) << "Multiplicity of p: " << m << std::endl;
          if (multiplicity == 0)  for (PrimeSeq qseq; *qseq < k; ++qseq)  { m = CoprimeFactor(m,*qseq); }
          multiplicity = gcd(m, multiplicity);
          VERBOSE(200) << "Overall multiplicity must divide " << multiplicity << std::endl;
          if (multiplicity == 1)  return N;
          lnKnownFactor += m*std::log(p);
          ++pseq;
          break;
        }
        VERBOSE(200) << "PROGRESS INFO: p = " << p << "  and  *pseq = " << *pseq << std::endl;
        if (lnN-lnKnownFactor < 0.6) // note that ln(2) > 0.6
        {
          // We get here only if we have completely factorized N
          CoCoA_ASSERT(!IsFalse3(IsPowerQuick(N,multiplicity)));
          VERBOSE(200) << "Factorization complete: multiplicity = " << multiplicity << std::endl;
          /*const*/ BigInt r = FloorRoot(N,multiplicity);
          CoCoA_ASSERT(N == power(r,multiplicity));
          return r;
        }
        // If N is a power then the exponent/multiplicity must be at least k,
        // also N must have a prime factor larger than p; so if (p*KnownFactor)^k > N
        // then we know that N cannot be a power.
        if (lnKnownFactor+k*std::log(p) > lnN)  return N;
        // Advance k (by at most 100 "steps") until we find a plausible next value for k
        bool plausible = false;
        for (int i=0; i < 100; ++i)
        {
          if (multiplicity%k == 0)
          {
            if (!IsFalse3(IsPowerQuick(N,k)))   { plausible = true; break; }
            if (multiplicity > 0)
            {
              multiplicity = CoprimeFactor(multiplicity,k);  MaxK = std::min(MaxK, multiplicity); 
              VERBOSE(200) << "Reduced multiplicity to " << multiplicity << std::endl;
            }
          }
          ++kseq; k = *kseq; // could be k = *++kseq;
          if (multiplicity != 0 && multiplicity/k < k)  { if (IsPrime(multiplicity)) { k = multiplicity; } else { multiplicity = 1; MaxK = 1; } }
          if (k > MaxK)  break;
        }
        if (!plausible)  { ++kseq; continue; }
        VERBOSE(200) << "Handling factor " << k << " in multiplicity" << std::endl;
        if (multiplicity != 0)
          multiplicity = CoprimeFactor(multiplicity, k);
        VERBOSE(200) << "New multiplicity is " << multiplicity << std::endl;
        if (!IsFalse3(IsPowerQuick(N,k)))
        {
          VERBOSE(200) << "Costly check" << std::endl;
          // Quick check permits N to be k-th power, so compute k-th root and check
          BigInt r = FloorRoot(N,k);
          if (N == power(r,k))  //  we expect this to be true with high probability
          {
            // N was perfect k-th power; now check if k^2-th power or k^3-th power etc.
            do
            {
              N = r;
              logN /= k;
              lnN /= k;
///          multiplicity /= k;  // exact division
              lnKnownFactor /= k; // still natural log of integer factor of N
              if (IsFalse3(IsPowerQuick(N,k)))  break;
              r = FloorRoot(N,k);
            }
            while (N == power(r,k));
            MaxK = static_cast<long>(std::ceil(lnN/std::log(static_cast<double>(*pseq))));
            VERBOSE(200) << "Update: MaxK = " << MaxK << std::endl;
          }
        }
        if (multiplicity == 1)  return N;
        ++kseq;
      }
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
    if (N < 2)  CoCoA_THROW_ERROR2(ERR::BadArg, "1st arg must be >= 2");
    if (k < 2 || !IsPrime(k))  CoCoA_THROW_ERROR2(ERR::BadArg, "arg 2 must be prime");

    const long step = (k==2)?k:2*k;
    long p = 1;
    while (true)
    {
      p += step;
      if (!IsPrime(p))  continue;
      if (IsDivisible(N,p))  continue; // could be clever: check multiplicity (& if multiplicity is mult of k then check if what's left is a k-th power)
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
    if (N < 2)  CoCoA_THROW_ERROR2(ERR::BadArg, "1st arg must be >= 2");
    if (k < 2 || !IsPrime(k))  CoCoA_THROW_ERROR2(ERR::BadArg, "2nd arg must be prime");
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




  BigInt StarRoot(BigInt N, long UPBexp/*=0*/)
  {
    return StarRoot_improved(N, UPBexp);
  }

  long StarRoot(long n)
  {
    return ConvertTo<long>(StarRoot(BigInt(n)));
  }

  // JUST FOR DEBUGGING
  BigInt StarRoot_OLD(BigInt N, long UPBexp/*=0*/)
  {
    return StarRoot_BUGGY(N, UPBexp);
  }


} // end of namespace CoCoA
