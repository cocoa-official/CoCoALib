//   Copyright (c)  1999,2009-2011,2018  John Abbott and Anna M. Bigatti

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


#include "CoCoA/NumTheory-prime.H"
#include "CoCoA/NumTheory-factor.H"
#include "CoCoA/NumTheory-modular.H"
#include "CoCoA/assert.H"
#include "CoCoA/error.H"
#include "CoCoA/interrupt.H"
#include "CoCoA/random.H"
/////?????#include "CoCoA/SignalWatcher.H"
#ifdef CoCoA_DEBUG
#include "CoCoA/NumTheory-gcd.H"
#endif


#include <iostream>
using std::ostream;
#include <limits>
using std::numeric_limits;
#include <vector>
using std::vector;

namespace CoCoA
{

  // Anonymous namespace for file local "static" variables.
  namespace
  {
    const int ProbPrimeIters = 25; // default num iters used by the single arg "ProbPrime" fns

    // Two tables used by NextPrime, NextProbPrime, PrevPrime, and PrevProbPrime
    // n+skip[n%30] is the next number after n not divisible by 2, 3 or 5
    // n-fall[n%30] is the largest number below n not divisible by 2, 3 or 5
    const int skip[30] = {1,6,5,4,3,2,1,4,3,2,1,2,1,4,3,2,1,2,1,4,3,2,1,6,5,4,3,2,1,2};
    const int fall[30] = {1,2,1,2,3,4,5,6,1,2,3,4,1,2,1,2,3,4,1,2,1,2,3,4,1,2,3,4,5,6};


    // Miller-Rabin strong pseudo-prime test (for machine integers).
    // Tests whether n is a strong pseudo-prime to base b.
    // Assumes 1 < b < n, and n odd & n*n < MaxULong.
    bool StrongPseudoPrime(unsigned long b, unsigned long n)
    {
      CoCoA_ASSERT(n <= MaxSquarableInteger<unsigned long>());
      CoCoA_ASSERT((n&1) == 1);
      CoCoA_ASSERT(1 < b && b < n);

      // Compute q,r such that  n-1 = q*2^r  with q odd.
      unsigned long q = n-1;
      int r = 0;
      while ((q&1) == 0)
      {
        ++r;
        q /= 2;
      }

      // Main loop
      const unsigned long n1 = n-1;
      unsigned long power = PowerMod(b, q, n);
      if (power == n1 || power == 1) return true;
      for (int i=0; i < r; i++)
      {
        power = (power*power)%n;
        if (power == n1) return true;
        if (power == 1) return false;
      }
      return false;
    }

    // SAME AS StrongPseudoPrime but uses BigInt for PowerMod;
    // necessary when n > MaxSquarableInteger!
    bool StrongPseudoPrime_SLOW(unsigned long b, unsigned long n)
    {
//      CoCoA_ASSERT(n > MaxSquarableInteger<unsigned long>());
      CoCoA_ASSERT((n&1) == 1);
      CoCoA_ASSERT(1 < b && b < n);

      // Compute q,r such that  n-1 = q*2^r  with q odd.
      unsigned long q = n-1;
      int r = 0;
      while ((q&1) == 0)
      {
        ++r;
        q /= 2;
      }

      // Main loop
      const BigInt N(n);
      const unsigned long n1 = n-1;
      unsigned long power = ConvertTo<unsigned long>(PowerMod(b, q, N));
      if (power == n1 || power == 1) return true;
      for (int i=0; i < r; i++)
      {
        power = ConvertTo<unsigned long>(PowerMod(power,2,N));
        if (power == n1) return true;
        if (power == 1) return false;
      }
      return false;
    }


    // Quick (but ugly) -- creams a std::bitset based version
    // Used only by IsSmallPrime.  Compare with alternative impl below.
    // WARNING:  See also fn StrongPseudoPrime64
    bool HasSmallFactor37(unsigned long n)
    {
      if ((n&1) == 0) return true;

      switch (n%33)
      {
      case 0: case 3: case 6: case 9: case 11: case 12:
      case 15: case 18: case 21: case 22: case 24: case 27: case 30:
        return true;
      }
      switch (n%35)
      {
      case 0: case 5: case 7: case 10: case 14: case 15:
      case 20: case 21: case 25: case 28: case 30:
        return true;
      }
      
      return (n%13 == 0) || (n%17 == 0) || (n%19 == 0) || (n%23 == 0) || (n%29 == 0) || (n%31 == 0) || (n%37 == 0);
    }

    // // This seems to be marginally faster than version above on my machine.  Which is cleaner???
    // // Compare with alternative impl above.
    // bool HasSmallFactor37(unsigned long n)
    // {
    //   if ((n&1) == 0) return true;

    //   const unsigned char n33 = n%33;
    //   static const unsigned char tbl33[33] = {1,0,0,1,0,0,1,0,0,1,0,1,1,0,0,1,0,0,1,0,0,1,1,0,1,0,0,1,0,0,1,0,0};
    //   if (tbl33[n33]) return true;

    //   const unsigned char n35 = n%35;
    //   static const unsigned char tbl35[35] = {1,0,0,0,0,1,0,1,0,0,1,0,0,0,1,1,0,0,0,0,1,1,0,0,0,1,0,0,1,0,1,0,0,0,0};
    //   if (tbl35[n35]) return true;
      
    //   return (n%13 == 0) || (n%17 == 0) || (n%23 == 0) || (n%29 == 0) || (n%31 == 0) || (n%37 == 0);
    // }


    // Tests whether n is prime (for 0 <= n <= 36)
    bool IsSmallPrime37(unsigned long n)
    {
      CoCoA_ASSERT(n <= 37);
      switch (n)
      {
      case 2: case 3: case 5: case 7:
      case 11: case 13: case 17: case 19:
      case 23: case 29:
      case 31: case 37:
      // case 41: case 43: case 47:
      // case 53: case 59:
      // case 61:
        return true;
      default:
        return false;
      }
    }


    // (1) "fast" machine integer implementation; see below for big integer implementation.
    // Use this defn if n^2 fits into an unsigned long.
    bool IsPrime_squarable(unsigned long n)
    {
      // Caller has already checked for small prime facs
      CoCoA_ASSERT(n <= MaxSquarableInteger<unsigned long>());
      CoCoA_ASSERT(n >= 1681); // 1681 = NextPrime(37)^2 where 37 is max small prime filtered out.
      if (n <= 4294967295UL) // ALWAYS TRUE if unsigned long has 32 bits
      {
        // Case n < 2^32.
        // exclude 4 exceptions (up to 2^32)
        if (n ==  746331041UL ||
            n == 2840871041UL ||
            n == 3014101261UL ||
            n == 3215031751UL)
          return false;

        // Enough to check bases 2, 5, 7; order 7,2,5 seems slightly better.
        return StrongPseudoPrime(7, n) &&
          StrongPseudoPrime(2, n) &&
          StrongPseudoPrime(5, n);
      }

// Use CPP to eliminate next block on 32-bitters to avoid warnings about "constant too big"
// #if ULONG_MAX > 4294967295UL
#if (ULONG_MAX/1000)/1000 > 4294
      // Get here only if unsigned long has more than 32 bits.
      // Cannot use "n < (1UL << 40)" as it is "malformed" on 32-bit platform.
      if ((n >> 20) < (1UL << 20)) // really just testing if n < 2^40
      {
        // Case n < 2^40
        // exclude  9 exceptions
        if (n == 321503175UL ||
            n == 118670087467UL ||
            n == 307768373641UL ||
            n == 315962312077UL ||
            n == 354864744877UL ||
            n == 457453568161UL ||
            n == 528929554561UL ||
            n == 546348519181UL ||
            n == 602248359169UL) return false;

        return StrongPseudoPrime(7, n) &&
          StrongPseudoPrime(2, n) &&
          StrongPseudoPrime(5, n) &&
          StrongPseudoPrime(3, n);

      }
#endif
      
      // NB we should get here only if unsigned long has more than 80 bits.
      // Case n < 2^64; see OEIS A014233.
      // Check for being strong pseudo-prime to bases 2,3,5,...,37
      CoCoA_ASSERT((n >> 32) < (1UL << 32)-1); // cannot portably shift by more than 32
      return StrongPseudoPrime(37, n) &&
        StrongPseudoPrime(2, n) &&
        StrongPseudoPrime(3, n) &&
        StrongPseudoPrime(5, n) &&
        StrongPseudoPrime(7, n) &&
        StrongPseudoPrime(11, n) &&
        StrongPseudoPrime(13, n) &&
        StrongPseudoPrime(17, n) &&
        StrongPseudoPrime(19, n) &&
        StrongPseudoPrime(23, n) &&
        StrongPseudoPrime(29, n) &&
        StrongPseudoPrime(31, n);
    }


    // (2) "slow" big integer implementation; see above for machine integer implementation.
    // Use this defn if n^2 DOES NOT FIT into an unsigned long.
    bool IsPrime_NotSquarable(unsigned long n)
    {
      // Caller has already checked for small prime facs
      CoCoA_ASSERT(n > MaxSquarableInteger<unsigned long>());
      // Since StrongPseudoPrime_SLOW is slow, we filter out a few more small primes...
      if (n%41 == 0 || n%43 == 0 || n%47 == 0 || n%53 == 0 || n%59 == 0 || n%61 == 0) return false;

      if (n < 4294967295UL) // ALWAYS TRUE if unsigned long has 32 bits
      {
        // Case n < 2^32.
        // exclude 4 exceptions (up to 2^32)
        if (n ==  746331041UL ||
            n == 2840871041UL ||
            n == 3014101261UL ||
            n == 3215031751UL)
          return false;

        return StrongPseudoPrime_SLOW(2, n) &&
          StrongPseudoPrime_SLOW(5, n) &&
          StrongPseudoPrime_SLOW(7, n);
      }

// Use CPP to eliminate next block on 32-bitters to avoid warnings about "constant too big"
// #if ULONG_MAX > 4294967295UL
#if (ULONG_MAX/1000)/1000 > 4294
      // Get here only if unsigned long has more than 32 bits.
      // Cannot use "n < (1UL << 40)" as it is "malformed" on 32-bit platform.
      if ((n >> 20) < (1UL << 20)) // really just testing if n < 2^40
      {
        // Case n < 2^40
        // exclude  9 exceptions
        if (n == 321503175UL ||
            n == 118670087467UL ||
            n == 307768373641UL ||
            n == 315962312077UL ||
            n == 354864744877UL ||
            n == 457453568161UL ||
            n == 528929554561UL ||
            n == 546348519181UL ||
            n == 602248359169UL) return false;

        return StrongPseudoPrime_SLOW(2, n) &&
          StrongPseudoPrime_SLOW(3, n) &&
          StrongPseudoPrime_SLOW(5, n) &&
          StrongPseudoPrime_SLOW(7, n);

      }
#endif
      
      // Case n < 2^64; see OEIS A014233.
      // Check for being strong pseudo-prime to bases 2,3,5,...,37
      CoCoA_ASSERT((n >> 32) < (1UL << 32)-1); // cannot portably shift by more than 32
      return StrongPseudoPrime_SLOW(37, n) &&
        StrongPseudoPrime_SLOW(2, n) &&
        StrongPseudoPrime_SLOW(3, n) &&
        StrongPseudoPrime_SLOW(5, n) &&
        StrongPseudoPrime_SLOW(7, n) &&
        StrongPseudoPrime_SLOW(11, n) &&
        StrongPseudoPrime_SLOW(13, n) &&
        StrongPseudoPrime_SLOW(17, n) &&
        StrongPseudoPrime_SLOW(19, n) &&
        StrongPseudoPrime_SLOW(23, n) &&
        StrongPseudoPrime_SLOW(29, n) &&
        StrongPseudoPrime_SLOW(31, n);
    }


    // We could create some intermediate fns: again see OEIS A014233 for various limits
    // Or we could add more "exceptions" with fewer StrongPseudoPrime tests.
    // Here are the composite (2,5,7)-SPSP up to 2^35 (about 34*10^9):
    // 7535192941, 12337298821, 27126501751, 29875480447, 30926647201, 31759334209

    // Here are the first few composite (2,3,5,7)-SPSP:
    // 3215031751  myFactors := [151,  751,  28351]
    // 118670087467 myFactors=[172243,  688969]
    // 307768373641 myFactors=[392281,  784561]
    // 315962312077 myFactors=[281053,  1124209]
    // 354864744877 myFactors=[297853,  1191409]
    // 457453568161 myFactors=[390493,  1171477]
    // 528929554561 myFactors=[419893,  1259677]
    // 546348519181 myFactors=[522661,  1045321]
    // 602248359169 myFactors=[347059,  1735291]
    // Checked up to 2^40 = 1024*2^30.


    // Assumes N >> 0.
    // Returns true if N is prime, and false otherwise.
    // Maybe slow for those large N where N-1 is hard to factorize.
    bool LucasTest(const BigInt& N)
    {
      if (N <= 0) CoCoA_THROW_ERROR(ERR::ReqPositive, "LucasTest(N)");
      if (N == 1) return false;
      if (N == 2) return true;
      const factorization<BigInt> facpows = factor(N-1);
      const vector<BigInt>& primes = facpows.myFactors();
      const int NumPrimes = len(primes); // overflow possible???
      using std::floor;
      const int Amax = static_cast<int>(floor(2*log(N)*log(N))); // experimentally checked up to 10^9: a better bound is 2.62*log(N)*log(log(N)) for N > 3
      for (int a=2; a <= Amax; ++a)
      {
        if (a==4 || a==8 || a==9 || a==16) continue; // skip some pure powers
        BigInt pwr;
        for (int i=0; i < NumPrimes; ++i)
        {
          pwr = PowerMod(a, (N-1)/primes[i], N);
          if (pwr == 1) break;
        }
        if (pwr != 1) return true;
      }
      return false;
    }

  } // end of anonymous namespace


  SmallPrime::SmallPrime(long p):
      myVal(p)
  {
    if (p < 1 || !IsPrime(p))
      CoCoA_THROW_ERROR(ERR::BadArg, "SmallPrime ctor");
  }


  std::ostream& operator<<(std::ostream& out, SmallPrime p)
  {
    if (!out) return out;  // short-cut for bad ostreams
    return out << long(p);
  }


  //------------------------------------------------------------------
  // Sieve of Eratosthenes: two versions: one starts from 1, the other takes a range
  
  std::vector<bool> eratosthenes(const MachineInt& n)
  {
    if (IsNegative(n) || IsZero(n)) CoCoA_THROW_ERROR(ERR::ReqPositive, "eratosthenes");
    if (!IsSignedLong(n)) CoCoA_THROW_ERROR(ERR::ArgTooBig, "eratosthenes");
    const long N = AsSignedLong(n)/2;
    vector<bool> sieve(N+1,true);
    sieve[0] = false; // 1 is not prime
    const long imax = FloorSqrt(N/2);
    for (long i=1; i <= imax; ++i)
    {
      if (!sieve[i]) continue;
      const long p = 2*i+1;  // cannot overflow (assuming MAXLONG is odd)
      for (long j=(p*p)/2; j <= N; j+=p) // (p*p)/2 int div; j+=p cannot overflow!
        sieve[j] = false;
    }
    return sieve; // wasteful copy
  }


  vector<bool> EratosthenesRange(const MachineInt& LWB, const MachineInt& UPB)
  {
    if (IsNegative(LWB) || IsNegative(UPB)) CoCoA_THROW_ERROR(ERR::ReqNonNegative, "EratosthenesRange");
    if (!IsSignedLong(LWB) || !IsSignedLong(UPB)) CoCoA_THROW_ERROR(ERR::ArgTooBig, "EratosthenesRange");
    long lwb = AsSignedLong(LWB);
    long upb = AsSignedLong(UPB);
    if (upb <= lwb) CoCoA_THROW_ERROR(ERR::BadArg, "EratosthenesRange");
    if (IsEven(lwb)) ++lwb;
    if (IsEven(upb)) ++upb;
    const long width = 1+(upb-lwb)/2; // exact division
    vector<bool> sieve(width, true); // slot k corr to integer lwb+2*k
    PrimeSeq PS;
    // NB first check in loop skips past 2
    while (!IsEnded(++PS))
    {
      const long p = CurrPrime(PS);
      if (p > upb/p) break;

      long k = (lwb-1)/p;
      if (k < p) k = p;
      else k += (IsEven(k))?1:2;
      long index = (p*k-lwb)/2; // index of first odd num greater than LWB div by p
      while (index < width)
      {
        sieve[index] = false;
        index += p;
      }
    }
    return sieve;
  }


  //------------------------------------------------------------------
  // PrimeSeqForCRT

  PrimeSeqForCRT::PrimeSeqForCRT():
      myCurrPrime(0),
      myIndex(0),
      myNearlyPrimeSeq(PrevPrime(ourTblStart + ourSieveRange))
  {
    InitTbl();
    myCurrPrime = ourTblStart + 2*ourPrimeDiffTbl[0];
  }

  // static data members
  const long PrimeSeqForCRT::ourTblStart = (1L << (numeric_limits<unsigned long>::digits/2-2))+1; // MUST BE ODD!!
  const long PrimeSeqForCRT::ourSieveRange = 1048576;
  int PrimeSeqForCRT::ourTblSize = 0; // will be set by PrimeSeqForCRT::InitTbl
  vector<unsigned char> PrimeSeqForCRT::ourPrimeDiffTbl; // will be set by PrimeSeqForCRT::InitTbl
  long PrimeSeqForCRT::ourLastPrimeInTbl = 0; // proper value set by InitTbl


  // BUG!  THIS NEEDS TO BE THREADSAFE
  long PrimeSeqForCRT::InitTbl()
  {
    if (!ourPrimeDiffTbl.empty()) return ourLastPrimeInTbl; // BUG not properly threadsafe!

    const vector<bool> sieve = EratosthenesRange(ourTblStart, ourTblStart+ourSieveRange);
    // IMPORTANT: we can use unsigned char because up to 2^32 the largest
    // prime gap is 2*168=336 between 3842610773 and 3842611109.
    vector<unsigned char> DiffTbl; // DiffTbl.reserve(?????);
    const long n = len(sieve);
    long LastVal = 0;
    for (int i=0; i < n; ++i)
      if (sieve[i]) { DiffTbl.push_back(i-LastVal); LastVal = i; }

    ourLastPrimeInTbl = ourTblStart + 2*LastVal;
    swap(ourPrimeDiffTbl, DiffTbl);
    ourTblSize = len(ourPrimeDiffTbl);
    return ourLastPrimeInTbl;
  }


  SmallPrime PrimeSeqForCRT::operator*() const
  {
    if (IamEnded()) CoCoA_THROW_ERROR(ERR::IterEnded, "PrimeSeqForCRT::op*");
    return SmallPrime(myCurrPrime, ArgIsPrime);
  }


  PrimeSeqForCRT& PrimeSeqForCRT::operator++()
  {
    if (IamEnded()) CoCoA_THROW_ERROR(ERR::IterEnded, "PrimeSeqForCRT::op++");
    ++myIndex;
    if (myIndex >= ourTblSize) myCurrPrime = NextPrime(myNearlyPrimeSeq);
    else myCurrPrime += 2*ourPrimeDiffTbl[myIndex];
    return *this;
  }

  PrimeSeqForCRT PrimeSeqForCRT::operator++(int)
  {
    if (IamEnded()) CoCoA_THROW_ERROR(ERR::IterEnded, "PrimeSeqForCRT::op++");
    PrimeSeqForCRT copy(*this);
    operator++();
    return copy;
  }


  std::ostream& operator<<(std::ostream& out, const PrimeSeqForCRT& PSeq)
  {
    if (!out) return out;  // short-cut for bad ostreams
//    if (IamEnded) return out << "PrimeSeqForCRT(ENDED)";
    out << "PrimeSeqForCRT(curr=" << *PSeq << ")";
    return out;
  }
  

  //------------------------------------------------------------------


  std::vector<unsigned char> FastFinitePrimeSeq::ourPrimeDiffTbl;
  long FastFinitePrimeSeq::ourTblSize = 0;
  const long FastFinitePrimeSeq::ourSieveSize = 1048576; // 1048576 gives a diff table of size about 82K
                                                         // 16290000 gives a diff table size about 1Mbyte
  long FastFinitePrimeSeq::ourLastPrime = 0; // proper value set by FastFinitePrimeSeq::InitTbl

  FastFinitePrimeSeq::FastFinitePrimeSeq():
      myCurrPrime(2),
      myIndex(0)
  {
    InitTbl();
  }

  // BUG!  THIS NEEDS TO BE THREADSAFE (or called by GlobalManager)
  void FastFinitePrimeSeq::InitTbl()
  {
    if (!ourPrimeDiffTbl.empty()) return;
    // Fill ourPrimeTbl: first make sieve, then fill ourPrimeTbl
    const vector<bool> sieve = eratosthenes(ourSieveSize);

    vector<unsigned char> PrimeDiffTbl; //??PrimeDiffTbl.reserve(82024);
    const long n = len(sieve);
    long LastVal = 1;
    for (int i=2; i < n; ++i)
      if (sieve[i])
      {
        PrimeDiffTbl.push_back(i-LastVal);// i-LastVal < 256, so narrowing is safe!
        LastVal = i;
      }

    ourLastPrime = 2*LastVal+1;
    swap(ourPrimeDiffTbl, PrimeDiffTbl); // really assignment
    ourTblSize = len(ourPrimeDiffTbl);
  }


  SmallPrime FastFinitePrimeSeq::operator*() const
  {
    CoCoA_ASSERT(!IamEnded());
//    if (!IsEnded(myPrimeSeq16)) return CurrPrime(myPrimeSeq16);
    return SmallPrime(myCurrPrime, ArgIsPrime);
  }

  FastFinitePrimeSeq& FastFinitePrimeSeq::operator++()
  {
    CoCoA_ASSERT(!IamEnded());
//    if (IamEnded()) CoCoA_THROW_ERROR(ERR::IterEnded, "PrimeSeq::op++");
    if (myCurrPrime == 2)
     myCurrPrime = 3;
    else
      myCurrPrime += 2*ourPrimeDiffTbl[myIndex++];

    return *this;
  }

  FastFinitePrimeSeq FastFinitePrimeSeq::operator++(int)
  {
    FastFinitePrimeSeq copy(*this);
    operator++();
    return copy;
  }


  bool FastFinitePrimeSeq::IamEnded() const noexcept
  {
    return (myIndex >= ourTblSize);
  }


  std::ostream& operator<<(std::ostream& out, const FastFinitePrimeSeq& seq)
  {
    if (!out) return out;  // short-cut for bad ostreams
    if (IsEnded(seq)) return out << "FastFinitePrimeSeq(ENDED)";
    out << "FastFinitePrimeSeq(curr=" << *seq << ")";
    return out;
  }

  //------------------------------------------------------------------
  // PrimeSeq
  

  PrimeSeq::PrimeSeq()
  {}
      
  
  SmallPrime PrimeSeq::operator*() const
  {
//    if (!IsEnded(myPrimeSeq16)) return CurrPrime(myPrimeSeq16);
    return SmallPrime(*myMostlyPrimeSeq, ArgIsPrime);
  }

  PrimeSeq& PrimeSeq::operator++()
  {
//    if (IamEnded()) CoCoA_THROW_ERROR(ERR::IterEnded, "PrimeSeq::op++");
    ++myMostlyPrimeSeq;
    if (myMostlyPrimeSeq.IamUsingPrimeTbl()) return *this;
    while (!IsPrime(*myMostlyPrimeSeq))
      ++myMostlyPrimeSeq;
    return *this;
  }

  PrimeSeq PrimeSeq::operator++(int)
  {
//    if (IamEnded()) CoCoA_THROW_ERROR(ERR::IterEnded, "PrimeSeq::op++");
    PrimeSeq copy(*this);
    operator++();
    return copy;
  }


  std::ostream& operator<<(std::ostream& out, const PrimeSeq& PSeq)
  {
    if (!out) return out;  // short-cut for bad ostreams
//    if (IsEnded(PSeq)) return out << "PrimeSeq(ENDED)";
    out << "PrimeSeq(curr=" << *PSeq << ")";
    return out;
  }
  


  //------------------------------------------------------------------
  // NoSmallFactorSeq

  // Check that n >= 0 and return next number which is not div by 2,3,5,7
  long CheckCtorArg(const MachineInt& n)
  {
    if (IsNegative(n) || !IsSignedLong(n))  CoCoA_THROW_ERROR(ERR::BadArg, "NoSmallFactorSeq ctor arg");
    long N = AsSignedLong(n);
    if (N > numeric_limits<long>::max() - 10) CoCoA_THROW_ERROR(ERR::ArgTooBig, "NoSmallFactorSeq ctor arg");
    static const unsigned char skip[30] = {1,0,5,4,3,2,1,0,3,2,1,0,1,0,3,2,1,0,1,0,3,2,1,0,5,4,3,2,1,0};
    N += skip[N%30];
    if (N%7 == 0) { ++N; N += skip[N%30]; }
    CoCoA_ASSERT(gcd(N,210)==1);
    return N;
  }

  unsigned char CalcIndex30(unsigned char mod30)
  {
    switch (mod30)
    {
    case 1: return 0;
    case 7: return 1;
    case 11: return 2;
    case 13: return 3;
    case 17: return 4;
    case 19: return 5;
    case 23: return 6;
    case 29: return 7;
    }
    CoCoA_THROW_ERROR(ERR::ShouldNeverGetHere, "CalcIndex");
    return 0; // just to keep compiler quiet
  }

  unsigned char CalcIndex210(unsigned char mod210)
  {
    // COULD REPLACE THIS TABLE BY ONE HALF THE SIZE (since tbl is symmetric)
    // 99 below means this entry should never be used.
    static unsigned char IndexTbl[210] =
      {99, 0, 99, 99, 99, 99, 99, 99, 99, 99, 99, 1, 99, 2, 99, 99, 99, 3, 99, 4, 99,
       99, 99, 5, 99, 99, 99, 99, 99, 6, 99, 7, 99, 99, 99, 99, 99, 8, 99, 99, 99, 9,
       99, 10, 99, 99, 99, 11, 99, 99, 99, 99, 99, 12, 99, 99, 99, 99, 99, 13, 99, 14, 99,
       99, 99, 99, 99, 15, 99, 99, 99, 16, 99, 17, 99, 99, 99, 99, 99, 18, 99, 99, 99, 19,
       99, 99, 99, 99, 99, 20, 99, 99, 99, 99, 99, 99, 99, 21, 99, 99, 99, 22, 99, 23, 99,
       99, 99, 24, 99, 25, 99, 99, 99, 26, 99, 99, 99, 99, 99, 99, 99, 27, 99, 99, 99, 99,
       99, 28, 99, 99, 99, 29, 99, 99, 99, 99, 99, 30, 99, 31, 99, 99, 99, 32, 99, 99, 99,
       99, 99, 33, 99, 34, 99, 99, 99, 99, 99, 35, 99, 99, 99, 99, 99, 36, 99, 99, 99, 37,
       99, 38, 99, 99, 99, 39, 99, 99, 99, 99, 99, 40, 99, 41, 99, 99, 99, 99, 99, 42, 99,
       99, 99, 43, 99, 44, 99, 99, 99, 45, 99, 46, 99, 99, 99, 99, 99, 99, 99, 99, 99, 47};

    CoCoA_ASSERT(/*mod210 >= 0 &&*/ mod210 < 210 && IndexTbl[mod210] != 99);
    return IndexTbl[mod210]; // just to keep compiler quiet
  }


  NoSmallFactorSeq::NoSmallFactorSeq(const MachineInt& StartVal):
      myCurrVal(CheckCtorArg(StartVal)),
      myIndex(CalcIndex210(myCurrVal%210)),
      myValMod11(myCurrVal%11),
      myValMod13(myCurrVal%13),
      myValMod17(myCurrVal%17),
      myValMod19(myCurrVal%19),
      myValMod23(myCurrVal%23)
  {
    if (myValMod11 == 0 || myValMod13 == 0 || myValMod17 == 0 || myValMod19 == 0 || myValMod23 == 0)
      operator++();
  }

  long NoSmallFactorSeq::operator*() const
  {
//    if (IamEnded()) CoCoA_THROW_ERROR(ERR::IterEnded, "PrimeSeq::op*");
    return myCurrVal;
  }


  NoSmallFactorSeq& NoSmallFactorSeq::operator++()
  {
//    if (IamEnded()) CoCoA_THROW_ERROR(ERR::IterEnded, "PrimeSeq::op++");
    static const unsigned char skip[48] = {10,  2,  4,  2,  4,  6,  2,  6,  4,  2,  4,  6,  6,  2,  6,  4,  2,  6,  4,  6,  8,  4,  2,  4,  2,  4,  8,  6,  4,  6,  2,  4,  6,  2,  6,  6,  4,  2,  4,  6,  2,  6,  4,  2,  4,  2,  10,2};
    do
    {
      const unsigned char delta = skip[myIndex];
      ++myIndex; if (myIndex == 48) myIndex = 0;
      myCurrVal += delta; // BUG BUG BUG check for overflow!!!
      myValMod11 += delta; if (myValMod11 >= 11) myValMod11 -= 11;
      myValMod13 += delta; if (myValMod13 >= 13) myValMod13 -= 13;
      myValMod17 += delta; if (myValMod17 >= 17) myValMod17 -= 17;
      myValMod19 += delta; if (myValMod19 >= 19) myValMod19 -= 19;
      myValMod23 += delta; if (myValMod23 >= 23) myValMod23 -= 23;
    }
    while (myValMod11 == 0 || myValMod13 == 0 || myValMod17 == 0 || myValMod19 == 0 || myValMod23 == 0);
    return *this;
  }

  NoSmallFactorSeq NoSmallFactorSeq::operator++(int)
  {
//    if (IamEnded()) CoCoA_THROW_ERROR(ERR::IterEnded, "PrimeSeq::op++");
    NoSmallFactorSeq copy(*this);
    operator++();
    return copy;
  }



//   void NoSmallFactorSeq::myStartFrom(long n)
//   {
//       CoCoA_ASSERT(IsCoprime(n,30));
//       myCurrVal = CheckCtorArg(n);
//       myIndex = CalcIndex(n);
// //      myValMod7 = myCurrVal%7;
//       myValMod11 = myCurrVal%11;
//       myValMod13 = myCurrVal%13;
//       myValMod17 = myCurrVal%17;
//       myValMod19 = myCurrVal%19;
//       myValMod23 = myCurrVal%23;
//     }


  SmallPrime NextPrime(NoSmallFactorSeq& seq)
  {
    do { ++seq; }
    while (!IsPrime(*seq));
    return SmallPrime(*seq, ArgIsPrime);
  }

  std::ostream& operator<<(std::ostream& out, const NoSmallFactorSeq& seq)
  {
    if (!out) return out;  // short-cut for bad ostreams
//    if (IsEnded(seq)) return out << "NoSmallFactorSeq(ENDED)";
    out << "NoSmallFactorSeq(curr=" << *seq << ")";
    return out;
  }


  //------------------------------------------------------------------

  FastMostlyPrimeSeq::FastMostlyPrimeSeq():
      myPrimeSeq(),
      myMostlyPrimeSeq(myPrimeSeq.myLastPrime()) // first elem is duplicate of last prime -- this is wanted!
  {}
    

  long FastMostlyPrimeSeq::operator*() const
  { if (!IsEnded(myPrimeSeq)) return *myPrimeSeq; else return *myMostlyPrimeSeq; }


  FastMostlyPrimeSeq& FastMostlyPrimeSeq::operator++()
  {
    if (!IsEnded(myPrimeSeq)) ++myPrimeSeq;
    else ++myMostlyPrimeSeq;
    return *this;
  }

  FastMostlyPrimeSeq FastMostlyPrimeSeq::operator++(int)
  {
    FastMostlyPrimeSeq copy(*this);
    operator++();
    return copy;
  }


  //------------------------------------------------------------------

  // This definition is for 64-bit machines, but remains valid for 32-bitters
  // ASSUMES n < 2^32  (automatically true on 32-bitters).
  // I do not recall the reference for this definition (H.Cohen's book?); verified independently (2018-03-09).
  bool IsPrime(const MachineInt& mi)
  {
    if (IsZero(mi) || IsNegative(mi))
      CoCoA_THROW_ERROR(ERR::ReqPositive, "IsPrime(n)");

    const unsigned long n = AsUnsignedLong(mi);
    if (n <= 37)  return IsSmallPrime37(n);
    if (HasSmallFactor37(n))  return false;
    if (n < 1681)  return true; // 1681 = square(NextPrime(37))
    // At this point n >= 1681, and n has no factor <= 37.
    if (n <= MaxSquarableInteger<unsigned long>())
      return IsPrime_squarable(n);
    return IsPrime_NotSquarable(n);
  }


  // If value happens to be small we use IsPrime(MachineInt)
  // o/w IsProbPrime followed by Lucas test.
  // HINT: should also make a version of StrongPseudoPrime64 for BigInt?!?
  bool IsPrime(const BigInt& N)
  {
    if (N <= 0) CoCoA_THROW_ERROR(ERR::ReqPositive, "IsPrime(N)");
    unsigned long n; // since we know N > 0
    if (IsConvertible(n, N))
      return IsPrime(n);

    CoCoA_ASSERT(N > numeric_limits<unsigned long>::max());
    if (!IsProbPrime(N)) return false;
    return LucasTest(N);
  }


  // "Probable prime" test.
  // According to GMP documentation uses Miller-Rabin (after a few trial divisions).
  bool IsProbPrime(const MachineInt& n) { return IsProbPrime(n, ProbPrimeIters); }
  bool IsProbPrime(const MachineInt& n, const MachineInt& NumIters)
  {
    if (IsZero(n) || IsNegative(n))
      CoCoA_THROW_ERROR(ERR::ReqPositive, "IsProbPrime(n,NumIters):  1st arg");
    if (IsZero(NumIters) || IsNegative(NumIters))
      CoCoA_THROW_ERROR(ERR::ReqPositive, "IsProbPrime(n,NumIters):  2nd arg");
    // Just call the slow version -- at least will guarantee coherent behaviour
    return IsProbPrime(BigInt(n), NumIters);
  }

  bool IsProbPrime(const BigInt& N) { return IsProbPrime(N, ProbPrimeIters); }
  bool IsProbPrime(const BigInt& N, const MachineInt& NumIters)
  {
    if (N <= 0)
      CoCoA_THROW_ERROR(ERR::ReqPositive, "IsProbPrime(N,NumIters): 1st arg");
    if (IsZero(NumIters) || IsNegative(NumIters))
      CoCoA_THROW_ERROR(ERR::ReqPositive, "IsProbPrime(N,NumIters): 2nd arg");
    return mpz_probab_prime_p(mpzref(abs(N)), AsUnsignedLong(NumIters));
  }



  SmallPrime NextPrime(const MachineInt& mi)
  {
    if (IsNegative(mi)) CoCoA_THROW_ERROR(ERR::ReqNonNegative, "NextPrime(n)");
    if (!IsSignedLong(mi)) CoCoA_THROW_ERROR(ERR::ArgTooBig, "NextPrime(n)");
    long n = AsSignedLong(mi);

    // Special cases if n < 5
    if (n < 2) return SmallPrime(2, ArgIsPrime);
    if (n == 2) return SmallPrime(3, ArgIsPrime);
    if (n < 5) return SmallPrime(5, ArgIsPrime);

    const long MaxLong = std::numeric_limits<long>::max();
    int n30 = n%30;
    while (true)
    {
      const int delta = skip[n30];
      if (n > MaxLong - delta) break; // break if n+delta would overflow
      n += delta;
      n30 += delta;
      if (n30 >= 30) n30 -= 30;
      if (n <= 37) { if (IsSmallPrime37(n)) return SmallPrime(n, ArgIsPrime); }
      else if (IsPrime(n)) return SmallPrime(n, ArgIsPrime);
    }
    // Reach here only if "overflow" has occurred.
    return SmallPrime(0, ArgIsPrime); // to signify "overflow"
  }


  SmallPrime PrevPrime(const MachineInt& mi)
  {
    if (IsNegative(mi)) CoCoA_THROW_ERROR(ERR::ReqNonNegative, "PrevPrime(n)");
    if (!IsSignedLong(mi)) CoCoA_THROW_ERROR(ERR::ArgTooBig, "PrevPrime(n)");
    long n = AsSignedLong(mi);

    // Special cases if n < 8
    if (n < 3) CoCoA_THROW_ERROR(ERR::OutOfRange, "PrevPrime(n): no previous prime");
    if (n == 3) return SmallPrime(2, ArgIsPrime);
    if (n < 6) return SmallPrime(3, ArgIsPrime);
    if (n < 8) return SmallPrime(5, ArgIsPrime);

    int n30 = n%30;
    do
    {
      const int delta = fall[n30];
      n -= delta;
      n30 -= delta;
      if (n30 < 0) n30 += 30;
    } while (!IsPrime(n));
    return SmallPrime(n, ArgIsPrime);
  }



  // Uniform random distribution on primes from 5 to MAX (incl).
  // Primes 2 and 3 are excluded!!

  // HINT: could also do RandomSmallPrime over a range provided that
  // the range is wide enough.  A quick search shows that the largest
  // gap between primes is at most 0.37*FloorLog2(p)^2
  // (with two exceptions being p=23 and p=113).  See on-line somewhere
  // Terry Tao comments that it is conjectured that the largest prime gaps
  // grow as const*log(p)^2.
  SmallPrime RandomSmallPrime(const MachineInt& MAX)
  {
    static const unsigned char shift30[8] = {1, 7, 11, 13, 17, 19, 23, 29};
    if (IsNegative(MAX) || !IsSignedLong(MAX)) CoCoA_THROW_ERROR(ERR::BadArg, "RandomSmallPrime");
    const long UPB = AsSignedLong(MAX);
    if (UPB < 5) CoCoA_THROW_ERROR(ERR::BadArg, "RandomSmallPrime");
    
    // Arbitrarily impose limit  MAX <= 2^31-1
    if (numeric_limits<long>::radix == 2 && numeric_limits<long>::digits > 32 && UPB > 2147483647L)
      CoCoA_THROW_ERROR(ERR::ArgTooBig, "RandomSmallPrime");

    // Divide interval into blocks of width 30; in each block consider only those numbers
    // not divisible by 2,3,5.  Handle the block from 0..29 specially.
    const long UPBOver30 = UPB/30;  // integer division!
    while (true)
    {
      const long candidate = 30*RandomLong(0, UPBOver30) + shift30[RandomLong(0,7)];
      if (candidate > UPB) continue;
      if (candidate > 30)
      {
        if (IsPrime(candidate)) return SmallPrime(candidate, ArgIsPrime);
        continue;
      }
      // Special Case: candidate < 30
      if (candidate == 1)
        return SmallPrime(5, ArgIsPrime);
      else
        return SmallPrime(candidate, ArgIsPrime);
    }
  }

  SmallPrime RandomNBitPrime(const MachineInt& N)
  {
    static const unsigned char shift30[8] = {1, 7, 11, 13, 17, 19, 23, 29};
    if (IsNegative(N) || !IsSignedLong(N))
      CoCoA_THROW_ERROR(ERR::BadArg, "RandomNBitPrime: arg must be between 10 and 31");
    const long n = AsSignedLong(N);
    if (n < 10 || n > 31)
      CoCoA_THROW_ERROR(ERR::BadArg, "RandomNBitPrime: arg must be between 10 and 31");
    
    const long pwr2 = 1L << (n-1);
    // Divide interval into blocks of width 30; in each block consider only those numbers
    // not divisible by 2,3,5.  Handle the block from 0..29 specially.
    const long nblocks = pwr2/30;  // integer division!
    const long shift = 30*nblocks;
    while (true)
    {
      const long candidate = shift + 30*RandomLong(0, nblocks) + shift30[RandomLong(0,7)]; // BUG??? OVERFLOW??? use ulong???
      if (candidate < pwr2 || candidate/2 > pwr2)  continue; // really necessary?
      if (IsPrime(candidate)) return SmallPrime(candidate, ArgIsPrime);
    }
  }


  BigInt NextProbPrime(const BigInt& N) { return NextProbPrime(N, ProbPrimeIters); }
  BigInt NextProbPrime(BigInt N, const MachineInt& NumIters)
  {
    if (N < 0)
      CoCoA_THROW_ERROR(ERR::ReqNonNegative, "NextProbPrime(N,NumIters): 1st arg");
    if (IsZero(NumIters) || IsNegative(NumIters))
      CoCoA_THROW_ERROR(ERR::ReqPositive, "NextProbPrime(N,NumIters): 2nd arg");

    // if N is small use NextPrime
    if (N < 2147483647L) // magic num is largest prime below 2^31, sure to fit into a long
      return BigInt(NextPrime(ConvertTo<long>(N)));

/////?????    SignalWatcher MonitorSIGINT(SIGINT);
    int N30 = N%30;
    do
    {
      CheckForInterrupt("NextProbPrime");
      const int delta = skip[N30];
      N += delta;
      N30 += delta;
      if (N30 >= 30) N30 -= 30;
    } while (!IsProbPrime(N, NumIters));
    return N;
  }

  BigInt PrevProbPrime(const BigInt& N) { return PrevProbPrime(N, ProbPrimeIters); }
  BigInt PrevProbPrime(BigInt N, const MachineInt& NumIters)
  {
    if (N < 0)
      CoCoA_THROW_ERROR(ERR::ReqNonNegative, "PrevProbPrime(N,NumIters): 1st arg");
    if (IsZero(NumIters) || IsNegative(NumIters))
      CoCoA_THROW_ERROR(ERR::ReqPositive, "PrevProbPrime(N,NumIters): 2nd arg");

    // if N is small use PrevPrime
    if (N <= 2147483647L)
      return BigInt(PrevPrime(ConvertTo<long>(N)));

/////?????    SignalWatcher MonitorSIGINT(SIGINT);
    int N30 = N%30;
    do
    {
      CheckForInterrupt("PrevProbPrime");
      const int delta = fall[N30];
      N -= delta;
      N30 -= delta;
      if (N30 < 0) N30 += 30;
    } while (!IsProbPrime(N, NumIters));
    return N;
  }


  // -------------------------------------------------------
  // Estimate of number of primes up to given (small) limit
  
  // Source https://math.stackexchange.com/questions/264544/how-to-find-number-of-prime-numbers-up-to-to-n
  // 2023-12-10 Not sure how accurate the estimate is.
  long EstNumPrimesUpto(long n)
  {
    if (n < 2)  CoCoA_THROW_ERROR("arg must be >= 2", "EstNumPrimesUpto");
    static const long tbl[23] = {0,0,1,2,2,3,3,4,4,4,4,5,5,6,6,6,6,7,7,8,8,8,8};
    if (n < 23) return tbl[n];
    return static_cast<long>(1+std::floor(1.0121*n/(std::log(n)-1)));
  }


} // end of namespace CoCoA
