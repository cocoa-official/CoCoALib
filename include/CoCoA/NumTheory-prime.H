#ifndef CoCoA_NumTheory_prime_H
#define CoCoA_NumTheory_prime_H

//   Copyright (c)  2009,2012,2018  John Abbott,  Anna M. Bigatti

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


#include <iosfwd>
//using std::ostream;
#include <vector>
// using std::vector;

namespace CoCoA
{
  
  class BigInt;     // fwd decl -- defined in BigInt.H
  class MachineInt; // fwd decl -- defined in MachineInt.H

  //------------------------------------------------------------------
  // Functions to do with generating and testing primes

  std::vector<bool> eratosthenes(const MachineInt& n);
  std::vector<bool> EratosthenesRange(const MachineInt& lwb, const MachineInt& upb);


  class FastFinitePrimeSeq; // fwd decl -- defined below
  class NoSmallFactorSeq;   // fwd decl -- defined below

  // --------------------------------------------
  // Class SmallPrime & related enum

  // SmallPrime: value is either 0 or a positive prime which fits into a long
  enum class SmallPrimeFlag { ArgIsPrime }; // used only by SmallPrime ctor
  constexpr SmallPrimeFlag ArgIsPrime = SmallPrimeFlag::ArgIsPrime; // safely inject into namespace CoCoA

  class SmallPrime
  {
  public:
    explicit SmallPrime(long p);  // WARNING!!  RELATIVELY SLOW (because checks for primality; throws if not!)
    SmallPrime(long p, SmallPrimeFlag) noexcept: myVal(p) { /*Allow also p==0; see NextPrime overflow*/ } // Caller's responsibility to ensure that p is prime (or 0)!
    operator long() const noexcept { return myVal; }
  private: // data member
    long myVal; // Either 0 or a (positive) prime
  };

  std::ostream& operator<<(std::ostream& out, SmallPrime p);
  

  //------------------------------------------------------------------
  // Prime seq starting from 2, but only limited length (reads from a table).
  // op++ is fast (just incr an index).
  
  class FastFinitePrimeSeq
  {
  public:
    FastFinitePrimeSeq();
    SmallPrime operator*() const; // NOT noexcept -- throws if IsEnded
    FastFinitePrimeSeq& operator++();
    FastFinitePrimeSeq operator++(int);
    bool IamEnded() const noexcept;
    long myLastPrime() const noexcept { return ourLastPrime; } // should be called only after ctor!!!
    static void InitTbl(); // this can be called by ctor for GlobalManager
  private: // data members
    long myCurrPrime;
    long myIndex;
    static long ourLastPrime;
    static const long ourSieveSize;
    static long ourTblSize;
    static std::vector<unsigned char> ourPrimeDiffTbl;  // all half-differences < 256 (first gap greater than 512 is between 304599509051-2*257 and 304599509051)
  };

  inline bool IsEnded(const FastFinitePrimeSeq& seq) noexcept
  { return seq.IamEnded(); }

  std::ostream& operator<<(std::ostream& out, const FastFinitePrimeSeq& seq);


  //------------------------------------------------------------------
  // Seq of positive ints having no small factor (see ourMaxPrime)
  // op++ is reasonably quick.  Used by SmoothFactor and IsSqFree (see NumTheory.C)

  class NoSmallFactorSeq
  {
  public:
    explicit NoSmallFactorSeq(const MachineInt& StartVal);
    long operator*() const; // NOT noexcept -- throws if IsEnded
    NoSmallFactorSeq& operator++();
    NoSmallFactorSeq operator++(int);
    bool IamEnded() const noexcept { return false; /*myCurrVal > std::numeric_limits<long>::max() - 1000;*/ }
//???    void myStartFrom(long n);
    static long ourMaxPrime() { return 23; }  // ???? constexpr ????
  private: // data members
    long myCurrVal;
    unsigned char myIndex;
//    unsigned char myValMod7;
    unsigned char myValMod11;
    unsigned char myValMod13;
    unsigned char myValMod17;
    unsigned char myValMod19;
    unsigned char myValMod23;
  };

  SmallPrime NextPrime(NoSmallFactorSeq& seq);

  inline bool IsEnded(const NoSmallFactorSeq& seq) noexcept { return seq.IamEnded(); }

  std::ostream& operator<<(std::ostream& out, const NoSmallFactorSeq& seq);


  //------------------------------------------------------------------
  // Seq of positive ints; seq contains all primes, and a "few" composites.
  // op++ is reasonably fast
  
  class FastMostlyPrimeSeq
  {
  public:
    FastMostlyPrimeSeq();
    long operator*() const; // NOT noexcept -- throws if IsEnded
    FastMostlyPrimeSeq& operator++();
    FastMostlyPrimeSeq operator++(int);
    bool IamUsingPrimeTbl() const noexcept { return !IsEnded(myPrimeSeq); }
    bool IamEnded() const noexcept { return IsEnded(myPrimeSeq) && IsEnded(myMostlyPrimeSeq); }
  private: // data members
    FastFinitePrimeSeq myPrimeSeq;     // use this until it runs out...
    NoSmallFactorSeq myMostlyPrimeSeq; // ... then use this.
  };

  inline bool IsEnded(const FastMostlyPrimeSeq& seq) noexcept
  { return seq.IamEnded(); }


  
  //------------------------------------------------------------------
  // PrimeSeqForCRT starts from a large number; contrast with PrimeSeq.
  // First "few" calls to op++ are fast (from a table), then slower.
  
  class PrimeSeqForCRT
  {
  public:
    PrimeSeqForCRT();
    SmallPrime operator*() const; // NOT noexcept -- throws if IamEnded
    PrimeSeqForCRT& operator++();
    PrimeSeqForCRT operator++(int);
    bool IamEnded() const noexcept { return myCurrPrime == 0; }
    static long InitTbl(); // this can be called by ctor for GlobalManager
  private: // data members
    long myCurrPrime;
    int myIndex;  // if greater than ourTblSize then we use fn NextPrime(.)
    NoSmallFactorSeq myNearlyPrimeSeq;
    static const long ourTblStart;
    static const long ourSieveRange;
    static long ourLastPrimeInTbl;
    static int ourTblSize;
    static std::vector<unsigned char> ourPrimeDiffTbl; // half-differences <= 250 < 256 for all primes up to 304599509051-2*257
  };

  inline SmallPrime CurrPrime(const PrimeSeqForCRT& PSeq)
  {
    return *PSeq;
  }
  
  inline SmallPrime NextPrime(PrimeSeqForCRT& PSeq)
  {
    return *++PSeq;
  }

  inline bool IsEnded(const PrimeSeqForCRT& PSeq) noexcept
  {
    return PSeq.IamEnded();
  }

  std::ostream& operator<<(std::ostream& out, const PrimeSeqForCRT& PSeq);


  //------------------------------------------------------------------
  // PrimeSeq starts from 2; contrast with PrimeSeqForCRT

  class PrimeSeq
  {
  public:
    PrimeSeq();
    SmallPrime operator*() const; // NOT noexcept -- throws if IamEnded
    PrimeSeq& operator++();
    PrimeSeq operator++(int);
    bool IamEnded() const { return IsEnded(myMostlyPrimeSeq); }
  private: // data members
    FastMostlyPrimeSeq myMostlyPrimeSeq;
  };

  inline SmallPrime CurrPrime(const PrimeSeq& PSeq)
  {
    return *PSeq;
  }
  
  inline SmallPrime NextPrime(PrimeSeq& PSeq)
  {
    return *++PSeq;
  }

  inline bool IsEnded(const PrimeSeq& PSeq)
  {
    return PSeq.IamEnded();
  }

  std::ostream& operator<<(std::ostream& out, const PrimeSeq& PSeq);



  bool IsPrime(const MachineInt& n);
  bool IsPrime(const BigInt& N); // WARNING: might be very slow!!
  bool IsProbPrime(const MachineInt& n); // probably not very useful
  bool IsProbPrime(const MachineInt& n, const MachineInt& NumIters); // probably not very useful
  bool IsProbPrime(const BigInt& N);
  bool IsProbPrime(const BigInt& N, const MachineInt& NumIters);

  SmallPrime NextPrime(const MachineInt& n); ///< next prime number after n or 0 if there is none, must have n>0
  SmallPrime PrevPrime(const MachineInt& n); ///< previous prime number before n or 0 if there is none, must have n>0
  SmallPrime RandomSmallPrime(const MachineInt& MAX);  ///< random prime from 5 to MAX (unif distrib)
  SmallPrime RandomNBitPrime(const MachineInt& n);     ///< random n-bit prime (may occasionally be just too small or too big)


  // !!No NextPrime for BigInt!!  Use NextProbPrime and call IsPrime if you want a (costly) guarantee of primality.
  BigInt NextProbPrime(const BigInt& N); ///< next probable prime number after N, requires N > 1
  BigInt NextProbPrime(BigInt N, const MachineInt& NumIters); ///< next probable prime number after N, requires N > 1
  BigInt PrevProbPrime(const BigInt& N); ///< previous probable prime number before N, requires N > 1
  BigInt PrevProbPrime(BigInt N, const MachineInt& NumIters); ///< previous probable prime number before N, requires N > 1


  long EstNumPrimesUpto(long n);  // estimate number of primes up to given limit (not sure how accurate)

} // end of namespace CoCoA

#endif
