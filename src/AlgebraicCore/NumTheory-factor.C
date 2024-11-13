//   Copyright (c)  1999,2009-2011,2019,2020  John Abbott and Anna M. Bigatti

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


#include <gmp.h>

#include "CoCoA/time.H"

#include "CoCoA/NumTheory-factor.H"
#include "CoCoA/NumTheory-gcd.H"
#include "CoCoA/NumTheory-CoprimeFactorBasis.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/assert.H"
#include "CoCoA/bool3.H"
#include "CoCoA/config.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/interrupt.H"
#include "CoCoA/utils.H"

#include <algorithm>
using std::min;
using std::sort;
#include <limits>
using std::numeric_limits;
//#include <vector>
using std::vector;

namespace CoCoA
{

  namespace // anonymous
  {

    long factor_TrialDiv_ArgCheckLimit(const MachineInt& TrialLimit/*, const char* const FnName*/)
    {
      static const long BiggestSmallPrime = PrevPrime(numeric_limits<long>::max());
      if (!IsSignedLong(TrialLimit) || AsSignedLong(TrialLimit) < 1)
        CoCoA_THROW_ERROR(ERR::BadArg, "factor_TrialDiv:  TrialLimit must be at least 1 and fit into a machine long");
      // Below Pmax is unsigned long so that the code will work even if input TrialLimit is numeric_limits<long>::max()
      return min(BiggestSmallPrime,  AsSignedLong(TrialLimit));
    }


    // Divide n by highest power of b possible (to get integer result).
    // b need not be prime.
    // Changes value of n; returns multiplicity of b in original value of n.
    // EXCEPTION SAFE
    long DivideOutMaxPower_NoArgCheck(long b, long& n)
    {
      CoCoA_ASSERT(n != 0);
      CoCoA_ASSERT(b > 1);

      long mult = 0;
      while (n%b == 0)
      {
        n /= b;
        ++mult;
      }
      return mult;
    }


    // Divide N by highest power of B possible (to get integer result).
    // B need not be prime.
    // Changes value of N; returns multiplicity of B in original value of N.
    // EXCEPTION SAFE UNLESS GMP RUNS OUT OF MEMORY
    // long DivideOutMaxPower_NoArgCheck_linear(const BigInt& B, BigInt& N)
    // {
    //   CoCoA_ASSERT(B > 1);
    //   CoCoA_ASSERT(!IsZero(N));
    //   long b;
    //   if (IsConvertible(b, B)) return DivideOutMaxPower_NoArgCheck(b,N);
    //   long mult = 0;
    //   BigInt r;
    //   while (true)
    //   {
    //     CheckForInterrupt("DivideOutMaxPower_NoArgCheck");
    //     quorem(N, r, N, B);
    //     if (!IsZero(r)) break;
    //     ++mult;
    //   }
    //   return mult;
    // }

    long DivideOutMaxPower_NoArgCheck(const BigInt& B, BigInt& N)
    {
      CoCoA_ASSERT(B > 1);
      CoCoA_ASSERT(!IsZero(N));
      if (B == 2)
      {
        const unsigned long m = mpz_scan1(mpzref(N), 0); // GMP can do this quickly
        mpz_tdiv_q_2exp(mpzref(N), mpzref(N), m);  // GMP can do this quickly
        return m;
      }

      long mult = 0;
      vector<BigInt> Bpwrs;  // Makes no difference:  Bpwrs.reserve(30);
      int k=0;
      Bpwrs.push_back(B);
      BigInt r;
      BigInt q;
      while (true)
      {
        CheckForInterrupt("DivideOutMaxPower_NoArgCheck (loop1)");
        const BigInt& Bpwr = Bpwrs.back();
        quorem(q, r, N, Bpwr);
        if (!IsZero(r)) break;
        swap(N,q); // really N = q;
        mult += 1<<k;
        Bpwrs.push_back(power(Bpwr,2));
        ++k;
      }

      while (!IsOne(N) && k > 0)
      {
        Bpwrs.pop_back();
        --k;
        CheckForInterrupt("DivideOutMaxPower_NoArgCheck (loop2)");
        const BigInt& Bpwr = Bpwrs.back();
        quorem(q, r, N, Bpwr);
        if (!IsZero(r)) { /*swap(N,r);*/ continue; }
        swap(N,q); // really N = q;
        mult += 1<<k;
//        Bpwrs.push_back(power(Bpwr,2));
      }
      return mult;
    }

    // Divide N by highest power of b possible (to get integer result).
    // b need not be prime.
    // Changes value of N; returns multiplicity of b in original value of N.
    // EXCEPTION SAFE UNLESS GMP RUNS OUT OF MEMORY
    long DivideOutMaxPower_NoArgCheck(long b, BigInt& N)
    {
      CoCoA_ASSERT(b > 1);
      if (b == 2)
      {
        const unsigned long m = mpz_scan1(mpzref(N), 0); // GMP can do this quickly
        mpz_tdiv_q_2exp(mpzref(N), mpzref(N), m);  // GMP can do this quickly
        return m;
      }

      return DivideOutMaxPower_NoArgCheck(BigInt(b), N);
      
      // long mult = 0;
      // {
      //   // Check if we can compute just with machine ints
      //   long n;
      //   if (IsConvertible(n, N))
      //   {
      //     mult = DivideOutMaxPower_NoArgCheck(b, n);
      //     N = n;
      //     return mult;
      //   }
      // }
      // // Case: N does not fit into a long
      // // Repeatedly divide by largest power of prime which fits in ulong:
      // long pwr = b;
      // int exp = 1;
      // const long limit = numeric_limits<long>::max()/b;
      // while (pwr <= limit)
      // {
      //   pwr *= b;
      //   ++exp;
      // }
      // // now pwr is largest power of base which fits into a long
      // BigInt q,r; // used for quotient and remainder in loop below
      // while (true)
      // {
      //   CheckForInterrupt("DivideOutMaxPower_NoArgCheck");
      //   quorem(q, r, N, pwr); // r could be unsigned long...
      //   if (!IsZero(r)) break;
      //   swap(N,q); // really just assignment N = q;
      //   mult += exp;
      // }
      // // N was not div by pwr, but might still be div by base...
      // long g = gcd(pwr, ConvertTo<long>(r)); // 0 < r < pwr
      // // Need to be a bit careful below since b need not be prime!
      // if (g < b) return mult;
      // const long extra = DivideOutMaxPower_NoArgCheck(b, g);
      // N /= SmallPower(b, extra);
      // return mult + extra;
    }



  } // end of namespace anonymous

  //------------------------------------------------------------------

  long radical(const MachineInt& n)
  {
    if (IsZero(n)) return 0;
    const factorization<long> FacInfo = factor(n);
    const vector<long>& facs = FacInfo.myFactors();
    long ans = 1;
    for (const long p: facs)
      ans *= p;
    return ans;
  }

  BigInt radical(const BigInt& N)
  {
    if (IsZero(N)) return N;
    if (IsProbPrime(abs(N))) return abs(N);
    const factorization<BigInt> FacInfo = factor(N);
    const vector<BigInt>& facs = FacInfo.myFactors();
    BigInt ans(1);
    for (const BigInt& p: facs)
      ans *= p;
    return ans;
  }
  

  // EXCLUDES most negative long!
  factorization<long> factor_TrialDiv(const MachineInt& n, const MachineInt& TrialLimit)
  {
    const long Pmax = factor_TrialDiv_ArgCheckLimit(TrialLimit);
    if (IsZero(n))
      CoCoA_THROW_ERROR(ERR::ReqNonZero, "factor_TrialDiv: arg 1");
    if (!IsSignedLong(n) || AsSignedLong(n) == numeric_limits<long>::min())
      CoCoA_THROW_ERROR(ERR::ArgTooBig, "factor_TrialDiv:  number to be factorized must fit into a machine long");
    long RemainingFactor = uabs(n); // cannot overflow since we exclude most neg long


    // Main loop: we simply do trial divisions by primes up to 31, & by numbers not divisible by 2,3,5,...,31.
    vector<long> factors;
    vector<long> exponents;

    // Use "mostly primes"; faster than using primes
    FastMostlyPrimeSeq TrialFactorSeq;
    
    long stop = Pmax; // highest prime we shall consider is "stop" (def'd below)
    if (RemainingFactor/stop < stop) stop = ConvertTo<long>(FloorSqrt(RemainingFactor));
    for (long p = *TrialFactorSeq; p <= stop; p = *++TrialFactorSeq)
    {
      if (RemainingFactor%p != 0) continue;
      const int exp = DivideOutMaxPower_NoArgCheck(p, RemainingFactor);
      CoCoA_ASSERT(exp > 0);
      factors.push_back(p);
      exponents.push_back(exp);
      if (RemainingFactor/stop < stop) stop = ConvertTo<long>(FloorSqrt(RemainingFactor));
    }
    // if RemainingFactor is non-triv & below limit, add it to the list of factors found.
    if (RemainingFactor > 1 && RemainingFactor <= Pmax)
    {
      factors.push_back(RemainingFactor);
      exponents.push_back(1);
      RemainingFactor = 1;
    }
    if (IsNegative(n))
      return factorization<long>(factors, exponents, -static_cast<long>(RemainingFactor));
    return factorization<long>(factors, exponents, RemainingFactor);
  }


  // This is very similar to the function above -- but I don't see how to share code.
  factorization<BigInt> factor_TrialDiv(const BigInt& N, const MachineInt& TrialLimit)
  {
    const long Pmax = factor_TrialDiv_ArgCheckLimit(TrialLimit);
    if (IsZero(N))
      CoCoA_THROW_ERROR(ERR::ReqNonZero, "factor_TrialDiv(N,TrialLimit): arg 1");
    {
      // if N fits into long: call fn above, then convert answer!
      long n;
      if (IsConvertible(n,N))
      {
        const auto facs = factor_TrialDiv(n, TrialLimit);
        factorization<BigInt> ans(BigInt(facs.myRemainingFactor()));
        for (int i=0; i < len(facs.myFactors()); ++i)
          ans.myAppend(BigInt(facs.myFactors()[i]), facs.myMultiplicities()[i]);
        return ans;
      }
    }

    // N is too large to fit into a long
    BigInt RemainingFactor = abs(N);

    // Main loop: we simply do trial divisions by primes up to 31 & then by numbers not divisible by 2,3,5,...,31.
    vector<BigInt> factors;
    vector<long> exponents;

    // Use "mostly primes"; faster than using primes
    FastMostlyPrimeSeq TrialFactorSeq;

    long stop = Pmax; // highest prime we shall consider is "stop"
    if (RemainingFactor/stop < stop) stop = ConvertTo<long>(FloorSqrt(RemainingFactor));
    long LogRemainingFactor = FloorLog2(RemainingFactor);
    long CountDivTests = 0;
    for (long p = *TrialFactorSeq; p <= stop; p = *++TrialFactorSeq)
    {
      CheckForInterrupt("factor_TrialDiv (main loop)");
      ++CountDivTests;
      // If several div tests have found no factors, check whether RemainingFactor is prime...
      if (p > 256 && CountDivTests == LogRemainingFactor && LogRemainingFactor < 64)
      {
        if (IsPrime(RemainingFactor)) break;
      }

      if (mpz_fdiv_ui(mpzref(RemainingFactor),p) != 0) continue;
      const long exp = DivideOutMaxPower_NoArgCheck(p, RemainingFactor); // updates RemainingFactor!!
      CoCoA_ASSERT(exp > 0);
      factors.push_back(BigInt(p));
      exponents.push_back(exp);
      LogRemainingFactor = FloorLog2(RemainingFactor);
      if (LogRemainingFactor < 64 && RemainingFactor/stop < stop)
        stop = ConvertTo<long>(FloorSqrt(RemainingFactor));
      CountDivTests = 0;
    }
    
    // if RemainingFactor is non-triv & below limit, add it to the list of factors found.
    if (RemainingFactor > 1 && RemainingFactor <= Pmax)
    {
      factors.push_back(RemainingFactor);
      exponents.push_back(1);
      RemainingFactor = 1;
    }
    if (N < 0)
      return factorization<BigInt>(factors, exponents, -RemainingFactor);
    return factorization<BigInt>(factors, exponents, RemainingFactor);
  }

  factorization<BigInt> factor_TrialDiv(const BigInt& N, const BigInt& TrialLimit)
  {
    if (IsZero(N))
      CoCoA_THROW_ERROR(ERR::ReqNonZero, "factor_TrialDiv(N,TrialLimit): 1st arg");
    if (TrialLimit < 1)
      CoCoA_THROW_ERROR(ERR::ReqPositive, "factor_TrialDiv(N,TrialLimit): 2nd arg");
    
    // Not implemented for large TrialLimit because it would be hideously slow...
    // A naive implementation could simply copy code from SmoothFactor(N,pmax) above.
    long pmax;
    if (!IsConvertible(pmax, TrialLimit))
      CoCoA_THROW_ERROR(ERR::NYI, "factor_TrialDiv(N,TrialLimit) with TrialLimit greater than largest signed long");
    return factor_TrialDiv(N, pmax);
  }


  factorization<long> factor(const MachineInt& n)
  {
    if (IsZero(n))
      CoCoA_THROW_ERROR(ERR::ReqNonZero, "factor(n)");
    if (!IsSignedLong(n))
      CoCoA_THROW_ERROR(ERR::ArgTooBig, "factor(n)");
    // Simple rather than efficient.
    if (uabs(n) < 2) return factor_TrialDiv(n,1);
    return factor_TrialDiv(n,uabs(n));
  }

  factorization<BigInt> factor(const BigInt& N)
  { return factor(N, CpuTimeLimit(20)); }  // arbitrarily limit to be at most about 20 secs

// //   factorization<BigInt> factor(const BigInt& N, const CpuTimeLimit& CheckTime)
// //   {
// //     // **NOTE** This fn should use StarRoot (when that fn is eventually impl'd)
// //     if (IsZero(N))
// //       CoCoA_THROW_ERROR(ERR::BadArg, "factor(N):  N must be non-zero");
// // //    const long PrimeLimit = FactorBigIntTrialLimit; // defined in config.H
// //     const long PrimeLimit = std::min(FactorBigIntTrialLimit, std::max(101L, 5*FloorLog2(abs(N))));  /// heuristic, why this formula???
// //     factorization<BigInt> ans = factor_TrialDiv(N, PrimeLimit);
// //     BigInt R = abs(ans.myRemainingFactor());
// // //std::clog<<"R="<<R<<std::endl;
// //     ans.myNewRemainingFactor(BigInt(sign(ans.myRemainingFactor())));
// //     if (R == 1) return ans;
// //     if (R < power(PrimeLimit,2) || IsProbPrime(R)) // ??? IsPrime or IsProbPrime ???
// //     {
// //       ans.myAppend(R,1);
// //       return ans;
// //     }
// // //DEBUG:   std::clog<<"Pollard loop:"<<std::endl;
// //     // ****IDEA**** Could check for abs(R) being a perfect power...  StarRoot???
// // //    const long PollardRhoIters = std::max(1000L, 1000*(100000000L/FloorLog2(R)));  // heuristic iter bound
// // //???    const long log2R = FloorLog2(R);
// // //???    const long HeuristicIterBlockSize = (log2R < 500)? 1000000 : static_cast<long>(2.0e+10/std::pow(log2R,1.6));  // no risk of overflow
// //     const long HeuristicIterBlockSize = static_cast<long>(2.0e+10/std::pow(FloorLog2(R),1.6));  // no risk of overflow
// // ///std::clog<<"log="<<FloorLog2(R)<<std::endl;
// // ///std::clog<<"BlockSize="<<HeuristicIterBlockSize<<std::endl;
// // const long PollardRhoIters = std::max(10L, std::min(1000000L, HeuristicIterBlockSize));
// //     CheckTime.myReset();
// //     while (true) // R is not prime, try to factorize it
// //     {
// // ///std::clog<<"J";
// //       if (CheckTime.IamTimedOut())
// //       { ans.myNewRemainingFactor(R*ans.myRemainingFactor()); return ans; } // Time's up!
// //       const factorization<BigInt> fac = factor_PollardRho(R, PollardRhoIters);
// //       if (fac.myFactors().empty())  continue;  // Pollard found no factor -- try more iters until time-out or success
// // /////      {ans.myNewRemainingFactor(R*ans.myRemainingFactor()); return ans;} // Pollard found no factor -- give up & return
// //       const BigInt& a = fac.myFactors().front(); 
// //       if (IsProbPrime(a))
// //       {
// //         const long m = FactorMultiplicity(a, R);
// //         ans.myAppend(a,m);
// //         R /= power(a,m);
// //         if (IsOne(R)) return ans;
// //         if (IsProbPrime(R)) { ans.myAppend(R,1); return ans; }
// //         continue;
// //       }
// // //DEBUG:   std::clog<<"SPLITTING COMPOSITE FACTOR"<<std::endl;
// // CoprimeFactorBasis_BigInt CFB;  CFB.myAddInfo(fac.myRemainingFactor()); CFB.myAddInfo(fac.myFactors()[0]);
// //       const vector<BigInt>& B = FactorBase(CFB);
// //       for (const BigInt& g: B)
// //       {
// // //DEBUG:   std::clog<<"g="<<g<<std::endl;
// //         const long m = FactorMultiplicity(g,R);
// //         R /= power(g,m);
// //         if (IsProbPrime(g))
// //         { ans.myAppend(g,m); continue; }
// // //DEBUG:   std::clog<<"g not prime: try recursing"<<std::endl;
// //         const factorization<BigInt> fac = factor(g, CheckTime);
// // //const factorization<BigInt> fac = factor_PollardRho(g, PollardRhoIters); // MUST DO CFB AGAIN!!!
// //         if (fac.myFactors().empty()) { /* composite factor*/ ans.myAppend(g,m); continue; }
// //         for (long i=0; i < len(fac.myFactors()); ++i)
// //           ans.myAppend(fac.myFactors()[i], m*fac.myMultiplicities()[i]);
// //       }
// //       return ans;
// // ///      CoCoA_THROW_ERROR(ERR::NYI, "factor(N): unimplemented -- PollardRho found composite factor");
// //     }
// // // NEVER GET HERE!!!
// // //    CoCoA_THROW_ERROR(ERR::NYI, "factor(N) unimplemented in this case -- too many large factors");
// // //    return ans;
// //   }


  namespace // anonymous
  {
struct PRS_Check
{
  PollardRhoSeq myPRS;
  int myCheck; // iter count when N inside PRS should be checked for primality
};
  } // end of namespace anonymous
  
  factorization<BigInt> factor(const BigInt& N, const CpuTimeLimit& CheckTime)
  {
    // **NOTE** BUG  BUG  This fn should use StarRoot (when that fn is eventually impl'd)
    if (IsZero(N))
      CoCoA_THROW_ERROR(ERR::ReqNonZero, "factor(N)");
//    const long PrimeLimit = FactorBigIntTrialLimit; // defined in config.H
    const long PrimeLimit = std::min(FactorBigIntTrialLimit, std::max(101L, 5*FloorLog2(abs(N))));  /// heuristic, why this formula???
////double t0=CpuTime();
    factorization<BigInt> ans = factor_TrialDiv(N, PrimeLimit);
////std::clog<<"Time for TrialDiv: "<<CpuTime()-t0<<std::endl;
    BigInt R = abs(ans.myRemainingFactor());
//std::clog<<"R="<<R<<std::endl;
    ans.myNewRemainingFactor(BigInt(sign(ans.myRemainingFactor())));
    if (R == 1) return ans;
    if (R < power(PrimeLimit,2) || IsProbPrime(R)) // ??? IsPrime or IsProbPrime ???
    {
      ans.myAppend(R,1);
      return ans;
    }
////std::clog<<"After IsProbPrime: "<<CpuTime()-t0<<std::endl;

///    long NumIters = 0;
//    vector<PRS_Check> L; L.push_back({PollardRhoSeq(R), 4});  // why 4???  Random magic number!?!
    vector<PollardRhoSeq> L;
    L.push_back(PollardRhoSeq(R));
    while (!L.empty() && !CheckTime.IamTimedOut())
    {
////      std::clog << "len(L)=" << len(L) << std::endl;
////double t1=CpuTime();
      for (auto& PRS: L)
        PRS.myAdvance(500000);  // MAGIC NUMBER -- BUG BUG BUG
////std::clog<<"Time for PRS: "<<CpuTime()-t1<<std::endl;
///      ++NumIters;
      
      vector<PollardRhoSeq> NewL;
      for (int i=0; i < len(L); ++i)
      {
        const BigInt& S = GetFactor(L[i]);
        if (IsOne(S))
        {
//???          if (L[i].myCheck == NumIters && IsProbPrime(GetN(L[i].myPRS)) { L[i].myCheck == 0; }
          NewL.push_back(L[i]);
          continue;
        }
////        std::clog << "FOUND A FACTOR: " << S << std::endl;
//        std::clog << "PRS is " << L[i] << std::endl;
        BigInt OtherFactor = GetN(L[i])/S;
//        std::clog << "OTHER FACTOR: " << OtherFactor << std::endl;
//        FactorsFound.push_back(S);
//        FactorsFound.push_back(OtherFactor);
        CoprimeFactorBasis_BigInt CFB;
        CFB.myAddInfo(vector<BigInt>{S, OtherFactor});
//        if (len(FactorBase(CFB))!=2) std::clog << "FactorBase is " << FactorBase(CFB) << std::endl;
        for (const BigInt& fac: FactorBase(CFB))
        {
          if (IsProbPrime(fac)) { long mult = FactorMultiplicity(fac,R); ans.myAppend(fac, mult); R /= power(fac,mult); continue; }
          NewL.push_back(PollardRhoSeq(L[i]).myReduceN(fac));
        }
      }
      swap(L, NewL);
    }

    if (!L.empty())
    {
      for (const auto& PRS: L)
        ans.myAppend(GetN(PRS), 1);
    }
      return ans;
  }


  //------------------------------------------------------------------

  PollardRhoSeq::PollardRhoSeq(const BigInt& N, long StartVal, long incr):
      myNumIters(1),
      myAnchorIterNum(1),
      myAnchorVal(StartVal),
      myN(abs(N)),
      myCurrVal(StartVal),
      myGCD(1),
      myStartVal(StartVal),
      myIncr(incr),
      myBlockSize(50)  // heuristic: 50-100 works well on my machine
  {
    if (myN < 2) CoCoA_THROW_ERROR(ERR::BadArg, "PollardRhoSeq ctor N must be >= 2");
    if (incr == 0) CoCoA_THROW_ERROR(ERR::ReqNonZero, "PollardRhoSeq ctor incr");
  }


  void PollardRhoSeq::myAdvance(long k)
  {
    CoCoA_ASSERT(k > 0);
    BigInt diff; // workspace used in loop below
    const long stop = k + myNumIters; // BUG  overflow???
    while (myNumIters < stop)
    {
      if (!IsOne(myGCD)) return;  // don't advance if myGCD != 1
      BigInt prod(1);
      for (int j=0; j < myBlockSize; ++j)
      {
        CheckForInterrupt("PollardRhoSeq::myAdvance");
        myCurrVal = (myCurrVal*myCurrVal+myIncr)%myN;
        if(myNumIters>=(3*myAnchorIterNum)/2) // check this!!!
        {
          diff = (myCurrVal-myAnchorVal);
          if (IsZero(diff)) // closed the cycle
          {
            myGCD = gcd(myN,prod);
            if (IsOne(myGCD)) { myNumIters=1; myAnchorIterNum=1;myStartVal=1;myAnchorVal=myStartVal; myCurrVal=myStartVal;++myIncr; myAdvance(k); return;}
//???myGCD = 0;
            return;
          }
        prod = (prod*diff)%myN;  // about 5% faster if reduce mod myN every 5 or 6 iters
        if (IsZero(prod)) { myGCD=gcd(diff,myN); return; }
        }
        ++myNumIters;
        if (myNumIters == 2*myAnchorIterNum)
        {
          myAnchorIterNum = myNumIters;
          myAnchorVal = myCurrVal;
        }
      }
      myGCD = gcd(myN, prod);
    }
  }


  PollardRhoSeq& PollardRhoSeq::myReduceN(const BigInt& fac)
  {
    CoCoA_ASSERT(myN%fac == 0);
    CoCoA_ASSERT(!IsProbPrime(fac));
    myAnchorVal %= fac;
    myN = fac;
    myCurrVal %= fac;
    myGCD = 1; //??? CORRECT BUG???
    return *this;
  }


  bool IsEnded(const PollardRhoSeq& PRS) noexcept
  {
    return !IsOne(PRS.myGCD);
  }


  const BigInt& GetN(const PollardRhoSeq& PRS) noexcept
  {
//    if (!IsEnded(PRS)) CoCoA_THROW_ERROR("","");
    return PRS.myN;
  }

  const BigInt& GetFactor(const PollardRhoSeq& PRS) noexcept
  {
//    if (!IsEnded(PRS)) CoCoA_THROW_ERROR("","");
    return PRS.myGCD;
  }


  long GetNumIters(const PollardRhoSeq& PRS) noexcept
  {
    return PRS.myNumIters;
  }


  std::ostream& operator<<(std::ostream& out, const PollardRhoSeq& PRS)
  {
    if (!out) return out;  // short-cut for bad ostreams
    if (out)
    {
      out << "PollardRhoSeq(N=" << PRS.myN
          << ",  start=" << PRS.myStartVal
          << ",  incr=" << PRS.myIncr
          << ",  NumIters=" << PRS.myNumIters
          << ",  gcd=" << PRS.myGCD << ")";
    }
    return out;
  }


  factorization<BigInt> factor_PollardRho(BigInt N, long NumIters)
  {
    factorization<BigInt> ans(BigInt(1));
    long a = 1; // pollard seq will use poly x^2+a; we may increase a
    do
    {
      PollardRhoSeq prs(N, a+1, a); // start from a+1 instead of 1 (to save 1 step)
      prs.myAdvance(NumIters);
      if (!IsOne(GetFactor(prs)) && !IsZero(GetFactor(prs)))
      {
        const long m = FactorMultiplicity(GetFactor(prs), N);
        ans.myAppend(GetFactor(prs),m);
        N /= power(GetFactor(prs),m);
        ans.myNewRemainingFactor(N);
        return ans;
      }
      NumIters -= GetNumIters(prs);
      ++a;
    } while (NumIters > 0);
    ans.myNewRemainingFactor(N);
    return ans;
  }


  //------------------------------------------------------------------

  // All positive factors of positive n; in incr order.
  // ??? Also worth having AllFactors_unsorted ???
  std::vector<long> AllFactors(long n)
  {
    if (n <= 0)
      CoCoA_THROW_ERROR(ERR::ReqPositive, "AllFactors");
    const auto facs = factor(n);
    const vector<long>& p = facs.myFactors();
    const vector<long>& mult = facs.myMultiplicities();
    long MaxSz = 1; for (long m: mult)  { MaxSz *= 1+m; }
    vector<long> AllFacs;  AllFacs.reserve(MaxSz);
    AllFacs.push_back(1);
    for (long k=0; k < len(p); ++k)
    {
      int stop = mult[k]*len(AllFacs);
      for (int i=0; i < stop; ++i)
        AllFacs.push_back(p[k]*AllFacs[i]);
    }
    std::sort(AllFacs.begin(), AllFacs.end());
    return AllFacs;
  }


  //------------------------------------------------------------------
  BigInt SumOfFactors(const MachineInt& n, long k)
  {
    if (IsNegative(n) || IsZero(n)) CoCoA_THROW_ERROR(ERR::ReqPositive, "SumOfFactors");
    if (IsNegative(k)) CoCoA_THROW_ERROR(ERR::ReqNonNegative, "SumOfFactors");
    const factorization<long> FacInfo = factor(n);
    const vector<long>& mult = FacInfo.myMultiplicities();
    const int s = len(mult);
    if (k == 0)
    {
      long ans = 1;
      for (int i=0; i < s; ++i)
        ans *= 1+mult[i];
      return BigInt(ans);
    }
    // Here k > 0
    BigInt ans(1);
    const vector<long>& fac = FacInfo.myFactors();
    for (int i=0; i < s; ++i)
      ans *= (power(fac[i], k*(mult[i]+1)) - 1)/(power(fac[i],k) - 1);
    return ans;
  }


  long SmallestNonDivisor(const MachineInt& N)
  {
    if (IsZero(N)) CoCoA_THROW_ERROR(ERR::ReqNonZero, "SmallestNonDivisor");
    unsigned long n = uabs(N);
    if (IsOdd(n)) return 2;
    FastFinitePrimeSeq TrialDivisorList;
    while (n%(*TrialDivisorList) == 0)
      ++TrialDivisorList;
    return *TrialDivisorList;
  }
  
  long SmallestNonDivisor(const BigInt& N)
  {
    if (IsZero(N)) CoCoA_THROW_ERROR(ERR::ReqNonZero, "SmallestNonDivisor");
    if (IsOdd(N)) return 2;
    // SLUG! simple rather than quick
    FastMostlyPrimeSeq TrialDivisorList;
    while (N%(*TrialDivisorList) == 0)
    {
      CheckForInterrupt("SmallestNonDivisor");
      ++TrialDivisorList;
    }
    return *TrialDivisorList;
  }


  bool IsSqFree(const MachineInt& N)
  {
    if (IsZero(N))
      CoCoA_THROW_ERROR(ERR::ReqNonZero, "IsSqFree(N):  N");
    if (!IsSignedLong(N))
      CoCoA_THROW_ERROR(ERR::ArgTooBig, "IsSqFree(N):  N must fit into long");
    long n = uabs(N); // implicit cast to long is safe
    // Handle primes 2,3,5 "manually":
    if (n < 4)  return true;
    if (n%4 == 0 || n%9 == 0 || n%25 == 0)  return false;
    n /= gcd(n,30);
    if (n < 49)  return true;
    
    // Main loop: we simply do trial divisions by the first few primes, then divisors from NoSmallFactorSeq
    long counter = 0;
    constexpr long TestIsPrime = 50; // Heuristic guess

    // Use "mostly primes"; much faster than using primes.
    FastMostlyPrimeSeq TrialFactorSeq;
    while (*TrialFactorSeq < 7) { ++TrialFactorSeq; } // since have already
    long p = *TrialFactorSeq;

    while (true)
    {
      ldiv_t qr = ldiv(n, p);
      if (p > qr.quot)  return true;  // test  equiv to p^2 > n
      if (qr.rem == 0)
      {
        n = qr.quot;  // equiv. to N /= p;
        if (n%p == 0)  return false;
      }
      else
      {
        if (qr.quot/p < p)  { const long sqrtn = FloorSqrt(n); return (n != sqrtn*sqrtn); }
        ++counter;
        if (counter == TestIsPrime && IsProbPrime(n))
          return true;
      }

      ++TrialFactorSeq;
      p = *TrialFactorSeq;
    }

    if (counter < TestIsPrime && IsProbPrime(n))  return true;
    return true;
  }


  bool3 IsSqFree(BigInt N)
  {
    if (IsZero(N))
      CoCoA_THROW_ERROR(ERR::ReqNonZero, "IsSqFree(N)");
    N = abs(N);
    if (N < 4) return true3;

    // Main loop: we simply do trial divisions by the first few primes, then divisors from NoSmallFactorSeq
    const long Pmax = FactorBigIntTrialLimit;
    long counter = 0;
    const long LogN = FloorLog2(N);
    const long TestIsPrime = (LogN>1000)?1000000:LogN*LogN;

    // Use "mostly primes"; much faster than using primes.
    FastMostlyPrimeSeq TrialFactorSeq;
    long p = *TrialFactorSeq;

    BigInt quot, rem;
    while (p <= Pmax)
    {
      quorem(quot, rem, N, p);
      if (p > quot) return true3; // test  equiv to p^2 > N
      if (IsZero(rem))
      {
        swap(N, quot);  // equiv. to N /= p;
        if (N%p == 0) return false3;
      }
      else
      {
        ++counter;
        if (counter == TestIsPrime && IsProbPrime(N))
          return true3;
      }

      ++TrialFactorSeq;
      p = *TrialFactorSeq;
    }

    if (IsPower(N)) return false3;
    if (counter < TestIsPrime && IsProbPrime(N)) return true3;
    return uncertain3;
  }


  long FactorMultiplicity(const MachineInt& b, const MachineInt& n)
  {
    const char* const FnName = "FactorMultiplicity";
    if (IsZero(n)) CoCoA_THROW_ERROR(ERR::ReqNonZero, FnName);
    if (IsNegative(b)) CoCoA_THROW_ERROR(ERR::ReqNonNegative, FnName);
    const long base = AsSignedLong(b);
    if (base == 0 || base == 1) CoCoA_THROW_ERROR(ERR::BadArg, FnName);
    long m = uabs(n);
    return DivideOutMaxPower_NoArgCheck(base, m);
  }


  long FactorMultiplicity(const MachineInt& b, BigInt N)
  {
    if (IsZero(N)) CoCoA_THROW_ERROR(ERR::ReqNonZero, "FactorMultiplicity");
    if (IsNegative(b)) CoCoA_THROW_ERROR(ERR::ReqNonNegative, "FactorMultiplicity");
    const unsigned long base = AsUnsignedLong(b);
    if (base == 0 || base == 1) CoCoA_THROW_ERROR(ERR::BadArg, "FactorMultiplicity");

    return DivideOutMaxPower_NoArgCheck(base, N); // can overwrite N
  }

  long FactorMultiplicity(const BigInt& B, BigInt N)
  {
    const char* const FnName = "FactorMultiplicity";
    if (IsZero(N)) CoCoA_THROW_ERROR(ERR::ReqNonZero, FnName);
    if (B < 2) CoCoA_THROW_ERROR(ERR::BadArg, FnName);
    return DivideOutMaxPower_NoArgCheck(B,N); // can overwrite N
  }


  // CoprimeFactor(N,b) = N/gcd(N,b^infty)
  long CoprimeFactor(long N, long b)
  {
    if (N == 0 || b == 0)
      CoCoA_THROW_ERROR("both args must be non-zero","CoprimeFactor");
    N = std::abs(N);
    long g = gcd(N,b);
    while (g != 1)
    {
      N /= g;
      g = gcd(g,N);
    }
    return N;
  }
  
  namespace // anonymous
  {
    BigInt CoprimeFactor2(BigInt N)
    {
      // N is non-zero
      N = abs(N);
      unsigned long shift = mpz_scan1(mpzref(N), 0);
      mpz_tdiv_q_2exp(mpzref(N), mpzref(N), shift);
      return N;
    }
  } // end of namespace anonymous

  BigInt CoprimeFactor(BigInt N, long b)
  {
    if (IsZero(N) || b == 0)
      CoCoA_THROW_ERROR("both args must be non-zero","CoprimeFactor");
    if (b == 2) return CoprimeFactor2(N);
    return CoprimeFactor(N, BigInt(b));
  }
  
  BigInt CoprimeFactor(BigInt N, const BigInt& b)
  {
    if (IsZero(N) || IsZero(b))
      CoCoA_THROW_ERROR("both args must be non-zero","CoprimeFactor");
    if (b == 2) return CoprimeFactor2(N);
// Is this commented out code worth keeping?
//    unsigned long logN = FloorLog2(N);
//    N = CoprimeFactor2(N);
//    unsigned long pwr2 = 0;
//    if (IsOdd(b)) pwr2 = logN - FloorLog2(N);
    N = abs(N);
    BigInt g = gcd(N,b);
    while (!IsOne(g))
    {
      N /= g;
      g = gcd(g*g,N);
    }
//    if (pwr2 > 0) mpz_mul_2exp(mpzref(N), mpzref(N), pwr2);
    return N;
  }



} // end of namespace CoCoA
