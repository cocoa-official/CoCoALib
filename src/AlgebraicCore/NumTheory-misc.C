//   Copyright (c)  1999,2009-2011,2024  John Abbott  &  Anna M. Bigatti
//   Some contributions from Nico Mexis


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


#include "CoCoA/NumTheory-misc.H"
#include "CoCoA/NumTheory-factor.H"
#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRatOps.H"
#include "CoCoA/error.H"
#include "CoCoA/interrupt.H"
#include "CoCoA/long32or64.H"


#include <cmath>
//using std::log;
#include <limits>
using std::numeric_limits;
//#include <vector>
using std::vector;

namespace CoCoA
{


  long EulerTotient(const MachineInt& n)
  {
    if (IsZero(n) || IsNegative(n))
      CoCoA_THROW_ERROR(ERR::ReqPositive, "EulerTotient(n)");
    if (!IsSignedLong(n)) CoCoA_THROW_ERROR(ERR::ArgTooBig, "EulerTotient(n)");
    const factorization<long> facpows = factor(n);
    const vector<long>& primes = facpows.myFactors(); // assume fully factorized
    long ans = AsSignedLong(n);
    for (const long p: primes)
    {
      ans = (ans/p)*(p-1); // exact division; cannot overflow
    }
    return ans;
  }

  BigInt EulerTotient(const BigInt& N)
  {
    if (N <= 0)
      CoCoA_THROW_ERROR(ERR::ReqPositive, "EulerTotient(N)");
    { long n; if (IsConvertible(n,N)) { return BigInt(EulerTotient(n)); } }
    const factorization<BigInt> facpows = factor(N);
    const vector<BigInt>& primes = facpows.myFactors();
    BigInt ans = N;
    for (const BigInt& p: primes)
    {
      if (!IsProbPrime(p))
        CoCoA_THROW_ERROR(ERR::ArgTooBig, "EulerTotient: unable to factorize arg");
      ans = (ans/p)*(p-1); // exact division
    }
    return ans;
  }


  //  true iff Eulertotient(n) == phi_n;  error if either arg is non-positive.
  bool IsEulerTotient(long n, long phi_n)
  {
    if (n < 1 || phi_n < 1)
      CoCoA_THROW_ERROR("Args must be positive", "IsEulerTotient");
    if (phi_n == 1) return (n == 1 || n == 2);
    if (IsOdd(phi_n)) return false;
    for (FastFinitePrimeSeq pseq; n > 1 && !IsEnded(pseq); ++pseq)
    {
      // if phi_n%4 == 2 then we can be cleverer (I think); means n must be prime power
      const long p = *pseq;
      if (n/p < p) { /*n is prime*/ return (phi_n == n-1); }
      if (n%p != 0) continue;
      const long m = FactorMultiplicity(p,n);
      // p_pwr = power(p,m-1)
      long p_pwr = 1; for (int i=1; i < m; ++i)  p_pwr *=p;  // cannot overflow!
      const long phi = (p-1) * p_pwr;   // cannot overflow!
      if (phi_n%phi != 0) return false;
      phi_n /= phi;
      n /= (p*p_pwr); // cannot overflow, exact division
      if (n != 1 && IsOdd(phi_n)) return false;
    }
    if (n == 1) return (phi_n == 1);
    CoCoA_THROW_ERROR(ERR::NYI, "IsEulerTotient: large factors");  // cannot happen???
  }



  namespace // anonymous
  {

    // Input plist is a list of primes such that (p-1) divides n
    // mode is either AllPreimages or SqFreePrimages (which allows a short-cut)
    // Algorithm is fairly obvious (but was apparently published somewhere)
    // We could make this faster by making a table of possible products of
    // the larger primes in plist; when recursion reaches the larger primes
    // we just run through the table -- is it worth it?
    void InvTotient_recursive(vector<long>& preimages, int i, const vector<long>& plist, long n, long fac, InvTotientMode mode)
    {
      if (n == 1)  { preimages.push_back(fac); return; }
      if (IsOdd(n))  return;
      if (i == len(plist))
      { if (n >= plist.back() && IsPrime(n+1))  preimages.push_back(fac*(n+1));  return; }
      const long p = plist[i];
      if (p-1 > n)  return;
//    /*??? worth it ???*/    if (p*p-p > n) { if (IsPrime(n+1)) preimages.push_back(fac*(n+1));  return; }
      InvTotient_recursive(preimages, i+1, plist, n, fac, mode);
      if (n%(p-1) != 0)  return;
      if (i&64)  CheckForInterrupt("InvTotient_recursive");  // Is there a better way???  This seems to cost about 10% time (amazing?!?)
      fac *= p;
      n /= (p-1);  // exact!
      while (true)
      {
        InvTotient_recursive(preimages, i+1, plist, n, fac, mode);
        if (mode == InvTotientMode::SqFreePreimages || n%p != 0)  return;
        n /= p;
        fac *= p;
      }
    }

  } // end of namespace anymous


  std::vector<long> InvTotient(long n, InvTotientMode mode /*=InvTotientMode::AllPreimages*/)
  {
    if (n < 1) CoCoA_THROW_ERROR("arg must be >= 1","InvTot");
   const unsigned long UPB = InvTotientBound_ulong(n);
   if (UPB == 1)
     CoCoA_THROW_ERROR(ERR::ArgTooBig, "InvTotient");
   if (UPB == 0) return vector<long>(); // no pre-images

    if (n == 1) return vector<long>{1,2};
    
    // plist will contain exactly those primes p such that p-1 divides n
    // No other prime can divide any preimage of n.
    vector<long> plist;
    const long sqrtn = static_cast<long>(std::floor(sqrt(n)));
    for (FastMostlyPrimeSeq PS; *PS <= 1+sqrtn; ++PS)
    {
      const long p = *PS;
      if (n%(p-1) == 0 && IsPrime(p))
        plist.push_back(p);
    }
    // Recursive fn puts values into preimages; not in incr order!
    vector<long> preimages;
    InvTotient_recursive(preimages, 0, plist, n, 1, mode);
    sort(preimages.begin(), preimages.end());   // ???sort???
    return preimages;
  }



  //--------------------------------------------
  // InvTotient (a bit faster than naive approach)
  // Sufficient for our current purposes.

  // Simple rather than efficient
  BigInt InvTotientBoundUpto(const BigInt& phi)
  {
    if (phi < 1)  CoCoA_THROW_ERROR("Arg must be > 0","InvTotientBoundUpto");
    if (IsOne(phi))  return BigInt(2);
    BigRat ratio(2,1);
    BigInt phi_reduced = phi;
    constexpr int pmax = 113; // PRIME!  Chosen because there is a big jump to the next prime
    const int p_next = NextPrime(pmax);
    for (int p=3; p <= pmax; p = NextPrime(p))
    {
      if (phi_reduced < p-1)
      {
        return (phi*num(ratio))/den(ratio); // integer division
      }
      ratio = ratio*BigRat(p,p-1);
      phi_reduced /= p-1;
    }
    const double log_excess = log(phi_reduced);
    const double log_base = std::log(p_next);
    const long pwr = static_cast<long>(std::floor(log_excess/log_base));
    ratio = ratio * power(BigRat(p_next, p_next-1), pwr);
    return (phi*num(ratio))/den(ratio); // integer division
  }


  namespace // anonymous
  {

    struct SmallRat
    {
      constexpr SmallRat(long n, long d): num(n), den(d) {}
      long num;
      long den;
    };

    long operator*(unsigned long n, SmallRat q)
    {
      // Assume result does not overflow -- caller checks this!
      return ConvertTo<unsigned long>((BigInt(n)*q.num)/q.den); // integer division!  Simple, but inefficient?
    }

    BigInt operator*(const BigInt& N, SmallRat q)
    {
      // Assume result does not overflow -- caller checks this!
      return (N*q.num)/q.den; // integer division!  Simple, but inefficient?
    }

    // BUG?  This table is probably far too big: maybe 10 entries suffice?
    // The gains are really not that great!
    // entry[k] is max{ n/phi(n) | phi(n)/2^k is odd }
    const SmallRat InvTotientRatioTable[] =
    { {0,1},
      {3,1},
      {7,2},
      {77,20},
      {35,8},
      {77,16},
      {1463,288},
      {3059,576},
      {19019,3456},
      {39767,6912},
      {1232777,207360},
      {7572773,1244160},
      {15474797,2488320},
      {913013023,144322560},
      {448769113,69672960},
      {913013023,139345920},
#ifndef CoCoA_32BIT_LONG
      {61171872541,9196830720},
      {4343202950411,643778150400},
      {7629074921,1114767360},
      {15521221391,2229534720},
      {1039921833197,147149291520}
#endif
      };
    constexpr long InvTotientRatioTableSize = (long)(sizeof(InvTotientRatioTable)/sizeof(SmallRat));

    // Power of 2 which divides n
    int PwrOf2(long n)
    {
      CoCoA_ASSERT(n > 0); // must have n != 0
      int pwr = 0;
      while (IsEven(n)) { n /= 2; ++pwr; }
      return pwr;
    }

    // Power of 2 which divides N; if greater than 20, return 21.
    int PwrOf2(BigInt N)
    {
      CoCoA_ASSERT(N > 0); // must have n != 0
      int pwr = 0;
      while (IsEven(N) && pwr <= 20) { N /= 2; ++pwr; }
      return pwr;
    }

  } // end of namespace anonymous


  /**
   * Provides a simple-to-compute bound for the inverse
   * of the Euler Totient function. All numbers k with
   * phi(k) = n will satisfy k < InversePhiBound(n).
   * See https://oeis.org/A355667
   */


  // Return B such that if EulerTotient(k) = n then k <= B
  // Returns 0 to mean no such k exists (not absolute check!)
  // Returns 1 if B would be too large for ulong!
  unsigned long InvTotientBound_ulong(const unsigned long n)
  {
    if (n == 1) return 2;
    if (IsOdd(n)) return 0;  // obv no pre-image
    const int pwr = PwrOf2(n);
    constexpr unsigned long UPB = numeric_limits<unsigned long>::max()/8;
    if (n <= UPB && pwr < InvTotientRatioTableSize)  return n*InvTotientRatioTable[pwr];
    return InvTotientBoundUpto_ulong(n);
  }

  
  // Return B such that if EulerTotient(k) = n then k <= B.
  // Returns 0 to mean no such k exists (not absolute check!).
  BigInt InvTotientBound(const BigInt& N)
  {
    constexpr unsigned long UPB = numeric_limits<unsigned long>::max()/8;
    if (N <= UPB)
      return BigInt(InvTotientBound_ulong(ConvertTo<unsigned long>(N)));
    const long pwr = PwrOf2(N);
    if (pwr < InvTotientRatioTableSize)
      return N*InvTotientRatioTable[pwr];
    return InvTotientBoundUpto(N); // fallback case.
  }


  // Return B such that B >= EulerTotient(k) for all k <= n
  // Returns 1 if B would be too large for ulong!
  unsigned long InvTotientBoundUpto_ulong(const unsigned long n)
  {
    if (n == 1)  return 2;
    const BigInt B = InvTotientBoundUpto(BigInt(n));
    unsigned long b;
    if (IsConvertible(b,B))  return b;
    return 1;  // to signify "would overflow"
  }

  
  // Return B such that B >= EulerTotient(k) for all k <= n
  BigInt InvTotientBoundUpto(const unsigned long n)
  {
    return InvTotientBoundUpto(BigInt(n));
  }


  //--------------------------------------------

  long MoebiusFn(const MachineInt &n)
  {
    if (IsZero(n) || IsNegative(n))
      CoCoA_THROW_ERROR(ERR::ReqPositive, "MoebiusFn(n)");
    if (!IsSignedLong(n))
      CoCoA_THROW_ERROR(ERR::ArgTooBig, "MoebiusFn(n):  n must fit into long");
    if (!IsSqFree(n))
      return 0;
    return SmallPower(-1, factor(n).myFactors().size());
  }






  //////////////////////////////////////////////////////////////////
  // Binomial repr of an integer.
  // A fair compromise between simplicity and speed.

  namespace  // anonymous -- file local fn
  {

    BigInt SearchUpwards(const BigInt& N, BigInt n, long r)
    {
      BigInt step(1);
      while (binomial(n+step,r) <= N)
      {
        n += step;
        step *= 2;
      }
      step /= 2; // step is power of 2 (may also be 1)
      while (step >= 1)
      {
        if (binomial(n+step, r) <= N)
          n += step;
        step /= 2;
      }
      return n;
    }

  } // end of anonymous namespace


  std::vector<BigInt> BinomialRepr(BigInt N, long r)
  {
    if (N < 0) CoCoA_THROW_ERROR(ERR::ReqNonNegative, "BinomialRepr(N,r): 1st arg");
    if (r < 1) CoCoA_THROW_ERROR(ERR::ReqPositive, "BinomialRepr(N,r): 2nd arg");
    vector<BigInt> ans(r+1);
    while (r > 0 && N > 0)
    {
      CheckForInterrupt("BinomialRepr");
      BigInt lwb = (r-1 + FloorRoot(power(2,r)*factorial(r)*N, r))/2;
      if (N <= r || lwb < r) lwb = r;
      if (lwb == r && N > r) lwb = r+1;
      ans[r] = SearchUpwards(N, lwb, r);
      N -= binomial(ans[r], r);
      --r;
    }
    return ans;
  }

  BigInt BinomialReprShift(const BigInt& N, long r, long shift1, long shift2)
  {
    const vector<BigInt> n = BinomialRepr(N, r);
    BigInt ans;
    for (long i=1; i <= r; ++i)
    {
      if (n[i] == 0) continue;
      if (i+shift2 < 0) continue;
      ans += binomial(n[i]+shift1, i+shift2);
    }
    return ans;
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/NumTheory-misc.C,v 1.21 2024/07/02 15:32:23 bigatti Exp $
// $Log: NumTheory-misc.C,v $
// Revision 1.21  2024/07/02 15:32:23  bigatti
// Summary: Changed error codes: LogZero into ReqNonZero
//
// and Not... into ReqNonNegative, ReqNonNegativeGrading, ReqPositive, ReqPositiveGrading
//
// Revision 1.20  2024/04/12 16:55:52  abbott
// Summary: Changed name of CPP symbol which signifies that longs are 32 bits (redmine 1661/1804)
//
// Revision 1.19  2024/03/25 19:08:51  abbott
// Summary: Changed to ifdef/ifndef
//
// Revision 1.18  2024/03/23 11:28:35  abbott
// Summary: Inserted missing include
//
// Revision 1.17  2024/03/23 11:28:07  abbott
// Summary: Revised to work on 32-bit platforms (redmine 1661)
//
// Revision 1.16  2024/03/08 19:59:59  abbott
// Summary: Used constexpr (redmine 1511)
//
// Revision 1.15  2024/01/25 08:41:31  abbott
// Summary: Removed lots of cruft
//
// Revision 1.14  2023/06/26 08:28:09  abbott
// Summary: Removed fprgotten debug print stmt
//
// Revision 1.13  2023/06/23 13:53:57  abbott
// Summary: Made InvTotient interruptible
//
// Revision 1.12  2023/05/11 19:26:15  abbott
// Summary: New impl and interface to InvTotient
//
// Revision 1.11  2023/03/15 21:31:50  abbott
// Summary: Added new fn InvTotientUpto
//
// Revision 1.10  2023/03/09 22:32:58  abbott
// Summary: Changed BinomialReprShift to take arg by const ref
//
// Revision 1.9  2023/03/06 22:08:44  abbott
// Summary: Added cast to signed long to avoid compiler warning
//
// Revision 1.8  2023/02/27 21:09:08  abbott
// Summary: Revised InvTotientBound (more complicated, but a bit better sometimes)
//
// Revision 1.7  2023/02/23 20:48:05  abbott
// Summary: Added InvTotientBoundUpto for BigInt; revised impls
//
// Revision 1.6  2023/02/16 20:46:33  abbott
// Summary: Many changes to InvTotient, InvTotientBound(Upto), InvTotientIter
//
// Revision 1.5  2023/01/31 12:29:15  abbott
// Summary: Added new fns InvTotient and InvTotientBound
//
// Revision 1.4  2022/11/30 15:15:13  abbott
// Summary: Added MoebiusFn (1st impl)
//
// Revision 1.3  2022/02/18 14:11:55  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.2  2020/06/17 15:49:24  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.1  2020/01/26 14:14:31  abbott
// Summary: Finished splitting NumTheory into smaller pieces (redming 1161)
//
// Revision 1.89  2019/09/16 17:35:15  abbott
// Summary: Changed name EulerPhi into EulerTotient; changed iroot into FloorRoot
//
// Revision 1.88  2019/03/18 11:24:58  abbott
// Summary: Split NumTheory into several smaller files
//
// Revision 1.87  2019/03/04 16:18:27  abbott
// Summary: Added KroneckerSymbol; improved SmoothFactor
//
// Revision 1.86  2018/05/22 14:16:39  abbott
// Summary: Split BigRat into BigRat (class defn + ctors) and BigRatOps
//
// Revision 1.85  2018/05/18 12:15:04  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.84  2018/05/05 15:23:46  abbott
// Summary: Added primality test in SmoothFactor if there were several tests without finding a factor, and the RemainingFactor is not too big
//
// Revision 1.83  2018/04/20 18:51:25  abbott
// Summary: Changed ctors for BigInt/BigRat from string or from MPZ/MPQ
//
// Revision 1.82  2018/04/18 14:27:57  abbott
// Summary: Removed dirty trick from SmoothFactor
//
// Revision 1.81  2018/03/13 16:15:06  abbott
// Summary: Removed ProbPrimeIters (it is already in NumTheory-prime.C)
//
// Revision 1.80  2018/03/13 14:52:33  abbott
// Summary: Cleaned impl of RatReconstructByContFrac; added (commented-out) code for verbose
//
// Revision 1.79  2018/02/27 17:30:22  abbott
// Summary: Renamed NumTheory_prime to NumTheory-prime; changed includes
//
// Revision 1.78  2018/02/27 10:50:56  abbott
// Summary: Moved stuff about primes into NumTheory_prime
//
// Revision 1.77  2018/02/19 10:17:33  abbott
// Summary: Added comment
//
// Revision 1.76  2018/02/15 17:26:11  abbott
// Summary: Added EratosthenesRange, and PrimeSeq
//
// Revision 1.75  2018/01/17 10:29:52  abbott
// Summary: Changed InvModNoCheck into InvModNoArgCheck
//
// Revision 1.74  2018/01/16 11:42:51  abbott
// Summary: Changed NoThrow into RtnZeroOnError
//
// Revision 1.73  2017/11/08 14:03:56  abbott
// Summary: Added new fn SumOfFactors
//
// Revision 1.72  2017/10/17 15:51:00  abbott
// Summary: Replaced gcd(...)==1 by IsCoprime(...)
//
// Revision 1.71  2017/10/17 15:44:26  abbott
// Summary: Added new fn IsCoprime
//
// Revision 1.70  2017/04/18 09:59:30  bigatti
// -- separated errors in myAddInfo
//
// Revision 1.69  2016/11/23 15:37:25  abbott
// Summary: Now SmallestNonDivisor returns 0 if NextPrime reaches itslimit
//
// Revision 1.68  2016/11/22 14:31:01  abbott
// Summary: Added SmallestNonDivisor
//
// Revision 1.67  2016/11/11 14:15:32  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.66  2016/10/25 20:54:09  abbott
// Summary: Added new fn IsSqFree (for BigInt and ringelem of PolyRing over field)
//
// Revision 1.65  2016/10/20 18:04:17  abbott
// Summary: Added radical; improved slightly factor and SmoothFactor
//
// Revision 1.64  2016/10/19 13:42:41  abbott
// Summary: Added radical for integers; fixed a bug in factor for MachineInt
//
// Revision 1.63  2016/07/04 11:05:39  abbott
// Summary: Removed two unnecessary assertions
//
// Revision 1.62  2016/06/29 13:11:37  abbott
// Summary: Added "NoThrow" option to InvMod; much cleaning inside NumTheory.C
//
// Revision 1.61  2016/04/19 08:53:26  bigatti
// -- fixed BinomialReprShift for 0 entries (case (45_2)^(+1) --> 55)
//
// Revision 1.60  2016/03/25 20:04:39  abbott
// Summary: Moved ExtendedEuclideanAlg into anon namespace
//
// Revision 1.59  2015/11/24 12:48:32  abbott
// Summary: Renamed "isqrt" --> "FloorSqrt"; also "ILogBase" --> "FloorLog"
//
// Revision 1.58  2015/11/21 19:17:53  abbott
// Summary: Added SimplestBinaryRatBetween; put InvModNoArgCheck_int into anon namespace
//
// Revision 1.57  2015/11/05 19:20:11  abbott
// Summary: Corrected two embarrassing bugs
//
// Revision 1.56  2015/11/05 14:18:19  abbott
// Summary: Changed InvMod to give error rather than return 0 if inv does not exist; added InvModNoArgCheck
//
// Revision 1.55  2015/10/09 18:26:36  abbott
// Summary: Corrected redmine reference
//
// Revision 1.54  2015/10/09 18:18:27  abbott
// Summary: Renamed "abs" to "uabs" for MachineInt; new fn "negate"; see redmine 783
//
// Revision 1.53  2015/07/30 15:57:43  abbott
// Summary: FactorMultiplicity now allows nonprime base
//
// Revision 1.52  2015/06/29 13:26:04  abbott
// Summary: Changed name "valuation" --> "FactorMultiplicity"; added new fns "BadMFactor"
// Author: JAA
//
// Revision 1.51  2014/10/28 15:12:08  abbott
// Summary: Renamed modulus --> CombinedModulus, residue --> CombinedResidue (for CRTMill)
// Author: JAA
//
// Revision 1.50  2014/09/16 10:41:41  abbott
// Summary: Added new fn eratosthenes (with doc, example, test)
// Author: JAA
//
// Revision 1.49  2014/09/01 16:25:52  abbott
// Summary: New condition for NextProbPrime/PrevProbPrime to use NextPrime/PrevPrime (previously allowed too large values); use new constant FactorBigIntTrialLimit; ExtGcd now gives error with args (0,0)
// Author: JAA
//
// Revision 1.48  2014/08/29 16:04:55  abbott
// Summary: Added optional 3rd arg to myAddInfo (so coprimality check is skipped)
// Author: JAA
//
// Revision 1.47  2014/05/06 13:13:36  abbott
// Summary: Removed useless fn PowerModLargeModulus
// Author: JAA
//
// Revision 1.46  2014/05/02 13:54:06  abbott
// Summary: Simplified ctor interface for RatReconstruct* (need explicit arg 0 for default behaviour)
// Author: JAA
//
// Revision 1.45  2014/04/24 09:54:50  abbott
// Summary: Corrected arg check in RatReconstructWithBounds
// Author: JAA
//
// Revision 1.44  2014/04/15 13:27:19  abbott
// Summary: Changed rtn type of PrimitiveRoot to long (for CoCoA-5/BuiltinOneLiners)
// Author: JAA
//
// Revision 1.43  2014/04/11 13:33:59  abbott
// Summary: Updated several assertions following revision of MaxSquarableInteger
// Author: JAA
//
// Revision 1.42  2014/04/04 10:15:17  abbott
// Summary: Updated to new interface for MaxSquarableInteger
// Author: JAA
//
// Revision 1.41  2014/03/24 12:09:21  abbott
// Summary: Major revision to public interface of factorization template class
// Author: JAA
//
// Revision 1.40  2014/01/16 16:09:54  abbott
// Added NumPartitions.
//
// Revision 1.39  2013/10/15 16:20:03  abbott
// Added valuation.
//
// Revision 1.38  2013/05/21 14:31:45  abbott
// Added BinomialRepr and BinomialReprShift to CoCoALib and CoCoA-5.
//
// Revision 1.37  2013/05/20 15:47:09  abbott
// Added new fn BinomialRepr (placed in NumTheory).
//
// Revision 1.36  2013/03/26 14:58:59  abbott
// Replaced calls to obsolete proc "convert" by calls to "ConvertTo<...>".
//
// Revision 1.35  2013/02/26 11:29:17  abbott
// Added impl of RatReconstructWithBounds
//
// Revision 1.34  2013/02/22 22:43:56  abbott
// Consequential change: new syntax for getting result from CRTMill.
//
// Revision 1.33  2013/02/22 18:56:50  abbott
// Added feature that RatReconstructByContFrac & RatReconstructByLattice
// ctors accept arg 0 to mean "use default value".
//
// Revision 1.32  2013/02/19 18:48:15  abbott
// Added printing for CRTMill and RatReconstructByContFrac and RatReconstructByLattice.
//
// Revision 1.31  2013/02/15 17:46:00  abbott
// Added RatReconstructByContFrac and RatReconstructByLattice.
//
// Revision 1.30  2012/12/12 18:25:06  abbott
// Added new fn IsFinal for ContFracIter.
// Corrected (subtle) bug in SimplestBigRatBetween.
//
// Revision 1.29  2012/12/12 10:39:06  abbott
// Corrected (embarrassing) typo in an assertion.
//
// Revision 1.28  2012/12/11 17:30:30  abbott
// Changed name from SimplestRationalInInterval to SimplestBigRatBetween.
// Also fixed a bug in the impl.
//
// Revision 1.27  2012/12/05 15:09:24  abbott
// Added new class ContFracApproximant.
// Added new fn SimplestBigRatBetween (NB name changed on 2012-12-11).
// Some minor cleaning.
//
// Revision 1.26  2012/12/04 20:14:11  abbott
// Added new class CRTMill.
// Improved impl of ContFracIter class.
// Fixed a bug in CFApproxIter class.
//
// Revision 1.25  2012/10/22 10:30:56  abbott
// SmoothFactor now allows limit to be 1 (previously it wanted limit >= 2).
//
// Revision 1.24  2012/10/05 09:30:35  abbott
// Changed myExponents into myMultiplicities.
//
// Revision 1.23  2012/09/26 12:51:26  abbott
// Clarified a comment about IsProbPrime.
//
// Revision 1.22  2012/06/29 15:17:27  abbott
// Fixed issue #199.
//
// Revision 1.21  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.20  2012/03/16 15:40:12  abbott
// Merged contents of NumTheoryQQ (continued fraction functions) into NumTheory.
// Merged the doc too.
//
// Revision 1.19  2011/11/09 14:09:53  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.18  2011/09/06 13:39:08  abbott
// Minor cleaning: now uses "int" for increments and decrements in NextPrime & PrevPrime.
// Fixed (embarrassing) overflow bug.
//
// Revision 1.17  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.16  2011/03/23 21:00:46  abbott
// Removed FindPrimRoot from NumTheory.H because it was already
// available as PrimitiveRoot (a better name).
// Updated documentation for NumTheory.
//
// Revision 1.15  2011/03/22 20:17:18  abbott
// Added fn FindPrimRoot.
// Merged impls from obsolescent SmallPrime.C.
//
// Revision 1.14  2011/03/16 13:26:36  abbott
// Removed all "unsigned" from fn interfaces, and many unsigned from inside fn impls.
//
// Revision 1.13  2011/01/19 16:12:18  bigatti
// -- added ERR::ModulusLT2
//
// Revision 1.12  2010/03/05 21:37:58  abbott
// Completed implementation of MultiplicativeOrderMod.
//
// Revision 1.11  2010/03/03 15:02:09  abbott
// Corrected bug in PowerMod: forgot to give error for negative power of non-invertible base.
// Added tests of PowerMod to test-NumTheory1.C (and consequent change to expected output).
//
// Revision 1.10  2010/03/03 10:43:34  abbott
// Added PrimitiveRoot for big primes (might be very slow).
// Added MultiplicativeOrderMod (currently very SLUGGY implementation).
//
// Revision 1.9  2009/12/29 22:44:32  abbott
// Removed buggy proxy class ZZ::rtn.
// Consequent changes for function prototypes also in NumTheory.
// Corrected some minor buglets in NumTheory.
//
// Revision 1.8  2009/12/11 11:46:32  abbott
// Changed fn  convert  into  IsConvertible.
// Added template procedure  convert.
// New version because change is not backward compatible.
//
// Revision 1.7  2009/12/03 17:41:47  abbott
// Minor correction to a comment.
//
// Revision 1.6  2009/10/08 13:39:47  abbott
// Renamed "round" into "RoundDiv".
// Added some new versions of "RoundDiv".
// Added a test for "RoundDiv".
//
// Revision 1.5  2009/09/29 12:26:13  abbott
// Added patch/workaround for a problem in gmp-4.3.1, so that exgcd gives
// the right answers (gmp-4.3.1 no longer produces reduced cofactors).
//
// Revision 1.4  2009/07/06 12:30:48  abbott
// Fixed incorrect assertion in IsBigPrime.
// Corrected a bug in SmoothFactor (had > instead of >=), and tidied up a little.
//
// Revision 1.3  2009/07/02 16:28:10  abbott
// Fairly extensive change to NumTheory (at least internally and philosophically).
// Improved and cleaned NumTheory.  Added documentation.
// Clarified the exact behaviour of most functions.
//
// Revision 1.2  2009/06/11 14:10:58  abbott
// Added commented out procedural forms for gcd/lcm, in case we should
// later want to activate them.
//
// Revision 1.1  2009/06/05 12:14:55  abbott
// Major change:
//   created new files NumTheory.H/C  which contain basic number theory operations
//   removed several basic number theory operations from ZZ.H/C
//   removed gcd from MachineInt.H/C
//   changed names of some basic fns:
//      IsPPrime -> IsProbPrime
//      invmod -> InvMod    (changed signature too)
//      powermod -> PowerMod  (changed signature too)
//   added new fns
//      NextProbPrime & PrevProbPrime
//   consequent changes to other code and tests and examples
//
//
