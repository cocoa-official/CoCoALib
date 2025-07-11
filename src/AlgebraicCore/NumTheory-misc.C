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
      CoCoA_THROW_ERROR1(ERR::ReqPositive);
    if (!IsSignedLong(n))  CoCoA_THROW_ERROR1(ERR::ArgTooBig);
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
      CoCoA_THROW_ERROR1(ERR::ReqPositive);
    { long n; if (IsConvertible(n,N)) { return BigInt(EulerTotient(n)); } }
    const factorization<BigInt> facpows = factor(N);
    const vector<BigInt>& primes = facpows.myFactors();
    BigInt ans = N;
    for (const BigInt& p: primes)
    {
      if (!IsProbPrime(p))
        CoCoA_THROW_ERROR2(ERR::ArgTooBig, "unable to factorize arg");
      ans = (ans/p)*(p-1); // exact division
    }
    return ans;
  }


  //  true iff Eulertotient(n) == phi_n;  error if either arg is non-positive.
  bool IsEulerTotient(long n, long phi_n)
  {
    if (n < 1 || phi_n < 1)
      CoCoA_THROW_ERROR1(ERR::ReqPositive);
    if (phi_n == 1)  return (n == 1 || n == 2);
    if (IsOdd(phi_n))  return false;
    for (FastFinitePrimeSeq pseq; n > 1 && !IsEnded(pseq); ++pseq)
    {
      // if phi_n%4 == 2 then we can be cleverer (I think); means n must be prime power
      const long p = *pseq;
      if (n/p < p)  { /*n is prime*/ return (phi_n == n-1); }
      if (n%p != 0)  continue;
      const long m = FactorMultiplicity(p,n);
      // p_pwr = power(p,m-1)
      long p_pwr = 1; for (int i=1; i < m; ++i)  p_pwr *=p;  // cannot overflow!
      const long phi = (p-1) * p_pwr;   // cannot overflow!
      if (phi_n%phi != 0)  return false;
      phi_n /= phi;
      n /= (p*p_pwr); // cannot overflow, exact division
      if (n != 1 && IsOdd(phi_n))  return false;
    }
    if (n == 1)  return (phi_n == 1);
    CoCoA_THROW_ERROR2(ERR::NYI, "for large factors");  // cannot happen???
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
    if (n < 1)  CoCoA_THROW_ERROR1(ERR::ReqPositive);
   const unsigned long UPB = InvTotientBound_ulong(n);
   if (UPB == 1)
     CoCoA_THROW_ERROR1(ERR::ArgTooBig);
   if (UPB == 0)  return vector<long>(); // no pre-images

    if (n == 1)  return vector<long>{1,2};
    
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
    if (phi < 1)  CoCoA_THROW_ERROR1(ERR::ReqPositive);
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
      CoCoA_THROW_ERROR1(ERR::ReqPositive);
    if (!IsSignedLong(n))
      CoCoA_THROW_ERROR2(ERR::ArgTooBig, "arg must fit into long");
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
    if (N < 0)  CoCoA_THROW_ERROR2(ERR::ReqNonNegative, "1st arg");
    if (r < 1)  CoCoA_THROW_ERROR2(ERR::ReqPositive, "2nd arg");
    vector<BigInt> ans(r+1);
    while (r > 0 && N > 0)
    {
      CheckForInterrupt("BinomialRepr");
      BigInt lwb = (r-1 + FloorRoot(power(2,r)*factorial(r)*N, r))/2;
      if (N <= r || lwb < r)  lwb = r;
      if (lwb == r && N > r)  lwb = r+1;
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
