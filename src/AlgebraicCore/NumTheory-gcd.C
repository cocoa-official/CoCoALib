//   Copyright (c)  1999,2009-2011  John Abbott  &  Anna M. Bigatti

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


#include "CoCoA/NumTheory-gcd.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRatOps.H"
#include "CoCoA/assert.H"
#include "CoCoA/utils.H"

#include <limits>
using std::numeric_limits;
// #include <vector>
using std::vector;

namespace CoCoA
{

  // Non-negative gcd of machine integers.
  // Emphasis is on simplicity rather than utmost speed.
  long gcd(const MachineInt& a, const MachineInt& b) noexcept
  {
    unsigned long A = uabs(a);
    unsigned long B = uabs(b);
    // Dispose of the trivial cases first.
    if (A == 0) return B;
    if (B == 0) return A;
    if (A == 1 || B == 1) return 1;
    if (A == B) return A;

    // General case.
    if (A < B) std::swap(A, B);
    while (B != 0)
    {
      A %= B;
      std::swap(A, B);
    }
    return IntegerCast<long>(A);
  }

  BigInt gcd(const BigInt& A, const MachineInt& b)
  { if (IsZero(b)) return abs(A); else return BigInt(gcd(A%uabs(b), b)); }

  BigInt gcd(const MachineInt& a, const BigInt& B)
  { if (IsZero(a)) return abs(B); else return BigInt(gcd(a, B%uabs(a))); }

  BigInt gcd(const BigInt& A, const BigInt& B)
  {
    BigInt ans;
    mpz_gcd(mpzref(ans), mpzref(A), mpzref(B));
    return ans;
  }


  long gcd(const std::vector<long>& L) noexcept
  {
    const int n = len(L);
    long g = 0;
    for (int i=0; i < n; ++i)
    {
      g = gcd(g, L[i]);
      if (g == 1) return g;
    }
    return g;
  }
    
  BigInt gcd(const std::vector<BigInt>& L)
  {
    // simple rather than smart
    const int n = len(L);
    BigInt g;
    for (int i=0; i < n; ++i)
    {
      g = gcd(g, L[i]);
      if (IsOne(g)) return g;
    }
    return g;
  }

  // -------------------------------------------------------
  // IsCoprime (trivial defns)
  
  bool IsCoprime(const MachineInt& a, const MachineInt& b) noexcept
  { return (gcd(a,b) == 1); }
  
  bool IsCoprime(const BigInt& A,     const MachineInt& b) // noexcept (with better impl)
  { return IsOne(gcd(A,b)); }

  bool IsCoprime(const MachineInt& a, const BigInt& B)     // noexcept (with better impl)
  { return IsOne(gcd(a,B)); }
  
  bool IsCoprime(const BigInt& A,     const BigInt& B)
  { return IsOne(gcd(A,B)); }


  namespace // anonymous
  {
    
    // The fn below is OK -- it is not obvious, but the longs cannot overflow: for
    // a proof see Theorem 4.3 & the comment after it in Shoup's book
    // "A Computational Introduction to Number Theory and Algebra".

    // InvModNoArgCheck is faster than this; perhaps mimic InvModNoArgCheck, and then compute directly
    // the other cofactor at the end??? see redmine!
    long ExtendedEuclideanAlg(long& CofacA, long& CofacB, unsigned long a, unsigned long b) noexcept
    {
      if (a < b) return ExtendedEuclideanAlg(CofacB, CofacA, b, a);

      // Now we are sure that a >= b.
      long m00 = 1;
      long m01 = 0;
      long m10 = 0;
      long m11 = 1;

      while (true)
      {
        if (b==0) { CofacA = m00; CofacB = m01; return IntegerCast<long>(a); }
        long q = a/b; // can overflow (harmlessly) only if b == 1 upon 1st iteration
        a -= q*b;
        m00 -= q*m10;
        m01 -= q*m11;

        if (a == 0) { CofacA = m10; CofacB = m11; return IntegerCast<long>(b); }
        q = b/a;
        b -= q*a;
        m10 -= q*m00;
        m11 -= q*m01;
      }
    }

  } // end of namespace anonymous


  long ExtGcd(long& CofacA, long& CofacB, const MachineInt& a, const MachineInt& b)
  {
    if (IsZero(a) && IsZero(b)) CoCoA_THROW_ERROR(ERR::ReqNonZero, "ExtGcd (machine int)");
    const long g = ExtendedEuclideanAlg(CofacA, CofacB, uabs(a), uabs(b));
    if (IsNegative(a)) CofacA = -CofacA;
    if (IsNegative(b)) CofacB = -CofacB;
    return g;
  }

  BigInt ExtGcd(BigInt& CofacA, BigInt& CofacB, const BigInt& A, const BigInt& B)
  {
    if (IsZero(A) && IsZero(B)) CoCoA_THROW_ERROR(ERR::ReqNonZero, "ExtGcd (BigInt)");
    BigInt ans;
    mpz_gcdext(mpzref(ans), mpzref(CofacA), mpzref(CofacB), mpzref(A), mpzref(B));
    // GMP guarantees only that abs(CofacA) < abs(B), but we want 2*abs(CofacA) <= abs(B)/ans
    // so the next few lines check, and if necessary tweak, the cofactors.
    // The code below is needed only for GMP-4.3.1; GMP-4.2.4 gave correctly reduced answers.
    if (A == 0 || B == 0 || abs(A) == abs(B)) return ans;
    const BigInt CofacAModulus = abs(B)/ans;
    const BigInt CofacBModulus = abs(A)/ans;
    const BigInt q = RoundDiv(abs(CofacA), CofacAModulus);
    if (q != 0)
    {
      if (CofacA < 0)
      { CofacA += q*CofacAModulus; if (sign(A) == sign(B)) CofacB -= q*CofacBModulus; else CofacB += q*CofacBModulus; }
      else
      { CofacA -= q*CofacAModulus; if (sign(A) == sign(B)) CofacB += q*CofacBModulus; else CofacB -= q*CofacBModulus; }
    }
    if (2*abs(CofacB) <= CofacBModulus) return ans; // NB already have 2*abs(CofacA) <= CofacAModulus
    if (CofacB < 0)
    { CofacB += CofacBModulus; if (sign(A) == sign(B)) CofacA -= CofacAModulus; else CofacA += CofacAModulus; }
    else
    { CofacB -= CofacBModulus; if (sign(A) == sign(B)) CofacA += CofacAModulus; else CofacA -= CofacAModulus; }
    return ans;
  }


  GcdAndCofacs::GcdAndCofacs(const BigInt& g, const std::vector<BigInt>& cofacs):
      myGcd(g), myCofacs(cofacs)
  {}


  GcdAndCofacs ExtGcd(const std::vector<BigInt>& L)
  {
    if (L.empty()) return GcdAndCofacs(BigInt(0), vector<BigInt>());
    const int n = len(L);
    if (n == 1) return GcdAndCofacs(abs(L[0]), vector<BigInt>(1, BigInt(-sign(L[0]))));
    vector<BigInt> cofacs(n);
    const BigInt g = gcd(L);
    cofacs[0] = 1;
    BigInt GcdSoFar = L[0];
    BigInt cofac1, cofac2;
    for (int j=1; j < n; ++j)
    {
      GcdSoFar = ExtGcd(cofac1, cofac2, GcdSoFar, L[j]);
      if (!IsOne(cofac1))
        for (int i=0; i < j; ++i)
          cofacs[i] *= cofac1;
      cofacs[j] = cofac2;
      if (GcdSoFar == g) break;
    }
    return GcdAndCofacs(g, cofacs);
  }


  //-------------------------------------------------------

  long lcm(const MachineInt& a, const MachineInt& b)
  {
    if (IsZero(a) || IsZero(b)) return 0;
    const unsigned long AoverG = uabs(a)/gcd(a,b);
    if (uabs(b) > numeric_limits<long>::max()/AoverG)
      CoCoA_THROW_ERROR(ERR::ArgTooBig, "lcm(MachineInt,MachineInt)");
    return AoverG*uabs(b); // implicit cast to long is safe (see overflow check above)
  }

  BigInt lcm(const BigInt& A, const MachineInt& b)
  {
    return lcm(A, BigInt(b));
  }

  BigInt lcm(const MachineInt& a, const BigInt& B)
  {
    return lcm(BigInt(a), B);
  }

  BigInt lcm(const BigInt& a, const BigInt& b)
  {
    BigInt ans;
    mpz_lcm(mpzref(ans), mpzref(a), mpzref(b));
    return ans;
  }


  BigInt lcm(const std::vector<long>& L)
  {
    const int n = len(L);
    BigInt LCM(1);
    for (int i=0; i < n; ++i)
      LCM = lcm(LCM, L[i]);
    return LCM;
  }
  
  BigInt lcm(const std::vector<BigInt>& L)
  {
    // simple rather than smart
    const int n = len(L);
    BigInt LCM(1);
    for (int i=0; i < n; ++i)
      LCM = lcm(LCM, L[i]);
    return LCM;
  }


} // end of namespace CoCoA
