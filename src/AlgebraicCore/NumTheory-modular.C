//   Copyright (c)  1999,2009-2011  John Abbott,  Anna M. Bigatti

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


#include "CoCoA/NumTheory-modular.H"
#include "CoCoA/NumTheory-gcd.H"
#include "CoCoA/NumTheory-factor.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/assert.H"
#include "CoCoA/error.H"
#include "CoCoA/utils.H"


#include <limits>
using std::numeric_limits;
#include <vector>
using std::vector;

namespace CoCoA
{

  // Anonymous namespace for file local "static" variables.
  namespace
  {

    // Function to compute base^exp mod modulus; assumes 0 <= base < modulus.
    // Unclever for exp >= modulus.  modulus need not be prime, but (modulus-1)^2
    // must fit inside an unsigned long.
    // This is an iterative implementation of binary powering; seems to be
    // usefully faster than the more obvious recursive implementation.
    // If modulus=1, always produces 0.  Produces 1 when asked to compute 0^0.
    unsigned long PowerModSmallModulus(unsigned long base, unsigned long exp, unsigned long modulus)
    {
      CoCoA_ASSERT(modulus != 0);
      CoCoA_ASSERT(modulus > 1);
      CoCoA_ASSERT(modulus-1 <= MaxSquarableInteger<unsigned long>());

      if (exp == 0 || base == 1) return 1;
      if (exp == 1 || base == 0) return base;
    
      unsigned long ans = 1;
      while (exp > 1)
      {
        if (exp&1) ans = (ans*base)%modulus;
        exp /= 2;
        base = (base*base)%modulus;
      }
      return (ans*base)%modulus;
    }


    // Common sanity check for the modulus
    void CheckModulus(const MachineInt& modulus, const char* const FnName)
    {
      if (IsNegative(modulus) || !IsSignedLong(modulus) || AsSignedLong(modulus) < 2)
        CoCoA_THROW_ERROR(ERR::BadModulus, FnName);
    }

    // Common sanity check for the modulus
    void CheckModulus(const BigInt& modulus, const char* const FnName)
    {
      if (modulus < 2)
        CoCoA_THROW_ERROR(ERR::BadModulus, FnName);
    }

  } // end of anonymous namespace


  // Compute (base^exp)%modulus  with exp >= 0  and modulus >= 2
  // The case 0^0 yields 1.
  long PowerMod(const MachineInt& base, const MachineInt& exponent, const MachineInt& modulus)
  {
    CheckModulus(modulus, "PowerMod(base,exponent,modulus)");
    if (IsNegative(exponent))
    {
      const long inv = InvMod(base, modulus); // Throws ERR::DivByZero if inverse does not exist.
      return PowerMod(inv, uabs(exponent), modulus);
    }
    const unsigned long e = AsUnsignedLong(exponent);
    if (e == 0) return 1; // Note that 0^0 gives 1
    const unsigned long m = AsUnsignedLong(modulus);
    unsigned long b = uabs(base)%m;
    if (b == 0) return 0;
    if (IsNegative(base) && (e%2 == 1))
      b = m-b;
    if (m <= MaxSquarableInteger<unsigned long>())
      return PowerModSmallModulus(b,e,m); // Call the fn in the anon namespace above.
    // Square of modulus doesn't fit into unsigned long, so compute with BigInts and convert answer back to unsigned long.
    return ConvertTo<long>(PowerMod(b,e,BigInt(m)));
  }

  long PowerMod(const MachineInt& base, const BigInt& exponent, const MachineInt& modulus)
  {
    // BUG/SLUG: horribly inefficient!!
    // What is the best way to implement this fn?????  Any volunteers?
    CheckModulus(modulus, "PowerMod(base,exponent,modulus)");
    return ConvertTo<long>(PowerMod(BigInt(base), exponent, BigInt(modulus)));
  }

  long PowerMod(const BigInt& base, const MachineInt& exponent, const MachineInt& modulus)
  {
    CheckModulus(modulus, "PowerMod(base,exponent,modulus)");
    return PowerMod(base%modulus, exponent, modulus);
  }

  long PowerMod(const BigInt& base, const BigInt& exponent, const MachineInt& modulus)
  {
    // BUG/SLUG: horribly inefficient!!
    // What is the best way to implement this fn?????  Any volunteers?
    CheckModulus(modulus, "PowerMod(base,exponent,modulus)");
    return ConvertTo<long>(PowerMod(base, exponent, BigInt(modulus)));
  }

  BigInt PowerMod(const MachineInt& base, const MachineInt& exponent, const BigInt& modulus)
  {
    return PowerMod(BigInt(base), BigInt(exponent), modulus);
  }

  BigInt PowerMod(const BigInt& base, const MachineInt& exponent, const BigInt& modulus)
  {
    CheckModulus(modulus, "PowerMod(BigInt, MachineInt, BigInt)");
    if (IsNegative(exponent))
    {
      const BigInt inv = InvMod(base, modulus); // Throws ERR::DivByZero if inverse does not exist.
      return PowerMod(inv, uabs(exponent), modulus);
    }
    BigInt ans;
    // ASSUME: line below gives 1 for 0^0
    mpz_powm_ui(mpzref(ans), mpzref(base), AsUnsignedLong(exponent), mpzref(modulus));
    return ans;
  }


  BigInt PowerMod(const MachineInt& base, const BigInt& exponent, const BigInt& modulus)
  {
    CheckModulus(modulus, "PowerMod(MachineInt, BigInt, BigInt)");
    if (exponent < 0)
    {
      const BigInt inv = InvMod(base, modulus); // Throws ERR::DivByZero if inverse does not exist.
      return PowerMod(inv, -exponent, modulus);
    }
    return PowerMod(BigInt(base), exponent, modulus);
  }


  BigInt PowerMod(const BigInt& base, const BigInt& exponent, const BigInt& modulus)
  {
    CheckModulus(modulus, "PowerMod(BigInt, BigInt, BigInt)");
    if (IsZero(exponent)) return BigInt(1); // x^0 --> 1 even for x==0
    if (exponent < 0)
    {
      const BigInt inv = InvMod(base, modulus); // Throws ERR::DivByZero if inverse does not exist.
      return PowerMod(inv, -exponent, modulus);
    }
    BigInt ans;
    mpz_powm(mpzref(ans), mpzref(base), mpzref(exponent), mpzref(modulus));
    return ans;
  }



  namespace // anonymous for file local fn
  {
    // On some platforms this "unsigned int" version is significantly faster
    // than the standard "unsigned long" version.  This version is "hidden
    // inside" the unsigned long version -- it calls this version if it can.

    // JAA believes that this routine always returns a reduced value,
    // but does not have a proof of this currently.
    unsigned int InvModNoArgCheck_int(unsigned int r, unsigned int m, const InvModErrorAction ErrorAction = ThrowOnError)
    {
//      CoCoA_ASSERT(r >= 0);
      CoCoA_ASSERT(r < m);
      CoCoA_ASSERT(m >= 2);
      const unsigned int p = m;
      unsigned int cr = 1;
      unsigned int cm = 0; // this is minus what you think it is!
      while (r != 0)
      {
        {
          const unsigned int q = m/r;
          m -= q*r;
          cm += q*cr;
        }
        if (m == 0) break;
        {
          const unsigned int q = r/m;
          r -= q*m;
          cr += q*cm;
        }
      }
      if (r+m != 1)
      {
        if (ErrorAction == ThrowOnError) CoCoA_THROW_ERROR(ERR::DivByZero, "InvMod");
        return 0;
      }
      if (r == 0) return p-cm;
      return cr;
    }

  } // end of anonymous namespace

  // JAA believes that this routine always returns a reduced value,
  // but does not have a proof of this currently.
  // If r has no inverse mod m then return 0 (if ErrorAction == RtnZeroOnError),
  // o/w throw ERR::DivByZero (if ErrorAction == ThrowOnError).
  unsigned long InvModNoArgCheck(unsigned long r, unsigned long m, const InvModErrorAction ErrorAction /*= ThrowOnError*/)
  {
//    CoCoA_ASSERT(r >= 0);
    CoCoA_ASSERT(r < m);
    CoCoA_ASSERT(m >= 2);
    CoCoA_ASSERT(m <= static_cast<unsigned long>(numeric_limits<long>::max())); // check fits into SIGNED long
    // special case if modulus fits into "unsigned int"
    if (m <= numeric_limits<unsigned int>::max())
      return InvModNoArgCheck_int(r, m);
    const unsigned long p = m;
    unsigned long cr = 1;
    unsigned long cm = 0; // this is minus what you think it is!
    while (r != 0)
    {
      {
        const unsigned long q = m/r;
        m -= q*r;
        cm += q*cr;
      }
      if (m == 0) break;
      {
        const unsigned long q = r/m;
        r -= q*m;
        cr += q*cm;
      }
    }
    if (r+m != 1)
    {
      if (ErrorAction == ThrowOnError) CoCoA_THROW_ERROR(ERR::DivByZero, "InvMod");
      return 0;
    }
    if (r == 0) return p-cm;
    return cr;
  }

  
  long InvMod(const MachineInt& r, const MachineInt& m, const InvModErrorAction ErrorAction /*= ThrowOnError*/)
  {
    CheckModulus(m, "InvMod(residue,modulus)");
    const unsigned long modulus = uabs(m);
    unsigned long residue = LeastNNegRemainder(r,m);
    return InvModNoArgCheck(residue, modulus, ErrorAction);
  }

  long InvMod(const BigInt& r, const MachineInt& m, const InvModErrorAction ErrorAction /*= ThrowOnError*/)
  {
    CheckModulus(m, "InvMod(residue,modulus)"); // throws if m < 2
    const long residue = LeastNNegRemainder(r,m);
    return InvModNoArgCheck(residue, AsUnsignedLong(m), ErrorAction);
  }

  BigInt InvMod(const MachineInt& r, const BigInt& m, const InvModErrorAction ErrorAction /*= ThrowOnError*/)
  {
    CheckModulus(m, "InvMod(residue,modulus)"); // throws if m < 2
    const BigInt residue(r);
    BigInt ans;
    const int InverseExists = mpz_invert(mpzref(ans), mpzref(residue), mpzref(m));
    if (InverseExists) return ans;
    // Not invertible, so handle error.
    if (ErrorAction == ThrowOnError) CoCoA_THROW_ERROR(ERR::DivByZero, "InvMod");
    return BigInt(0);
  }

  BigInt InvMod(const BigInt& r, const BigInt& m, const InvModErrorAction ErrorAction /*= ThrowOnError*/)
  {
    CheckModulus(m, "InvMod(residue,modulus)"); // throws if m < 2
    BigInt ans;
    const int InverseExists = mpz_invert(mpzref(ans), mpzref(r), mpzref(m));
    if (InverseExists) return ans;
    // Not invertible, so handle error.
    if (ErrorAction == ThrowOnError) CoCoA_THROW_ERROR(ERR::DivByZero, "InvMod");
    return BigInt(0);
  }



  //------------------------------------------------------------------

  // Taken from Cohen's book (page 24).
  long MultiplicativeOrderModPrime(unsigned long residue, unsigned long p)
  {
    CoCoA_ASSERT(IsPrime(p));
    CoCoA_ASSERT(p-1 <= MaxSquarableInteger<unsigned long>());
    CoCoA_ASSERT(IsCoprime(residue, p));

    if (p ==  2) return 1;
    const factorization<long> facpows = factor(p-1);
    const vector<long>& q = facpows.myFactors();
    const vector<long>& pwr = facpows.myMultiplicities();
    const int n = len(q);
    long e = p-1;
    for (int i=0; i < n; ++i)
    {
      e /= SmallPower(q[i], pwr[i]);
      unsigned long rpower = PowerMod(residue, e, p);
      while (rpower != 1)
      {
        rpower = PowerMod(rpower, q[i], p);
        e *= q[i];
      }
    }
    return e;
  }

  long MultiplicativeOrderModPrimePower(unsigned long residue, unsigned long p, int e)
  {
    CoCoA_ASSERT(e > 0);
    CoCoA_ASSERT(residue >= 1 && residue < p);
    CoCoA_ASSERT(IsPrime(p));
    CoCoA_ASSERT(p-1 <= MaxSquarableInteger<unsigned long>());
    long ord =  MultiplicativeOrderModPrime(residue, p);
    if (e == 1) return ord;
    const unsigned long q = SmallPower(p,e);
    unsigned long rpower = PowerMod(residue, ord, q);
    while (rpower != 1)
    {
      rpower = PowerMod(rpower, p, q);
      ord *= p;
    }
    return ord;
  }


  // Taken from Cohen's book (page 24).
  BigInt MultiplicativeOrderModPrime(const BigInt& residue, const BigInt& p)
  {
    CoCoA_ASSERT(IsPrime(p));
    CoCoA_ASSERT(IsCoprime(residue, p));

// ???BUG NYI: use machine int code if possible???
    if (p == 2) return BigInt(1);
    const factorization<BigInt> facpows = factor(p-1);
    const vector<BigInt>& q = facpows.myFactors();
    const vector<long>& pwr = facpows.myMultiplicities();
    const int n = len(q);
    BigInt e = p-1;
    for (int i=0; i < n; ++i)
    {
      e /= power(q[i], pwr[i]);
      BigInt rpower = PowerMod(residue, e, p);
      while (rpower != 1)
      {
        rpower = PowerMod(rpower, q[i], p);
        e *= q[i];
      }
    }
    return e;
  }

  BigInt MultiplicativeOrderModPrimePower(const BigInt& residue, const BigInt& p, long e)
  {
    CoCoA_ASSERT(e > 0);
    CoCoA_ASSERT(residue >= 1 && residue < p);
    CoCoA_ASSERT(IsPrime(p));
    BigInt ord =  MultiplicativeOrderModPrime(residue, p);
    if (e == 1) return ord;
    const BigInt q = power(p,e);
    BigInt rpower = PowerMod(residue, ord, q);
    while (rpower != 1)
    {
      rpower = PowerMod(rpower, p, q);
      ord *= p;
    }
    return ord;
  }


  long MultiplicativeOrderMod(const MachineInt& residue, const MachineInt& modulus)
  {
    if (IsNegative(modulus) || AsUnsignedLong(modulus) < 2)
      CoCoA_THROW_ERROR(ERR::BadArg, "MultiplicativeOrderMod: must have modulus >= 2");
    if (gcd(residue, modulus) != 1)
      CoCoA_THROW_ERROR(ERR::BadArg, "MultiplicativeOrderMod: residue must be coprime to modulus");
    const unsigned long m = AsUnsignedLong(modulus);
    if (m-1 > MaxSquarableInteger<unsigned long>())
    {
      // return ConvertTo<long>(MultiplicativeOrderMod(residue, BigInt(modulus)));
      return ConvertTo<long>(MultiplicativeOrderMod(residue, BigInt(modulus)));
    }

    unsigned long r = uabs(residue)%m;
    if (r != 0 && IsNegative(residue)) r = m-r;

    const factorization<long> mfacs = factor(m);
    const vector<long>& factor = mfacs.myFactors();
    const vector<long>& exponent = mfacs.myMultiplicities();
    const int n = len(factor);
    long ord = 1;
    for (int i=0; i < n; ++i)
    {
      unsigned long facpow = SmallPower(factor[i], exponent[i]);
      ord = lcm(ord, MultiplicativeOrderModPrimePower(r%facpow, factor[i], exponent[i]));
    }
    return ord;
  }

  long MultiplicativeOrderMod(const BigInt& residue, const MachineInt& modulus)
  {
    return MultiplicativeOrderMod(residue%uabs(modulus), modulus);
  }

  BigInt MultiplicativeOrderMod(const MachineInt& residue, const BigInt& modulus)
  {
    return MultiplicativeOrderMod(BigInt(residue), modulus);
  }

  BigInt MultiplicativeOrderMod(const BigInt& residue, const BigInt& modulus)
  {
    if (modulus < 2)
      CoCoA_THROW_ERROR(ERR::BadArg, "MultiplicativeOrderMod: must have modulus >= 2");
    if (gcd(residue, modulus) != 1)
      CoCoA_THROW_ERROR(ERR::BadArg, "MultiplicativeOrderMod: residue must be coprime to modulus");
    unsigned long m;
    if (IsConvertible(m, modulus) && m-1 <= MaxSquarableInteger<unsigned long>())
    {
      return BigInt(MultiplicativeOrderMod(residue%m, m));
    }

    const factorization<BigInt> mfacs = factor(modulus);
    const vector<BigInt>& factor = mfacs.myFactors();
    const vector<long>& exponent = mfacs.myMultiplicities();
    const int n = len(factor);
    BigInt ord(1);
    for (int i=0; i < n; ++i)
    {
      const BigInt facpow = power(factor[i], exponent[i]);
      ord = lcm(ord, MultiplicativeOrderModPrimePower(residue%facpow, factor[i], exponent[i]));
    }
    return ord;
  }

  // BigInt version of MultiplicativeOrder???  Could be VERY SLOW!!!

  long PrimitiveRoot(const MachineInt& pp)
  {
    if (IsNegative(pp) || !IsSignedLong(pp) || !IsPrime(pp))
      CoCoA_THROW_ERROR(ERR::BadArg, "PrimitiveRoot(p):  p must be a (positive) prime");

    unsigned long p = AsUnsignedLong(pp);
    if (p == 2) return 1;
    const factorization<long> facpows = factor(p-1);
    const vector<long>& primes = facpows.myFactors();
    const int NumPrimes = len(primes);
    for (unsigned long root = 2; /*empty*/; ++root)
    {
      CoCoA_ASSERT(root < p);
      if (root == 4 || root == 8 || root == 9 || root == 16) continue; // skip the first few prime powers
      bool failed = false;
      for (int i=0; i < NumPrimes; ++i)
      {
	if (PowerMod(root, (p-1)/primes[i], p) == 1)
        { failed = true; break; } // effectively a "continue;" for the outer loop
      }
      if (!failed) return root;
    }
  }

  // Essentially identical to the function above.
  // !!!WARNING: calls factor(P-1) so could be VERY SLOW!!!
  long PrimitiveRoot(const BigInt& P)
  {
    if (P <= 0 || !IsPrime(P))
      CoCoA_THROW_ERROR(ERR::BadArg, "PrimitiveRoot(P):  P must be a (positive) prime");

    if (P == 2) return 1;
    const factorization<BigInt> facpows = factor(P-1);
    const vector<BigInt>& primes = facpows.myFactors();
    const int NumPrimes = len(primes);
    for (unsigned long root = 2; /*empty*/; ++root)
    {
      CoCoA_ASSERT(root < P);
      if (root == 4 || root == 8 || root == 9 || root ==16) continue; // skip the first few prime powers
      bool failed = false;
      for (int i=0; i < NumPrimes; ++i)
      {
	if (PowerMod(root, (P-1)/primes[i], P) == 1)
        { failed = true; break; } // effectively a "continue;" for the outer loop
      }
      if (!failed) return root;
    }
  }


  //------------------------------------------------------------------
  // Quadratic residues...

  long KroneckerSymbol(const MachineInt& residue, const MachineInt& modulus)
  {
    if (IsNegative(modulus) || !IsSignedLong(modulus) || AsSignedLong(modulus) < 2) CoCoA_THROW_ERROR(ERR::BadModulus, "KroneckerSymbol");
    return mpz_kronecker_si(mpzref(BigInt(residue)), AsSignedLong(modulus));
  }

  long KroneckerSymbol(const BigInt& residue,     const MachineInt& modulus) 
  {
    if (IsNegative(modulus) || !IsSignedLong(modulus) || AsSignedLong(modulus) < 2) CoCoA_THROW_ERROR(ERR::BadModulus, "KroneckerSymbol");
    return mpz_kronecker_si(mpzref(residue), AsSignedLong(modulus));
  }
  
  long KroneckerSymbol(const MachineInt& residue, const BigInt& modulus)
  {
    if (modulus < 2) CoCoA_THROW_ERROR(ERR::BadModulus, "KroneckerSymbol");
    return mpz_si_kronecker(AsSignedLong(residue), mpzref(modulus));
  }

  long KroneckerSymbol(const BigInt& residue,     const BigInt& modulus)
  {
    if (modulus < 2) CoCoA_THROW_ERROR(ERR::BadModulus, "KroneckerSymbol");
    return mpz_kronecker(mpzref(residue), mpzref(modulus));
  }




} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/NumTheory-modular.C,v 1.5 2022/02/18 14:11:55 abbott Exp $
// $Log: NumTheory-modular.C,v $
// Revision 1.5  2022/02/18 14:11:55  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.4  2021/03/03 22:09:32  abbott
// Summary: New enum class (redmine 894)
//
// Revision 1.3  2021/02/10 19:40:00  abbott
// Summary: Added noexcept (sometimes instead of throw()) -- redmine 1572
//
// Revision 1.2  2020/06/17 15:49:24  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.1  2020/01/26 14:14:31  abbott
// Summary: Finished splitting NumTheory into smaller pieces (redming 1161)
//
//
//
