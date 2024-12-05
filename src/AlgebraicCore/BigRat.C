//   Copyright (c)  2009-2010,2013  John Abbott and Anna M. Bigatti

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


#include "CoCoA/BigRat.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/FloatApprox.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/error.H"
#include "CoCoA/utils-gmp.H"
#include "CoCoA/utils.H"
#ifdef CoCoA_DEBUG
#include "CoCoA/NumTheory-gcd.H"
#endif

#include <cmath>
using std::abs;
using std::floor;
#include <iostream>
using std::ostream;
using std::istream;
#include <limits>
using std::numeric_limits;
#include <sstream>
using std::istringstream;
#include <string>
using std::string;
#include <vector>
using std::vector;

namespace CoCoA
{
  
  BigRat::BigRat()
  {
    mpq_init(myRepr);
  }


  BigRat::BigRat(const mpq_t q, CopyFromMPQ /*NotUsed*/)
  {
    if (q == nullptr)
      CoCoA_THROW_ERROR(ERR::NullPtr, "ctor BigRat(mpq_t)");

    mpq_init(myRepr);
    mpq_set(myRepr, q);
  }


  BigRat::BigRat(const MachineInt& n)
  {
    mpq_init(myRepr);
    if (IsNegative(n))
      mpq_set_si(myRepr, AsSignedLong(n), 1UL);
    else
      mpq_set_ui(myRepr, AsUnsignedLong(n), 1UL);
  }
  
  BigRat::BigRat(const BigInt& N)
  {
    mpq_init(myRepr);
    mpq_set_z(myRepr, mpzref(N));
  }

  BigRat::BigRat(const MachineInt& n1, const MachineInt& n2, ReduceFlag status)
  {
    if (IsZero(n2))
      CoCoA_THROW_ERROR(ERR::DivByZero, "BigRat(n1,n2)");
    mpq_init(myRepr);
    myAssign(BigInt(n1), BigInt(n2), status);
//    BELOW IS ORIGINAL CODE -- slightly better 'cos does not create temporaries.
//     const bool IsNegativeFraction = IsNegative(n1) ^ IsNegative(n2);
//     mpq_set_ui(myRepr, uabs(n1), uabs(n2));
//     if (status == NotReduced)
//       mpq_canonicalize(myRepr);
//     else
//       CoCoA_ASSERT(IsCoprime(n1,n2));
//     if (IsNegativeFraction)
//       mpq_neg(myRepr, myRepr);
  }

  BigRat::BigRat(const MachineInt& n1, const BigInt& N2, ReduceFlag status)
  {
    if (IsZero(N2))
      CoCoA_THROW_ERROR(ERR::DivByZero, "BigRat(n1,N2)");
    mpq_init(myRepr);
    myAssign(BigInt(n1), N2, status);
  }

  BigRat::BigRat(const BigInt& N1, const MachineInt& n2, ReduceFlag status)
  {
    if (IsZero(n2))
      CoCoA_THROW_ERROR(ERR::DivByZero, "BigRat(N1,n2)");
    mpq_init(myRepr);
    myAssign(N1, BigInt(n2), status);
  }

  BigRat::BigRat(const BigInt& N1, const BigInt& N2, ReduceFlag status)
  {
    if (IsZero(N2))
      CoCoA_THROW_ERROR(ERR::DivByZero, "BigRat(N1,N2)");
    mpq_init(myRepr);
    myAssign(N1, N2, status);
  }


  BigRat::BigRat(const std::string& str, ReadFromString /*NotUsed*/, ReduceFlag status)
  {
    mpq_init(myRepr);
//     if (base != 0 && (base < 2 || base > 36))
//       CoCoA_THROW_ERROR(ERR::BadNumBase, "BigRat(string,int)");
    constexpr int base = 10;
    if (mpq_set_str(myRepr, str.c_str(), base) != 0)
      CoCoA_THROW_ERROR(ERR::BadArg, "BigRat(string)");
    if (status == NotReduced)
      mpq_canonicalize(myRepr);
  }


  BigRat::BigRat(const MantExp2& ME)
  {
    mpq_init(myRepr);
    if (IsZero(ME.myMantissa)) return;
    const long exp = ME.myExponent-ME.myNumDigits+1;
    mpz_set(mpq_numref(myRepr), mpzref(ME.myMantissa));
    if (ME.mySign == -1) mpq_neg(myRepr, myRepr);
    if (exp >= 0)
      mpq_mul_2exp(myRepr, myRepr, exp);
    else
      mpq_div_2exp(myRepr, myRepr, -exp);
  }
  
  BigRat::BigRat(const MantExp10& ME)
  {
    mpq_init(myRepr);
    if (IsZero(ME.myMantissa)) return;
    const long exp = ME.myExponent-ME.myNumDigits+1;
    mpz_set(mpq_numref(myRepr), mpzref(ME.myMantissa));
    if (ME.mySign == -1) mpq_neg(myRepr, myRepr);
    if (exp >= 0)
    {
      const BigInt scale = power(10, exp);
      mpz_mul(mpq_numref(myRepr), mpq_numref(myRepr), mpzref(scale));
    }
    else
    {
      const BigInt scale = power(10,exp);
      mpz_set(mpq_denref(myRepr), mpzref(scale));
      mpq_canonicalize(myRepr);
    }
  }


  BigRat::BigRat(OneOverZero_t /*NotUsed*/)
  {
    mpq_init(myRepr);
    myAssign(BigInt(1), BigInt(0), AlreadyReduced);  // AlreadyReduced disables check that denom is non-zero!!
  }


  BigRat::BigRat(const BigRat& from) // std copy ctor
  {
    mpq_init(myRepr);
    mpq_set(myRepr, from.myRepr);
  }


  BigRat::BigRat(BigRat&& from)  /*noexcept*/ // std move ctor
  {
    mpq_init(myRepr);
    mpq_swap(myRepr, from.myRepr);
  }


  BigRat::~BigRat()
  {
    mpq_clear(myRepr);
  }


  // NOTE: This is NOT EXCEPTION CLEAN if the GMP fns can throw.
  void BigRat::myAssign(const BigInt& N1, const BigInt& N2, ReduceFlag status/*=NotReduced*/)
  {
    CoCoA_ASSERT(!IsZero(N2) || status == AlreadyReduced);
    const bool IsNegativeFraction = (N1 < 0) ^ (N2 < 0);
    mpz_abs(mpq_numref(myRepr), mpzref(N1));
    mpz_abs(mpq_denref(myRepr), mpzref(N2));
    if (status == NotReduced)
      mpq_canonicalize(myRepr);
    else
      CoCoA_ASSERT(IsCoprime(N1,N2));
    if (IsNegativeFraction)
      mpq_neg(myRepr, myRepr);
  }


  BigRat& BigRat::operator=(const BigRat& rhs)
  {
    mpq_set(myRepr, rhs.myRepr);
    return *this;
  }

  BigRat& BigRat::operator=(BigRat&& rhs)
  {
    mpq_swap(myRepr, rhs.myRepr);
    return *this;
  }

    // -------- functions that modify at least one argument or `*this' ----------

  BigRat& BigRat::operator+=(const BigRat& rhs)
  {
    mpq_add(myRepr, myRepr, rhs.myRepr);
    return *this;
  }

  BigRat& BigRat::operator-=(const BigRat& rhs)
  {
    mpq_sub(myRepr, myRepr, rhs.myRepr);
    return *this;
  }

  BigRat& BigRat::operator*=(const BigRat& rhs)
  {
    mpq_mul(myRepr, myRepr, rhs.myRepr);
    return *this;
  }

  BigRat& BigRat::operator/=(const BigRat& rhs)
  {
    if (mpz_sgn(mpq_numref(rhs.myRepr)) == 0)
      CoCoA_THROW_ERROR(ERR::DivByZero, "q1 /= q2");
    mpq_div(myRepr, myRepr, rhs.myRepr);
    return *this;
  }
                        
  // Same but with RHS a BigInt...
  BigRat& BigRat::operator=(const BigInt& rhs)
  {
    mpq_set_z(myRepr, mpzref(rhs));
    return *this;
  }


  // SLUG: impl is needlessly inefficient: makes useless copy in D
  BigRat& BigRat::operator+=(const BigInt& rhs)
  {
    const BigInt D = BigIntFromMPZ(mpq_denref(myRepr));
    const BigInt tmp = rhs*D;
    mpz_add(mpq_numref(myRepr), mpq_numref(myRepr), mpzref(tmp));
    // no need to call mpq_canonicalize
    return *this;
  }

  // SLUG: impl is needlessly inefficient: makes useless copy in D
  BigRat& BigRat::operator-=(const BigInt& rhs)
  {
    const BigInt D = BigIntFromMPZ(mpq_denref(myRepr));
    const BigInt tmp = rhs*D;
    mpz_sub(mpq_numref(myRepr), mpq_numref(myRepr), mpzref(tmp));
    // no need to call mpq_canonicalize
    return *this;
  }

  BigRat& BigRat::operator*=(const BigInt& rhs)
  {
    return operator*=(BigRat(rhs,1));
  }

  BigRat& BigRat::operator/=(const BigInt& rhs)
  {
    if (IsZero(rhs))
      CoCoA_THROW_ERROR(ERR::DivByZero, "Q /= N");
    // Could be more efficient if "*this" is 0.
    return operator/=(BigRat(rhs,1));
  }
                        
  // Same but with RHS a MachineInt...
  BigRat& BigRat::operator= (const MachineInt& rhs)
  {
    if (IsNegative(rhs))
      mpq_set_si(myRepr, AsSignedLong(rhs), 1);
    else
      mpq_set_ui(myRepr, AsUnsignedLong(rhs), 1);
    return *this;
  }

  BigRat& BigRat::operator+=(const MachineInt& rhs)
  {
    return operator+=(BigInt(rhs));
  }

  BigRat& BigRat::operator-=(const MachineInt& rhs)
  {
    return operator-=(BigInt(rhs));
  }

  BigRat& BigRat::operator*=(const MachineInt& rhs)
  {
    return operator*=(BigInt(rhs));
  }

  BigRat& BigRat::operator/=(const MachineInt& rhs)
  {
    if (IsZero(rhs))
      CoCoA_THROW_ERROR(ERR::DivByZero, "Q /= n");
    return operator/=(BigInt(rhs));
  }


  const BigRat& BigRat::operator++()
  {
    mpz_add(mpq_numref(myRepr), mpq_numref(myRepr), mpq_denref(myRepr)); // no need to reduce
    return *this;
  }

  const BigRat& BigRat::operator--()
  {
    mpz_sub(mpq_numref(myRepr), mpq_numref(myRepr), mpq_denref(myRepr));
    return *this;
  }

  const BigRat BigRat::operator++(int) // INEFFICIENT
  {
    BigRat ans(*this);
    operator++();
    return ans;
  }

  const BigRat BigRat::operator--(int) // INEFFICIENT
  {
    BigRat ans(*this);
    operator--();
    return ans;
  }



  // I/O FUNCTIONS

  string ConvertToString(const BigRat& src, int base/*=10*/)
  {
    if (base < 2 || base > 36)
      CoCoA_THROW_ERROR(ERR::BadNumBase, "IsConvertible(string,BigRat,int)");
    const std::size_t digits = mpz_sizeinbase(mpq_numref(mpqref(src)),base) +
                               mpz_sizeinbase(mpq_denref(mpqref(src)),base);
    vector<char> buffer(digits+3); // +3 to allow for minus sign, "/" character and terminating NUL
    mpq_get_str(&buffer[0], base, mpqref(src));
    return &buffer[0];
  }


  std::ostream& operator<<(std::ostream& out, const BigRat& Q)
  {
    if (!out) return out;  // short-cut for bad ostreams
    CoCoA_ASSERT(IsDecimal(out)); // this is also checked by output for BigInt
    out << num(Q);
    if (!IsOneDen(Q))
      out << "/" << den(Q);
    return out;
  }

  std::istream& operator>>(std::istream& in, BigRat& Q)
  {
    static const char* const FnName = "operator>> for BigRat";
    if (!in.good()) CoCoA_THROW_ERROR("istream is not good", FnName);
    if (!IsDecimal(in)) CoCoA_THROW_ERROR("istream is not in \"decimal\" mode", FnName);
    in.peek(); if (in.eof()) CoCoA_THROW_ERROR("EOF", FnName); // so that err mesg refers to correct fn
    BigInt N;
    in >> N;
    Q = N;
    const char SlashOrDot = in.peek();  // might trigger EOF
    if (in.eof()) { in.clear(); return in; }
    if (SlashOrDot != '/' && SlashOrDot != '.')
    {
      return in;
    }
    in.ignore(); // cannot trigger EOF
    if (SlashOrDot == '.')
    {
      const string AfterDot = ScanUnsignedIntegerLiteral(in);
      const long NumPlaces = len(AfterDot);
      if (NumPlaces == 0) return in;

      istringstream FracDigits(AfterDot);
      FracDigits >> N; // N now contains the "fractional decimal part"
      constexpr int base = 10;
      Q += BigRat(N, power(base, NumPlaces));
      return in;
    }

    // Found a slash
    const char AfterSlash = in.peek(); // might trigger EOF
    if (!in.good() || !isdigit(AfterSlash))
      CoCoA_THROW_ERROR("Missing denominator in rational", FnName);
    BigInt D;
    in >> D;
    Q /= D; // Might throw DivByZero
    return in;
  }


  OpenMathOutput& operator<<(OpenMathOutput& OMOut, const BigRat& Q)
  {
    OMOut->mySendApplyStart();
    OMOut->mySendSymbol("nums1","rational");
    OMOut->mySend(num(Q));
    OMOut->mySend(den(Q));
    OMOut->mySendApplyEnd();
    return OMOut;
  }


  OpenMathInput& operator>>(OpenMathInput& OMIn, BigRat& /*Q*/)
  {
    CoCoA_THROW_ERROR(ERR::NYI, "OpenMathInput fn for BigRat");
    return OMIn;
  }


  void swap(BigRat& a, BigRat& b)
  {
    mpq_swap(mpqref(a), mpqref(b));
  }

  BigInt num(const BigRat& Q)
  {
    return BigIntFromMPZ(mpq_numref(mpqref(Q)));
  }

  BigInt den(const BigRat& Q)
  {
    return BigIntFromMPZ(mpq_denref(mpqref(Q)));
  }


} // end of namespace CoCoA
