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
    myAssign(BigInt(1), BigInt(0), AlreadyReduced);  // AlreadyReduced diables check that denom is non-zero!!
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



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/BigRat.C,v 1.44 2022/10/10 09:53:34 abbott Exp $
// $Log: BigRat.C,v $
// Revision 1.44  2022/10/10 09:53:34  abbott
// Summary: Added move-assign op
//
// Revision 1.43  2022/08/05 11:54:19  abbott
// Summary: Changed ConvertToString to use internally size_t
//
// Revision 1.42  2022/02/18 14:11:53  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.41  2021/09/01 08:16:34  abbott
// Summary: Added comment
//
// Revision 1.40  2021/08/04 19:07:31  abbott
// Summary: Removed const (redmine 1606)
//
// Revision 1.39  2021/03/03 22:09:32  abbott
// Summary: New enum class (redmine 894)
//
// Revision 1.38  2021/02/17 17:50:36  abbott
// Summary: Now asserts ostream is in decimal mode (redmine 1547)
//
// Revision 1.37  2021/01/15 16:59:33  abbott
// Summary: Redmine 1563: new ctor for BigRat direct from integer
//
// Revision 1.36  2021/01/07 15:07:02  abbott
// Summary: Corrected copyright
//
// Revision 1.35  2020/12/05 13:09:06  abbott
// Summary: Cleaning
//
// Revision 1.34  2020/12/04 13:51:57  abbott
// Summary: Moved basic query fns to BigRat.H; updated BigRat.C (redmine 1529)
//
// Revision 1.33  2020/10/30 19:17:20  abbott
// Summary: Throw exc if istream is not good (redmine 1523)
//
// Revision 1.32  2020/10/05 19:23:43  abbott
// Summary: Change data member name to myRepr; added move ctor
//
// Revision 1.31  2020/06/17 15:49:21  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.30  2020/01/26 14:41:58  abbott
// Summary: Revised includes after splitting NumTheory (redmine 1161)
//
// Revision 1.29  2019/03/19 11:07:07  abbott
// Summary: Replaced 0 by nullptr where appropriate
//
// Revision 1.28  2019/03/18 11:09:56  abbott
// Summary: Added include (after changing NumTheory)
//
// Revision 1.27  2018/05/22 14:16:39  abbott
// Summary: Split BigRat into BigRat (class defn + ctors) and BigRatOps
//
// Revision 1.26  2018/05/18 12:13:36  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.25  2018/04/20 18:51:25  abbott
// Summary: Changed ctors for BigInt/BigRat from string or from MPZ/MPQ
//
// Revision 1.24  2018/04/18 14:16:18  abbott
// Summary: Ctor from mpq_t now checks if arg is nullptr, and if so, throws helpful error.
//
// Revision 1.23  2017/10/17 15:51:00  abbott
// Summary: Replaced gcd(...)==1 by IsCoprime(...)
//
// Revision 1.22  2016/11/11 14:15:32  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.21  2016/10/08 19:46:35  abbott
// Summary: Correct handling of oct/hex base when reading a "decimal"
//
// Revision 1.20  2016/08/02 12:49:34  abbott
// Summary: Renamed NumDigits to SizeInBase; updated doc.
//
// Revision 1.19  2016/07/21 14:14:22  abbott
// Summary: Changed op>> so it can read decimals
//
// Revision 1.18  2016/03/25 20:42:37  abbott
// Summary: Renamed utils_gmp to utils-gmp
//
// Revision 1.17  2016/03/25 19:58:40  abbott
// Summary: Added BigRat ctors from MantExp2 and MantExp10 structures
//
// Revision 1.16  2015/11/23 18:21:10  abbott
// Summary: Renamed ILogBase -> FloorLogBase; added FloorLog2, FloorLog10
//
// Revision 1.15  2015/10/09 18:26:36  abbott
// Summary: Corrected redmine reference
//
// Revision 1.14  2015/10/09 18:18:27  abbott
// Summary: Renamed "abs" to "uabs" for MachineInt; new fn "negate"; see redmine 783
//
// Revision 1.13  2015/09/17 18:10:53  abbott
// Summary: Cleaned impl of ILogBase; fixed redmine 776
//
// Revision 1.12  2014/07/09 11:41:35  abbott
// Summary: Corrected (embarassing) bug in ILogBase
// Author: JAA
//
// Revision 1.11  2014/07/07 12:11:12  abbott
// Summary: Corrected operator>> (forgot to ignore the '/')
// Author: JAA
//
// Revision 1.10  2014/06/14 19:25:11  abbott
// Summary: Added new fn CmpAbs (for BigRat)
// Author: JAA
//
// Revision 1.9  2014/05/16 12:02:28  abbott
// Summary: Changed comment about fn "round"
// Author: JAA
//
// Revision 1.8  2014/01/28 09:58:30  abbott
// Revised impl of ctor from std::string so that it accepts and respects 2nd arg saying whether the fraction should be canonicalized.
//
// Revision 1.7  2013/05/20 15:50:20  abbott
// Added new ctor for BigRat from mpq_t.
//
// Revision 1.6  2013/03/26 14:56:06  abbott
// Updated the conversion fns (in ptic removed procedure "convert");
// numerous consequential changes.
//
// Revision 1.5  2012/12/12 10:38:35  abbott
// Changed assertion to allow creation of 1/0 if marked as AlreadyReduced.
//
// Revision 1.4  2012/12/04 20:14:49  abbott
// Modified BigRat ctor to allow one to create 1/0 (if specified as AleadyReduced).
//
// Revision 1.3  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.2  2011/11/09 14:03:40  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.1  2011/09/23 13:20:35  bigatti
// -- QQ.C renamed into BigRat.C
//
// Revision 1.18  2011/09/06 15:21:53  abbott
// Changed "cmp" functions so that the return value is in {-1,0,+1}.
//
// Revision 1.17  2011/08/24 10:28:49  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.16  2011/08/23 16:18:38  abbott
// Simplified defn of round; added comment about rounding halves towards zero.
//
// Revision 1.15  2011/08/17 11:57:39  abbott
// Added static_cast to keep compiler quiet.
//
// Revision 1.14  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.13  2011/06/23 16:01:07  abbott
// Removed single arg ctor QQ(MachineInteger), & consequential changes.
//
// Revision 1.12  2011/03/01 15:26:10  abbott
// Improved impl of ILogBase -- faster in most cases.
//
// Revision 1.11  2011/02/25 12:06:51  abbott
// Added new fn IsOneNum; also some minor code cleaning in QQ.C
//
// Revision 1.10  2011/01/14 17:23:19  abbott
// Fixed a minor bug in power.
//
// Revision 1.9  2010/12/26 13:03:16  abbott
// Added ILogBase function (to ZZ & QQ).
//
// Revision 1.8  2010/05/07 14:57:52  abbott
// Two main changes:
//   power(QQ,ZZ) now allows negative exponent
//   renamed QQ::AlreadyNormalized to QQ::AlreadyReduced
//           (and allowed denoms to be negative; the ctor then makes them positive).
//
// Revision 1.7  2010/03/22 11:49:28  abbott
// Added ctor from a string.
//
// Revision 1.6  2010/03/18 16:40:42  abbott
// Added missing include directive.
//
// Revision 1.5  2010/03/18 16:34:10  abbott
// Added new pseudo-ctors for QQ with optional flag to indicate that value is already normalized.
// Added OpenMath I/O operators.
//
// Revision 1.4  2009/10/26 15:39:24  bigatti
// -- added CopyFromMPZ in ZZ ctor
//
// Revision 1.3  2009/07/08 12:26:53  abbott
// Added floor and ceil functions for QQs.
// Added example program for QQs.
// Minor correction to convert.C; minor cleaning to ex-ZZ1.C
//
// Revision 1.2  2009/07/06 12:31:26  abbott
// Commented out two unused function arguments (to keep compiler quiet
// when compiling a debugging version).
//
// Revision 1.1  2009/07/02 16:29:42  abbott
// Added new class QQ to represent rational numbers.
// Consequent change to the Makefile.
//
//
