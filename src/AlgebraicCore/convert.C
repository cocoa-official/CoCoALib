//   Copyright (c)  2007,2009  John Abbott and Anna M. Bigatti

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

#include "CoCoA/convert.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRatOps.H"
#include "CoCoA/ring.H"

#include <limits>
using std::numeric_limits;
#include <cmath>
using std::ldexp;
using std::floor;

namespace CoCoA
{

  //-------------------------------------------------------
  // from BigInt

  bool IsConvertible(long& lhs, const BigInt& src) noexcept
  {
    if (!mpz_fits_slong_p(mpzref(src))) return false;
    lhs = mpz_get_si(mpzref(src));
    return true;
  }

  bool IsConvertible(int& lhs, const BigInt& src) noexcept
  {
    if (!mpz_fits_sint_p(mpzref(src))) return false;
    lhs = mpz_get_si(mpzref(src));
    return true;
  }

  bool IsConvertible(unsigned long& lhs, const BigInt& src) noexcept
  {
    if (!mpz_fits_ulong_p(mpzref(src))) return false;
    lhs = mpz_get_ui(mpzref(src));
    return true;
  }

  bool IsConvertible(unsigned int& lhs, const BigInt& src) noexcept
  {
    if (!mpz_fits_uint_p(mpzref(src))) return false;
    lhs = mpz_get_ui(mpzref(src));
    return true;
  }


  //-------------------------------------------------------
  // from BigRat

  bool IsConvertible(long& lhs, const BigRat& src)
  {
    // Simple rather than fast -- makes unnecessary copy of num(src)
    return IsOneDen(src) && IsConvertible(lhs, num(src));
  }

  bool IsConvertible(int& lhs, const BigRat& src)
  {
    // Simple rather than fast -- makes unnecessary copy of num(src)
    return IsOneDen(src) && IsConvertible(lhs, num(src));
  }

  bool IsConvertible(unsigned long& lhs, const BigRat& src)
  {
    // Simple rather than fast -- makes unnecessary copy of num(src)
    return IsOneDen(src) && IsConvertible(lhs, num(src));
  }

  bool IsConvertible(unsigned int& lhs, const BigRat& src)
  {
    // Simple rather than fast -- makes unnecessary copy of num(src)
    return IsOneDen(src) && IsConvertible(lhs, num(src));
  }


  //-------------------------------------------------------
  // double with long, BigInt, BigRat

  bool IsConvertible(long& n, double z) noexcept
  {
    //??? BUG should handle infinity and Nan???
    if (z != std::floor(z)) return false;
    constexpr double MIN = static_cast<double>(numeric_limits<long>::min());
    constexpr double MAX = static_cast<double>(numeric_limits<long>::max());
    if (z < MIN || z > MAX) return false;
    n = static_cast<long>(z); // already checked that z is integer and in range.
    return true;
  }


  // NB ALL sufficiently large doubles are considered to be integers
  bool IsConvertible(BigInt& N, double z)
  {
    //??? BUG should handle infinity and Nan???
    if (z != std::floor(z)) return false;
    mpz_set_d(mpzref(N), z);
    return true;
  }




  bool IsConvertible(BigRat& Q, double z)
  {
    //??? BUG should handle infinity and Nan???
    if (z == 0) { Q = 0; return true; }
    mpq_set_d(mpqref(Q), z);  // BUG Don't know what happens if this fails.
    return true;
  }


  // Sets z = N if possible.
  bool IsConvertible(double& z, const BigInt& N)
  {
    if (IsZero(N)) { z = 0.0; return true; }
    // I prefer to do this manually because mpz_get_d does not give
    // good guarantees about what happens when overflow occurs.
    long ExpN;
    double MantissaN = mpz_get_d_2exp(&ExpN, mpzref(N));
    if (ExpN > numeric_limits<double>::max_exponent)
    {
      return false;
    }
    z = ldexp(MantissaN, ExpN);
    return true;
  }


  // Sets z = num/den if possible.
  // Underflow sets z to 0 and returns true.
  // Overflow returns false (and does not change z).
  bool IsConvertible(double& z, const BigRat& Q)
  {
    if (IsZero(Q)) { z = 0.0; return true; }
    // I prefer to do this manually rather than use mpq_get_d because
    // num and den may not be coprime, and the cost of computing the
    // gcd could be high.  Also mpq_get_d does not give nice guarantees
    // about what happens in the case of over-/under-flow.
    long Nexp;
    double N = mpz_get_d_2exp(&Nexp,mpq_numref(mpqref(Q)));
    long Dexp;
    double D = mpz_get_d_2exp(&Dexp, mpq_denref(mpqref(Q)));
    long ExpAns = Nexp-Dexp;
    double MantissaAns = N/D;
    if (MantissaAns >= 1 || MantissaAns <= -1)
    { MantissaAns /= 2; ++ExpAns; }

    // I believe this is the right interpretation of the exponent limits...
    if (ExpAns < numeric_limits<double>::min_exponent)
    {
      // We have underflow, so convert to 0 and indicate "success"
      z = 0.0;
      return true;
    }

    if (ExpAns > numeric_limits<double>::max_exponent)
    {
      // We have overflow, so indicate "failure".
      return false;
    }

    z = ldexp(MantissaAns, ExpAns); // Cannot overflow or underflow.
    return true;
  }


  //-------------------------------------------------------
  // RingElem to long, BigInt, BigRat

  bool IsConvertible(long& n, ConstRefRingElem x)
  {
    BigInt N;
    return (IsConvertible(N, x) && IsConvertible(n, N));
  }

  bool IsConvertible(BigInt& N, ConstRefRingElem x)
  {
    return IsInteger(N, x);
  }

  bool IsConvertible(BigRat& N, ConstRefRingElem x)
  {
    return IsRational(N, x);
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/convert.C,v 1.29 2023/03/27 14:01:12 abbott Exp $
// $Log: convert.C,v $
// Revision 1.29  2023/03/27 14:01:12  abbott
// Summary: Added explicit casts to silence compiler warnings; better readability also?
//
// Revision 1.28  2022/02/18 14:12:02  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.27  2021/02/10 19:40:00  abbott
// Summary: Added noexcept (sometimes instead of throw()) -- redmine 1572
//
// Revision 1.26  2021/01/07 15:23:38  abbott
// Summary: Corrected copyright
//
// Revision 1.25  2018/05/22 14:16:40  abbott
// Summary: Split BigRat into BigRat (class defn + ctors) and BigRatOps
//
// Revision 1.24  2018/05/18 12:25:54  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.23  2014/04/02 15:55:45  abbott
// Summary: Reorganized file contents (more sensible order now, I hope)
// Author: JAA
//
// Revision 1.22  2013/03/26 14:56:06  abbott
// Updated the conversion fns (in ptic removed procedure "convert");
// numerous consequential changes.
//
// Revision 1.21  2013/02/15 16:32:44  abbott
// Removed IsInRange (moved to MachineInt).
// Added IsConvertible(BigInt, RingElem) and IsConvertible(BigRat, RingElem).
//
// Revision 1.20  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.19  2011/11/09 15:35:00  bigatti
// -- moved include MachineInt into convert.H
//
// Revision 1.18  2011/11/09 14:29:37  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.17  2011/08/27 20:50:47  abbott
// Added include directive for MachineInt.
//
// Revision 1.16  2011/08/24 10:32:04  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.15  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.14  2011/08/12 16:08:32  abbott
// Hastily added conversion fns for BigInt (like those for ZZ).
//
// Revision 1.13  2011/03/22 20:26:12  abbott
// Added  IsConvertible(long&, double).
// Activated commented out ConvertTo template.
// Checking in because some other files need the new code.
//
// Revision 1.12  2010/04/23 13:21:45  abbott
// Conversion from string to ZZ now defaults to decimal (and not the quirky C standard).
//
// Revision 1.11  2010/03/22 11:51:37  abbott
// Simplified convert fn from string to BigRat -- now uses GMP's implementation.
//
// Revision 1.10  2009/12/23 22:27:28  abbott
// Added conversions from BigRat to machine integers.
//
// Revision 1.9  2009/12/23 18:53:51  abbott
// Major change to conversion functions:
//   convert(..) is now a procedure instead of a function
//   IsConvertible(..) replaces the former convert(..) function
//   Added new NumericCast conversion function (placeholder for BOOST feature)
//   Consequent changes in code which uses these features.
//
// Revision 1.8  2009/12/11 11:46:32  abbott
// Changed fn  convert  into  IsConvertible.
// Added template procedure  convert.
// New version because change is not backward compatible.
//
// Revision 1.7  2009/07/08 12:26:53  abbott
// Added floor and ceil functions for BigRats.
// Added example program for BigRats.
// Minor correction to convert.C; minor cleaning to ex-ZZ1.C
//
// Revision 1.6  2009/07/02 16:32:11  abbott
// Consequential changes stemming from new class BigRat, and modified interface to the member
// function RingBase::myIsRational.  Also some new conversion functions.
//
// Revision 1.5  2008/10/07 15:45:22  abbott
// Changed ErrorInfo objects so they include the name of their own error ID.
// Changed catch statements to catch const objects.
// Removed calls to the member fn which accessed the error ID member of an
// ErrorInfo; now you simply compare directly with the error ID (makes the
// code easier to read).
//
// Revision 1.4  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.3  2007/05/22 22:51:39  abbott
// Changed name of fn ndigits to NumDigits.
// Changed return type of exponent and NumDigits.
// Changed some exceptions from nonstandard to the appropriate standard one.
//
// Revision 1.2  2007/03/28 15:34:03  abbott
// Underflow during conversion to a double no longer causes "failure" but now
// converts the value to zero and yields "success".
//
// Revision 1.1  2007/03/23 18:38:42  abbott
// Separated the "convert" functions (and CheckedCast) into their own files.
// Many consequential changes.  Also corrected conversion to doubles.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.3  2007/03/03 15:24:39  cocoa
// Changed 2006 to 2007.
//
// Revision 1.2  2006/10/06 14:04:14  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
