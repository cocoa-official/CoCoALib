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
