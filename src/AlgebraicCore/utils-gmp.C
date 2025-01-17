//   Copyright (c)  2014,2016,2021  John Abbott and Anna M. Bigatti

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


#include "CoCoA/utils-gmp.H"
#include "CoCoA/convert.H"

#include <cstddef>
//using std::size_t;

namespace CoCoA
{

  // Equiv to cmp(abs(q1), abs(q2))
  int mpq_cmpabs(const mpq_t q1, const mpq_t q2)
  {
    if (mpq_sgn(q1) == 0) { if (mpq_sgn(q2) == 0) return 0; else return -1; }
    if (mpq_sgn(q2) == 0) return 1;

    if (mpz_cmp_ui(mpq_denref(q1),1) == 0 && mpz_cmp_ui(mpq_denref(q2),1) == 0)
      return mpz_cmpabs(mpq_numref(q1), mpq_numref(q2));

    const std::size_t logn1 = mpz_sizeinbase(mpq_numref(q1),2); // logn1-1 <= log(n1) < logn1
    const std::size_t logd1 = mpz_sizeinbase(mpq_denref(q1),2); //
    const std::size_t logn2 = mpz_sizeinbase(mpq_numref(q2),2); //
    const std::size_t logd2 = mpz_sizeinbase(mpq_denref(q2),2); //

    if (logn1+logd2-2 >= logn2+logd1) return 1; //??? overflow???
    if (logn1+logd2 <= logn2+logd1-2) return -1;//??? overflow???
    // IDEA: could compare more accurate logs of n1*d2 and n2*d1.
    mpz_t n1d2; mpz_init(n1d2);
    mpz_t n2d1; mpz_init(n2d1);
    mpz_mul(n1d2,mpq_numref(q1),mpq_denref(q2));
    mpz_mul(n2d1,mpq_numref(q2),mpq_denref(q1));
    const int ans = mpz_cmpabs(n1d2,n2d1);
    mpz_clear(n2d1);
    mpz_clear(n1d2);
    return ans;
  }


  // Sets Q = N/D using rounding;
  // rounding is symmetric wrt 0: either towards 0 or away from 0
  bool mpz_rounddiv(mpz_t Q, const mpz_t N, const mpz_t D)
  {
    mpz_t remainder; mpz_init(remainder);
    mpz_fdiv_qr(Q, remainder, N, D);
    if (mpz_sgn(remainder) == 0) { mpz_clear(remainder); return true; }
    mpz_mul_2exp(remainder, remainder, 1);
    const int cmp = mpz_cmpabs(remainder, D);
    mpz_clear(remainder);
    if ((cmp == 1) || (cmp == 0 && mpz_sgn(Q) >= 0))  // halves round away from 0
//    if ((cmp == 1) || (cmp == 0 && mpz_sgn(Q) < 0))  // halves round towards 0
      mpz_add_ui(Q, Q, 1);
    return false;
  }

  // N = rounded(Q); rounding of halves is decided by mpz_rounddiv.
  bool mpq_round(mpz_t N, const mpq_t Q)
  {
    return mpz_rounddiv(N, mpq_numref(Q), mpq_denref(Q));
  }
  

  double mpq_get_d_2exp(signed long int* exp, const mpq_t Q)
  {
    constexpr long MantBits = 63; // DBL_MANT_DIG, or std::numeric_limits<double>::digits (assuming std::numerics_limits<double>:: radix == 2)
    constexpr long WordSize = 64; // CHAR_BIT*sizeof(long), or std::numeric_limits<unsigned long>::digits
    mpz_t N;
    mpz_t D;
    // Short-cut for triv cases where Q is an integer
///???    if (mpq_sgn(Q) == 0) {*exp = 0; return 0.0; }
    if (mpz_cmp_ui(mpq_denref(Q), 1) == 0)
      return mpz_get_d_2exp(exp, mpq_numref(Q));
    mpz_init(N);  mpz_init(D);
    long log2num = IntegerCast<long>(mpz_sizeinbase(mpq_numref(Q),2));
    long log2den = IntegerCast<long>(mpz_sizeinbase(mpq_denref(Q),2));
    // We know that abs(Q) > 2^r  where r = log2num-log2den-1
    // Find pwr of 2 s.t. abs(Q*2^shift) >= 2^WordSize (but not too much bigger)
    long BinaryShift = 1+MantBits+(log2den-log2num);  // overflow safe???
    if (BinaryShift <= 0)
    {
      BinaryShift = -WordSize*((-BinaryShift)/WordSize); // make BinaryShift a mult of WordSize (optional)
      mpz_set(N, mpq_numref(Q));
      mpz_mul_2exp(D, mpq_denref(Q), -BinaryShift);
    }
    else if (BinaryShift > 0)
    {
      BinaryShift = WordSize*(1+(BinaryShift-1)/WordSize); // make BinaryShift a mult of WordSize (optional)
      mpz_set(D, mpq_denref(Q));
      mpz_mul_2exp(N, mpq_numref(Q), BinaryShift);
    }
    mpz_tdiv_q(N,N,D);
    long e;
    const double m = mpz_get_d_2exp(&e,N);
    mpz_clear(D);  mpz_clear(N);
    *exp = e - BinaryShift;
    return m;
  }


  // Shorter: less accurate, but does not alloc memory!
  // THIS WILL GIVE (slightly) WRONG ANS for e.g. BigRat(2^98-1, 2^99-1)
  // double mpq_get_d_2exp(signed long int* e, mpq_srcptr q)
  // {
  //   long expnum, expden;
  //   double mantnum, mantden;
  //   if (mpq_sgn(q) == 0) { e = 0; return 0.0; }
  //   if (mpz_cmp_ui(mpq_denref(Q), 1) == 0)
  //     return mpz_get_d_2exp(exp, mpq_numref(Q));
  //   mantnum = mpz_get_d_2exp(&expnum, mpq_numref(q)); // ???expnum could overflow
  //   mantden = mpz_get_d_2exp(&expden, mpq_denref(q)); // ???expden could overflow
  //   *e = expnum-expden;
  //   double ans = mantnum/mantden;
  //   // Renormalize ans; special case if ans == 1
  //   if (abs(ans) == 1 && mpz_cmpabs(mpq_numref(q), mpq_denref(q)) < 0)
  //   { ans = nextafter(ans,0.0); }
  //   else if (abs(ans) >= 1) { ans /= 2; ++*e; } //???overflow in e???
  //   return ans;
  // }


  // 2021-09-13 NOT SURE THIS IS SUCH A GOOD IDEA
  // bool mpz_mul_check(mpz_t ans, const mpz_t A, const mpz_t B)
  // {
  //   const long logA = mpz_sizeinbase(A,2);
  //   const long logB = mpz_sizeinbase(B,2);
  //   if (logA > OVERFLOW_BITS || logB > OVERFLOW_BITS || logA > OVERFLOW_BITS-logB)
  //   {
  //     return false; // computation failed, would overflow
  //   }
  //   mpz_mul(ans, A, B);
  //   return true;
  // }



} // end of namespace CoCoA
