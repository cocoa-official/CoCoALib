//   Copyright (c)  2009-2010,2013,2018  John Abbott and Anna Bigatti
//   Contributions from: Vibhash Singh

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


#include "CoCoA/BigRatOps.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/FloatApprox.H"
#include "CoCoA/NumTheory-gcd.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/error.H"
#include "CoCoA/utils.H"
#include "CoCoA/utils-gmp.H"

#include <cmath>
using std::abs;
using std::floor;
#include <limits>
using std::numeric_limits;
#include <vector>
using std::vector;

namespace CoCoA
{
  

  // STANDARD ARITHMETIC OPERATIONS

  BigRat abs(const BigRat& Q)
  {
    BigRat ans;
    mpq_abs(mpqref(ans), mpqref(Q));
    return ans;
  }

  BigRat operator-(const BigRat& Q)
  {
    BigRat ans;
    mpq_neg(mpqref(ans), mpqref(Q));
    return ans;
  }

  BigRat operator+(const BigRat& Q1, const BigRat& Q2)
  {
    BigRat ans;
    mpq_add(mpqref(ans), mpqref(Q1), mpqref(Q2));
    return ans;
  }

  BigRat operator-(const BigRat& Q1, const BigRat& Q2)
  {
    BigRat ans;
    mpq_sub(mpqref(ans), mpqref(Q1), mpqref(Q2));
    return ans;
  }

  BigRat operator*(const BigRat& Q1, const BigRat& Q2)
  {
    BigRat ans;
    mpq_mul(mpqref(ans), mpqref(Q1), mpqref(Q2));
    return ans;
  }

  BigRat operator/(const BigRat& Q1, const BigRat& Q2)
  {
    if (IsZero(Q2))
      CoCoA_THROW_ERROR1(ERR::DivByZero);
    BigRat ans;
    mpq_div(mpqref(ans), mpqref(Q1), mpqref(Q2));
    return ans;
  }


  BigRat operator+(const BigRat& Q, const BigInt& N)
  {
    BigRat ans = Q;
    return ans += N;
//   THE LINES BELOW SHOULD BE MORE EFFICIENT (but they're not very readable).
//     BigRat ans;
//     mpz_mul(mpq_denref(ans), mpzref(N), mpq_denref(mpqref(Q)));
//     mpz_add(mpq_numref(mpqref(ans)), mpq_denref(mpqref(ans)), mpq_numref(mpqref(Q)));
//     mpz_set(mpq_denref(mpqref(ans)), mpq_denref(mpqref(Q)));
//     return ans;
  }

  BigRat operator-(const BigRat& Q, const BigInt& N)
  {
    BigRat ans(Q);
    return ans -= N;
  }

  BigRat operator*(const BigRat& Q, const BigInt& N)
  {
    BigRat ans(Q);
    return ans *= N;
  }

  BigRat operator/(const BigRat& Q, const BigInt& N)
  {
    if (IsZero(N))
      CoCoA_THROW_ERROR1(ERR::DivByZero);
    BigRat ans(Q);
    return ans /= N;
  }

  BigRat operator+(const BigInt& N, const BigRat& Q)
  {
    return Q+N;
  }

  BigRat operator-(const BigInt& N, const BigRat& Q)
  {
    BigRat ans = Q-N;
    return -ans;
  }

  BigRat operator*(const BigInt& N, const BigRat& Q)
  {
    return Q*N;
  }

  BigRat operator/(const BigInt& N, const BigRat& Q)
  {
    if (IsZero(Q))
      CoCoA_THROW_ERROR1(ERR::DivByZero);
    return BigRat(N,1)/Q;
  }


  BigRat operator+(const BigRat& Q, const MachineInt& n)
  {
    return Q + BigRat(n,1);
  }

  BigRat operator-(const BigRat& Q, const MachineInt& n)
  {
    return Q - BigRat(n,1);
  }

  BigRat operator*(const BigRat& Q, const MachineInt& n)
  {
    return Q * BigRat(n,1);
  }

  BigRat operator/(const BigRat& Q, const MachineInt& n)
  {
    if (IsZero(n))
      CoCoA_THROW_ERROR1(ERR::DivByZero);
    return Q / BigRat(n,1);
  }


  BigRat operator+(const MachineInt& n, const BigRat& Q)
  {
    return BigRat(n,1) + Q;
  }

  BigRat operator-(const MachineInt& n, const BigRat& Q)
  {
    return BigRat(n,1) - Q;
  }

  BigRat operator*(const MachineInt& n, const BigRat& Q)
  {
    return BigRat(n,1) * Q;
  }

  BigRat operator/(const MachineInt& n, const BigRat& Q)
  {
    if (IsZero(Q))
      CoCoA_THROW_ERROR1(ERR::DivByZero);
    if (IsOne(n))  { BigRat ans; mpq_inv(mpqref(ans), mpqref(Q)); return ans; }
    return BigRat(n,1) / Q;
  }

  BigRat power(const BigRat& base, const BigInt& exponent)
  {
    if (exponent >= 0)
      return BigRat(power(num(base), exponent), power(den(base), exponent), BigRat::AlreadyReduced);
    if (IsZero(base))
      CoCoA_THROW_ERROR1(ERR::BadPwrZero);
    return BigRat(power(den(base), -exponent), power(num(base), -exponent), BigRat::AlreadyReduced);
  }

  BigRat power(const BigRat& base, const MachineInt& exponent)
  {
    if (!IsNegative(exponent))
      return BigRat(power(num(base), exponent), power(den(base), exponent), BigRat::AlreadyReduced);
    if (IsZero(base))
      CoCoA_THROW_ERROR1(ERR::BadPwrZero);
    return BigRat(power(den(base), uabs(exponent)), power(num(base), uabs(exponent)), BigRat::AlreadyReduced);
  }



//--------------------------------------------
// SumBigRat code developed from code by Vibhash Singh (Uni Passau)
  
  inline std::size_t SumBigRat::ourIndex(const BigRat& Q) noexcept
  {
    const std::size_t sz = mpz_size(mpq_denref(mpqref(Q)));
    return FloorLog2(sz); 
  }


  SumBigRat& SumBigRat::operator+=(long n)
  {
    if (n != 0)
      myBigIntPart += n;
    return *this;
  }

  SumBigRat& SumBigRat::operator+=(const BigInt& N)
  {
    if (!IsZero(N))
      myBigIntPart += N;
    return *this;
  }

  SumBigRat& SumBigRat::operator+=(const BigRat& Q)
  {
    if (IsZero(Q))  return *this;
    BigInt IntPart;  // Values computed via the mpz_ fns right below
    BigRat FracPart; //
    mpz_set(mpq_denref(mpqref(FracPart)), mpq_denref(mpqref(Q)));
    mpz_tdiv_qr(mpzref(IntPart), mpq_numref(mpqref(FracPart)),
                mpq_numref(mpqref(Q)), mpq_denref(mpqref(Q)));
    myBigIntPart += IntPart;
    const std::size_t i = ourIndex(FracPart);
    const bool MustExtend = (i >= myGeobucket.size());
    if (MustExtend)  myGeobucket.resize(i+1);
    myGeobucket[i] += FracPart;
    if (!MustExtend)  myCascade(i);
    return *this;
  }


  // void SumBigRat::myAdd(const BigInt& N) // threadsafe version of op+=
  // {
  //   //    std::scoped_lock lock(myMutex); // needs C++17
  //   std::lock_guard<std::mutex> lock(myMutex); // cannot use auto :-(
  //   this->operator+=(N);
  // }

  // void SumBigRat::myAdd(const BigRat& Q) // threadsafe version of op+=
  // {
  //   //    std::scoped_lock lock(myMutex); // needs C++17
  //   std::lock_guard<std::mutex> lock(myMutex); // cannot use auto :-(
  //   this->operator+=(Q);
  // }


  void SumBigRat::myCascade(std::size_t i)
  {
    while (true)
    {
      const std::size_t index = ourIndex(myGeobucket[i]);
      if (index == i) return;
      if (index >= myGeobucket.size())
      {
        myGeobucket.resize(index+1);
        swap(myGeobucket[i], myGeobucket[index]); // sets myGeobucket[i]=0
        return;
      }
      myGeobucket[index] += myGeobucket[i];
      myGeobucket[i] = 0;
      i = index;
    }
  }


  BigRat SumBigRat::myTotal() const
  {
//???    //    std::scoped_lock lock(myMutex); // needs C++17
//???    std::lock_guard<std::mutex> lock(myMutex); // cannot use auto :-(
    BigRat sum;
    for (BigRat& q: myGeobucket)
    {
      sum += q;
      q = 0;
    }
    if (IsZero(sum))
      return BigRat(myBigIntPart.myTotal());

    const std::size_t i = ourIndex(sum);
    if (i >= myGeobucket.size())  myGeobucket.resize(i+1);
    myGeobucket[i] = sum;
    return sum + total(myBigIntPart);
  }


  std::ostream& operator<<(std::ostream& out, const SumBigRat& S)
  {
    if (!out)  return out;
    out << "SumBigRat(myBigIntPart=" << S.myBigIntPart << " bits, myGeobucket="<< len(S.myGeobucket) <<"buckets)";
    return out;
  }


// End of code about SumBigRat
//--------------------------------------------

  // COMPARISON FUNCTIONS

  int cmp(const BigRat& Q1, const BigRat& Q2)
  {
    return sign(mpq_cmp(mpqref(Q1), mpqref(Q2)));
  }


  int cmp(const BigRat& Q, const BigInt& N)
  {
    return cmp(num(Q), N*den(Q));
  }


  int cmp(const BigInt& N, const BigRat& Q)
  {
    return cmp(N*den(Q), num(Q));
  }


  int cmp(const BigRat& Q, const MachineInt& n)
  {
    if (IsNegative(n))
      return sign(mpq_cmp_si(mpqref(Q), AsSignedLong(n),1));
    return sign(mpq_cmp_ui(mpqref(Q), AsUnsignedLong(n),1));
  }


  int cmp(const MachineInt& n, const BigRat& Q)
  {
    return -cmp(Q, n);
  }


  int CmpAbs(const BigRat& Q1, const BigRat& Q2)
  {
    return mpq_cmpabs(mpqref(Q1), mpqref(Q2));
  }

  // The next 4 prefer simplicity over speed.
  int CmpAbs(const BigRat& Q, const BigInt& N)
  {
    return CmpAbs(Q, BigRat(N,1));
  }

  int CmpAbs(const BigInt& N, const BigRat& Q)
  {
    return CmpAbs(BigRat(N,1), Q);
  }

  int CmpAbs(const BigRat& Q, const MachineInt& n)
  {
    return CmpAbs(Q, BigRat(n,1));
  }

  int CmpAbs(const MachineInt& n, const BigRat& Q)
  {
    return CmpAbs(BigRat(n,1), Q);
  }



  bool operator==(const BigRat& Q1, const BigRat& Q2)
  {
    return mpq_equal(mpqref(Q1), mpqref(Q2));
  }

  bool operator!=(const BigRat& Q1, const BigRat& Q2)
  {
    return !(Q1 == Q2);
  }

  bool operator> (const BigRat& Q1, const BigRat& Q2)
  {
    return cmp(Q1,Q2) > 0;
  }

  bool operator>=(const BigRat& Q1, const BigRat& Q2)
  {
    return cmp(Q1,Q2) >= 0;
  }

  bool operator< (const BigRat& Q1, const BigRat& Q2)
  {
    return cmp(Q1,Q2) < 0;
  }

  bool operator<=(const BigRat& Q1, const BigRat& Q2)
  {
    return cmp(Q1,Q2) <= 0;
  }

                        
  bool operator==(const BigRat& Q, const BigInt& N)
  {
    return IsOneDen(Q) && num(Q) == N;
  }

  bool operator!=(const BigRat& Q, const BigInt& N)
  {
    return !(Q == N);
  }

  bool operator> (const BigRat& Q, const BigInt& N)
  {
    return cmp(Q,N) > 0;
  }

  bool operator>=(const BigRat& Q, const BigInt& N)
  {
    return cmp(Q,N) >= 0;
  }

  bool operator< (const BigRat& Q, const BigInt& N)
  {
    return cmp(Q,N) < 0;
  }

  bool operator<=(const BigRat& Q, const BigInt& N)
  {
    return cmp(Q,N) <= 0;
  }

                        
  bool operator==(const BigInt& N, const BigRat& Q)
  {
    return Q == N;
  }

  bool operator!=(const BigInt& N, const BigRat& Q)
  {
    return !(Q == N);
  }

  bool operator> (const BigInt& N, const BigRat& Q)
  {
    return cmp(N,Q) > 0;
  }

  bool operator>=(const BigInt& N, const BigRat& Q)
  {
    return cmp(N,Q) >= 0;
  }

  bool operator< (const BigInt& N, const BigRat& Q)
  {
    return cmp(N,Q) < 0;
  }

  bool operator<=(const BigInt& N, const BigRat& Q)
  {
    return cmp(N,Q) <= 0;
  }

                        
  bool operator==(const BigRat& Q, const MachineInt& n)
  {
    return IsOneDen(Q) && num(Q) == n;
  }

  bool operator!=(const BigRat& Q, const MachineInt& n)
  {
    return !(Q == n);
  }

  bool operator> (const BigRat& Q, const MachineInt& n)
  {
    return cmp(Q,n) > 0;
  }

  bool operator>=(const BigRat& Q, const MachineInt& n)
  {
    return cmp(Q,n) >= 0;
  }

  bool operator< (const BigRat& Q, const MachineInt& n)
  {
    return cmp(Q,n) < 0;
  }

  bool operator<=(const BigRat& Q, const MachineInt& n)
  {
    return cmp(Q,n) <= 0;
  }

                
  bool operator==(const MachineInt& n, const BigRat& Q)
  {
    return IsOneDen(Q) && num(Q) == n;
  }

  bool operator!=(const MachineInt& n, const BigRat& Q)
  {
    return !(n == Q);
  }

  bool operator> (const MachineInt& n, const BigRat& Q)
  {
    return cmp(n,Q) > 0;
  }

  bool operator>=(const MachineInt& n, const BigRat& Q)
  {
    return cmp(n,Q) >= 0;
  }

  bool operator< (const MachineInt& n, const BigRat& Q)
  {
    return cmp(n,Q) < 0;
  }

  bool operator<=(const MachineInt& n, const BigRat& Q)
  {
    return cmp(n,Q) <= 0;
  }

                        

  // MISCELLANEOUS FUNCTIONS

  double mantissa(long& exp, const BigRat& Q) noexcept
  {
    return mpq_get_d_2exp(&exp, mpqref(Q));
  }

  // double mantissa(const BigRat& Q)
  // {
  //   long DiscardExponent;
  //   return mpq_get_d_2exp(&DiscardExponent, mpqref(Q));
  // }

  long exponent(const BigRat& Q)
  {
    long exp;
    /*Discard mantissa*/ mpq_get_d_2exp(&exp, mpqref(Q));
    return exp;
  }


  double log(const BigRat& Q)
  {
    if (sign(Q) <= 0)
      CoCoA_THROW_ERROR1(ERR::ReqPositive);
    long exp;
    double mant = mpq_get_d_2exp(&exp, mpqref(Q));
    const double log2 = std::log(2.0);
    return std::log(mant) + exp*log2;
  }

  double LogAbs(const BigRat& Q)
  {
    const int s = sign(Q);
    if (s == 0)  CoCoA_THROW_ERROR1(ERR::ReqNonZero);
    if (s > 0)  return log(Q);
    return log(-Q);
  }


  long FloorLog2(const BigRat& Q) { return FloorLogBase(Q,2); }
  long FloorLog10(const BigRat& Q)  { return FloorLogBase(Q,10); }

  long FloorLogBase(const BigRat& Q, const MachineInt& base)
  {
    return FloorLogBase(Q, BigInt(base));
  }

  // BUGLY: Can we merge the two almost identical sections?
  long FloorLogBase(const BigRat& Q, const BigInt& base)
  {
    if (base < 2)  CoCoA_THROW_ERROR2(ERR::BadArg, "base must be at least 2");
    if (IsZero(Q))  CoCoA_THROW_ERROR1(ERR::ReqNonZero);
    if (IsOne(abs(Q)))  return 0;
    const double ApproxLog = LogAbs(Q)/log(base);
    const double delta = 16 * std::abs(ApproxLog) * numeric_limits<double>::epsilon(); // ***ASSUME*** numerical error in ApproxLog is less than delta; 16 is "magic number" (but seems to be large enough)
///???BUG not yet fully implemented    if (ApproxLog > numeric_limits<long>::max())  CoCoA_THROW_ERROR1(ERR::ArgTooBig);
    const long NearestInt = static_cast<long>(std::floor(ApproxLog+0.5)); // probably right but could be too big by 1

    // Easy case if ApproxLog is far from an integer...
    if (std::abs(ApproxLog - NearestInt) > delta)
    {
      const long ans = static_cast<long>(std::floor(ApproxLog));
      CoCoA_ASSERT(abs(Q) >= power(BigRat(base), ans)); // allow for ans < 0
      return ans;
    }
    
    // Hard case: ApproxLog is very close to integer NearestInt
    if (NearestInt >= 0)
    {
      const BigInt pwr = power(base, NearestInt);
      const int test = cmp(abs(Q), pwr);
      if (test >= 0) return NearestInt;
      return NearestInt-1;
    }
    // else NearestInt < 0
    const BigInt pwr = power(base, -NearestInt);
    const BigRat shifted = abs(Q)*pwr;
    const int test = cmp(shifted, 1);
    if (test >= 0) return NearestInt;
    return NearestInt-1;
  }


  bool IsPowerOf2(const BigRat& Q)
  {
    if (sign(Q) <= 0) return false;
    if (IsOneDen(Q)) return IsPowerOf2(num(Q));
    if (IsOneNum(Q)) return IsPowerOf2(den(Q));
    return false;
  }


  BigInt floor(const BigRat& Q)
  {
    BigInt ans;
    mpz_fdiv_q(mpzref(ans), mpq_numref(mpqref(Q)), mpq_denref(mpqref(Q)));
    return ans;
  }

  BigInt ceil(const BigRat& Q)
  {
    BigInt ans;
    mpz_cdiv_q(mpzref(ans), mpq_numref(mpqref(Q)), mpq_denref(mpqref(Q)));
    return ans;
  }


  // This impl clearly guarantees that rounding is compatible with RoundDiv!
  BigInt round(const BigRat& Q)
  {
    if (IsOneDen(Q)) return num(Q);
    return RoundDiv(num(Q), den(Q));
  }


  BigInt CommonDenom(const std::vector<BigRat>& L)
  {
    const long n = len(L);
    vector<BigInt> D(n);
    for (long i=0; i < n; ++i)
      D[i] = den(L[i]);
    return lcm(D);
  }

} // end of namespace CoCoA
