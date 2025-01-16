#ifndef CoCoA_BigRat_H
#define CoCoA_BigRat_H

//   Copyright (c)  2009-2010, 2018  John Abbott,  Anna M. Bigatti

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


#include "gmp.h"

#include <iosfwd>
//using std::ostream;
//using std::istream;

namespace CoCoA
{

  class BigInt;          // fwd decl -- defined in BigInt.H
  class MachineInt;      // fwd decl -- defined in MachineInt.H
  class MantExp10;       // fwd decl -- defined in FloatApprox.H
  class MantExp2;        // fwd decl -- defined in FloatApprox.H
  class OpenMathOutput;  // fwd decl -- defined in OpenMath.H
  class OpenMathInput;   // fwd decl -- defined in OpenMath.H

  class BigRat
  {
  private:
    friend mpq_t& mpqref(BigRat& N);
    friend const mpq_t& mpqref(const BigRat& N);
    mpq_t myRepr;

  public:
    BigRat();
    enum class ReduceFlag { NotReduced, AlreadyReduced };
    static constexpr ReduceFlag NotReduced = ReduceFlag::NotReduced;
    static constexpr ReduceFlag AlreadyReduced = ReduceFlag::AlreadyReduced;
    explicit BigRat(const MachineInt& n);
    explicit BigRat(const BigInt& N);
    BigRat(const MachineInt& n1, const MachineInt& n2, ReduceFlag status=NotReduced);  ///< creates the fraction n1/n2
    BigRat(const MachineInt& n1, const BigInt& N2, ReduceFlag status=NotReduced);  ///< creates the fraction n1/N2
    BigRat(const BigInt& N1, const MachineInt& n2, ReduceFlag status=NotReduced);  ///< creates the fraction N1/n2
    BigRat(const BigInt& N1, const BigInt& N2, ReduceFlag status=NotReduced);  ///< creates the fraction N1/N2
    explicit BigRat(const MantExp2& ME2);
    explicit BigRat(const MantExp10& ME10);
    enum class OneOverZero_t { OneOverZero };
    static constexpr OneOverZero_t OneOverZero = OneOverZero_t::OneOverZero;
    explicit BigRat(OneOverZero_t);
    BigRat(const BigRat& from); // std copy ctor
    BigRat(BigRat&& from) /*noexcept*/; // std move ctor (cannot be noexcept because GMP allocs)
    ~BigRat();

  private: // pseudo-ctor
    enum class CopyFromMPQ { marker };
    enum class ReadFromString { marker };
    explicit BigRat(const mpq_t q, CopyFromMPQ); // copies the mpq_t value
    explicit BigRat(const std::string& str, ReadFromString, ReduceFlag status=NotReduced);
    friend BigRat BigRatFromMPQ(const mpq_t q) { return BigRat(q, BigRat::CopyFromMPQ::marker); }
    friend BigRat BigRatFromString(const std::string& str, ReduceFlag ReducedStatus/*=NotReduced*/) { return BigRat(str, BigRat::ReadFromString::marker, ReducedStatus); }

  private:
    void myAssign(const BigInt& N1, const BigInt& N2, ReduceFlag status=NotReduced); ///< does the work for the ctors
  public:

    BigRat& operator=(const BigRat& rhs);
    BigRat& operator=(BigRat&& rhs);

    // -------- functions that modify at least one argument or `*this' ----------

    BigRat& operator+=(const BigRat& rhs);
    BigRat& operator-=(const BigRat& rhs);
    BigRat& operator*=(const BigRat& rhs);
    BigRat& operator/=(const BigRat& rhs);
  			
    // Same but with RHS a BigInt...
    BigRat& operator= (const BigInt& rhs);
    BigRat& operator+=(const BigInt& rhs);
    BigRat& operator-=(const BigInt& rhs);
    BigRat& operator*=(const BigInt& rhs);
    BigRat& operator/=(const BigInt& rhs);
  			
    // Same but with RHS a MachineInt...
    BigRat& operator= (const MachineInt& rhs);
    BigRat& operator+=(const MachineInt& rhs);
    BigRat& operator-=(const MachineInt& rhs);
    BigRat& operator*=(const MachineInt& rhs);
    BigRat& operator/=(const MachineInt& rhs);


    const BigRat& operator++( );
    const BigRat& operator--( );
    const BigRat operator++(int); // INEFFICIENT
    const BigRat operator--(int); // INEFFICIENT

  };


  // pseudo-ctors (impl as friends)
  BigRat BigRatFromMPQ(const mpq_t q);
  BigRat BigRatFromString(const std::string& str, BigRat::ReduceFlag ReducedStatus=BigRat::NotReduced);

  // I/O FUNCTIONS

  std::ostream& operator<<(std::ostream& out, const BigRat& Q);
  std::istream& operator>>(std::istream& in, BigRat& Q);   // throws or leaves "in" in good state

  OpenMathOutput& operator<<(OpenMathOutput& OMOut, const BigRat& Q);
  OpenMathInput& operator>>(OpenMathInput& OMIn, BigRat& Q);

  bool IsZero(const BigRat& Q);
  bool IsOne(const BigRat& Q);
  bool IsMinusOne(const BigRat& Q);
  bool IsOneNum(const BigRat& Q);
  bool IsOneDen(const BigRat& Q);
  int sign(const BigRat& Q);


  void swap(BigRat& a, BigRat& b);
  BigInt num(const BigRat& Q);
  BigInt den(const BigRat& Q);


  // Inline fns

  inline mpq_t& mpqref(BigRat& Q)
  {
    return Q.myRepr;
  }

  inline const mpq_t& mpqref(const BigRat& Q)
  {
    return Q.myRepr;
  }


  inline bool IsZero(const BigRat& Q)
  {
    return (mpq_sgn(mpqref(Q)) == 0);
//    return IsZero(num(Q));
  }


  inline bool IsOne(const BigRat& Q)
  {
    return mpq_cmp_ui(mpqref(Q), 1,1) == 0;
  }


  inline bool IsMinusOne(const BigRat& Q)
  {
    return mpq_cmp_si(mpqref(Q), -1,1) == 0;
  }


  inline bool IsOneNum(const BigRat& Q)
  {
    return mpz_cmp_ui(mpq_numref(mpqref(Q)), 1) == 0;
  }


  inline bool IsOneDen(const BigRat& Q)
  {
    return mpz_cmp_ui(mpq_denref(mpqref(Q)), 1) == 0;
  }


  inline int sign(const BigRat& Q)
  {
    return mpq_sgn(mpqref(Q));
  }


} // end of namespace CoCoA

#endif
