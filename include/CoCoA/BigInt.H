#ifndef CoCoA_BigInt_H
#define CoCoA_BigInt_H

//   Copyright (c)  2003-2011  John Abbott,  Anna M. Bigatti

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


#include <iosfwd> // already included in MachineInt.H
//using std::ostream;
//using std::istream;
#include <string>

#include "gmp.h"

namespace CoCoA
{

  class MachineInt;      // fwd decl -- defined in MachineInt.H
  class OpenMathOutput;  // fwd decl -- defined in OpenMath.H
  class OpenMathInput;   // fwd decl -- defined in OpenMath.H

  class BigInt
  {
  private:  // data member
    friend mpz_t& mpzref(BigInt& N) noexcept;
    friend const mpz_t& mpzref(const BigInt& N) noexcept;
    mpz_t myRepr;

  public:
    BigInt();
    explicit BigInt(const MachineInt& n);
    // The following ctors cause ambiguity problems: see redmine 586 comment 13
    // BigInt(signed char n);
    // BigInt(unsigned char n);
    // BigInt(short n);
    // BigInt(unsigned short n);
    // BigInt(int n);
    // BigInt(unsigned int n);
    // BigInt(long n);
    // BigInt(unsigned long n);
    BigInt(const BigInt& from);  // std copy ctor
    BigInt(BigInt&& from) /*noexcept*/;  // std move ctor (cannot be noexcept because allocs mem)
    ~BigInt();

  private: // pseudo-ctors
    enum class CopyFromMPZ { marker };
    enum class ReadFromString { marker };
    BigInt(const std::string& str, ReadFromString);
    BigInt(const mpz_t N, CopyFromMPZ); // copies the mpz_t value!
    friend BigInt BigIntFromMPZ(const mpz_t N) { return BigInt(N, BigInt::CopyFromMPZ::marker); }
    friend BigInt BigIntFromString(const std::string& str) { return BigInt(str, BigInt::ReadFromString::marker); }
  public:

    BigInt& operator=(const BigInt& rhs);
    BigInt& operator=(BigInt&& rhs);

    // -------- functions that modify at least one argument or `*this' ----------

    BigInt& operator+=(const BigInt& rhs);
    BigInt& operator-=(const BigInt& rhs);
    BigInt& operator*=(const BigInt& rhs);
    BigInt& operator/=(const BigInt& rhs);
    BigInt& operator%=(const BigInt& rhs);
  			
    // Same but with RHS a MachineInt...
    BigInt& operator= (const MachineInt& rhs);
    BigInt& operator+=(const MachineInt& rhs);
    BigInt& operator-=(const MachineInt& rhs);
    BigInt& operator*=(const MachineInt& rhs);
    BigInt& operator/=(const MachineInt& rhs);
    BigInt& operator%=(const MachineInt& rhs);


    BigInt& operator++( );
    BigInt& operator--( );
    const BigInt operator++(int); // INEFFICIENT
    const BigInt operator--(int); // INEFFICIENT

  };


  // pseudo-ctors (impl as friends)
  BigInt BigIntFromMPZ(const mpz_t N);
  BigInt BigIntFromString(const std::string& str);

  // I/O FUNCTIONS

  std::ostream& operator<<(std::ostream& out, const BigInt& N);
  std::istream& operator>>(std::istream& in, BigInt& N);   // throws or leaves "in" in good state
////  bool IsDigitBase(char c, int base); // only for base = 8, 10 or 16
  std::string ScanUnsignedIntegerLiteral(std::istream& in);  // throws or leaves "in" in good state

  OpenMathOutput& operator<<(OpenMathOutput& OMOut, const BigInt& N);
  OpenMathInput& operator>>(OpenMathInput& OMIn, BigInt& N);


  void swap(BigInt& a, BigInt& b) noexcept;

  /////////////////////////////////////////////////////////////////////////////
  // Inline functions

  inline mpz_t& mpzref(BigInt& N) noexcept
  {
    return N.myRepr;
  }

  inline const mpz_t& mpzref(const BigInt& N) noexcept
  {
    return N.myRepr;
  }


  inline BigInt::~BigInt()
  {
    mpz_clear(myRepr);
  }


  inline void swap(BigInt& a, BigInt& b) noexcept
  {
    if (&a == &b) return;
    mpz_swap(mpzref(a), mpzref(b));
  }


} // end of namespace CoCoA

#endif
