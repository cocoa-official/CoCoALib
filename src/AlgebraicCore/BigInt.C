//   Copyright (c)  2003-2011  John Abbott and Anna M. Bigatti

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

#include "CoCoA/BigInt.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/error.H"
#include "CoCoA/utils.H"


#include <iostream>
using std::ostream;
using std::istream;
using std::ios;  // obsolete???

#include <string>
using std::string;

#include <vector>
using std::vector;

namespace CoCoA
{

  BigInt::BigInt()
  {
    mpz_init(myRepr);
  }


  BigInt::BigInt(const MachineInt& n)
  {
    if (IsNegative(n))
      mpz_init_set_si(myRepr, AsSignedLong(n));
    else
      mpz_init_set_ui(myRepr, AsUnsignedLong(n));
  }


  BigInt::BigInt(const std::string& str, ReadFromString /*NotUsed*/)
  {
    mpz_init(myRepr);
//     if (base != 0 && (base < 2 || base > 36))
//       CoCoA_THROW_ERROR1(ERR::BadNumBase);
    if (mpz_set_str(myRepr, str.c_str(), 10) != 0)
      CoCoA_THROW_ERROR2(ERR::BadArg, "BigIntFromString");
  }


  BigInt::BigInt(const mpz_t N, CopyFromMPZ /*NotUsed*/)
  {
    if (N == nullptr)
      CoCoA_THROW_ERROR2(ERR::NullPtr, "BigIntFromMPZ");

    mpz_init_set(myRepr, N);
  }


// COPY CONSTRUCTOR

  BigInt::BigInt(const BigInt& from)
  {
    mpz_init_set(myRepr, from.myRepr);
  }


  BigInt::BigInt(BigInt&& from) /*noexcept*/ // std move ctor
  {
////std::clog<<"BigInt MOVE"<<std::endl;
    mpz_init(myRepr);
    mpz_swap(myRepr, from.myRepr);
  }


  BigInt& BigInt::operator=(const BigInt& rhs)
  {
    if (this == &rhs) return *this;
    mpz_set(myRepr, rhs.myRepr);
    return *this;
  }

  BigInt& BigInt::operator=(BigInt&& rhs)
  {
//????    if (this == &rhs) return *this;
    mpz_swap(myRepr, rhs.myRepr);
    return *this;
  }



  namespace // anonymous
  {

    // BUG/SLUG slow and wastes space -- better to use ostringstream?
    std::string ConvertToString(const BigInt& src, int base)
    {
      if (base < 2 || base > 36)
        CoCoA_THROW_ERROR1(ERR::BadNumBase);
      const size_t digits = mpz_sizeinbase(mpzref(src), base);
      string ans; ans.reserve(digits+1); // +1 to allow for minus sign (overflow check?)
      vector<char> buffer(digits+2); // +2 to allow for minus sign and terminating NUL
      mpz_get_str(&buffer[0], base, mpzref(src));
      ans = &buffer[0]; // Won't throw as space is already reserved.
      return ans;
    }

  } // end of namespace anonymous


  ostream& operator<<(ostream& out, const BigInt& N)
  {
    if (!out) return out;  // short-cut for bad ostreams
    CoCoA_ASSERT(IsDecimal(out));
//???redmine 1547    if (!IsDecimal(out)) CoCoA_THROW_ERROR1("ostream is not in \"decimal\" mode");
//    Using GMPXX next two lines should do it all...
//    mpz_class Ncopy(N.myRepr);
//    return out << Ncopy;

    // Dispose of the exceptional case n=0 immediately.
    if (IsZero(N)) return out << '0';
///???    const int base = 10;
    const int base = (out.flags() & ios::oct) ? 8 :
                     (out.flags() & ios::hex) ? 16 : 10;

    // prefix will contain sign and/or base specification
    string prefix;

    // Firstly put the sign in the prefix...
    if (N < 0) prefix += '-';
    else if (out.flags() & ios::showpos) prefix += '+';

    // Next put base indication in the prefix if required...
    if (out.flags() & ios::showbase)
    {
      if (base == 8 || base == 16) prefix += '0';
      if (base == 16) prefix += 'x';
    }

    out << prefix << ConvertToString(abs(N), base);
    return out;
  }


  // // Returns true if c is a valid digit for the given base
  // // (both upper and lower cases are allowed for hexadecimal).
  // // Note: base must be 8, 10 or 16 (otherwise always gives false)
  // bool IsDigitBase(char c, int base) // base = 8, 10 or 16
  // {
  //   if (!(base == 10 || base == 8 || base == 16))  CoCoA_THROW_ERROR2(ERR::BadArg, "Base must be 10, 8 or 16");
  //   switch (base)
  //   {
  //   case 10:
  //     return isdigit(c);
  //   case 16:
  //     return isxdigit(c);
  //   case 8:
  //     return (isdigit(c) &&  c != '8' && c != '9');
  //   default:
  //     return false; // should never get here
  //   }
  // }


  // This fn does not "expect" leading whitespace!
  // Leaves "in" in good state; if "in" was already not good, an exc is thrown
  std::string ScanUnsignedIntegerLiteral(std::istream& in)
  {
    if (!in.good())
      CoCoA_THROW_ERROR1("istream is not good");
    if (!IsDecimal(in))
      CoCoA_THROW_ERROR1("istream is not in \"decimal\" mode");

    // Read in as many digits as possible.
    string digits;
    while (true)
    {
      char ch;
      in.get(ch);
      if (in.eof()) { in.clear(); break; }
      if (!isdigit(ch)) { in.unget(); break; }
      digits += ch;
    }
    return digits; // BUG???: better to use  std::move???
  }

  
  std::istream& operator>>(std::istream& in, BigInt& N)
  {
    if (!in.good())  CoCoA_THROW_ERROR1("istream is not good");
    if (!IsDecimal(in))  CoCoA_THROW_ERROR1("istream is not in \"decimal\" mode");
    ws(in);

    // Look for sign of number.
    const char FirstChar = in.peek(); // this may set eofbit
    if (in.eof())  CoCoA_THROW_ERROR1("EOF");
    if (FirstChar == '-' || FirstChar == '+')  in.ignore();
    const string digits = ScanUnsignedIntegerLiteral(in); // leaves "in" in good state
    if (digits.empty())
      CoCoA_THROW_ERROR1("No decimal digits in input");
    // We found some digits, so convert them into a number.
    N = BigIntFromString(digits); // could throw if number is huge.
    if (FirstChar == '-')  negate(N);
    return in;
  }


  OpenMathOutput& operator<<(OpenMathOutput& OMOut, const BigInt& N)
  {
    OMOut->mySend(N);
    return OMOut;
  }


  OpenMathInput& operator>>(OpenMathInput& OMIn, BigInt& /*N*/)
  {
    CoCoA_THROW_ERROR1(ERR::NYI);
    return OMIn;
  }


  //---------------------------------------------------------------------------
  // Assignment and assignment arithmetic functions

  BigInt& BigInt::operator+=(const BigInt& rhs)
  {
    mpz_add(myRepr, myRepr, rhs.myRepr);
    return *this;
  }


  BigInt& BigInt::operator-=(const BigInt& rhs)
  {
    mpz_sub(myRepr, myRepr, rhs.myRepr);
    return *this;
  }


  BigInt& BigInt::operator*=(const BigInt& rhs)
  {
    mpz_mul(myRepr, myRepr, rhs.myRepr);
    return *this;
  }


  BigInt& BigInt::operator/=(const BigInt& rhs)
  {
    if (IsZero(rhs))
      CoCoA_THROW_ERROR1(ERR::DivByZero);
//???    if (rhs < 0 && !mpz_divisible_p(myRepr, mpzref(rhs)))
//???      CoCoA_THROW_ERROR1(ERR::IntDivByNeg);
    mpz_tdiv_q(myRepr, myRepr, rhs.myRepr);
    return *this;
  }


  BigInt& BigInt::operator%=(const BigInt& rhs)
  {
    if (IsZero(rhs))
      CoCoA_THROW_ERROR1(ERR::ReqNonZeroModulus);
//???    if (rhs < 0 && !mpz_divisible_p(myRepr, mpzref(rhs)))
//???      CoCoA_THROW_ERROR1(ERR::IntDivByNeg);
    mpz_tdiv_r(myRepr, myRepr, rhs.myRepr);
    return *this;
  }


  //---------------------------------------------------------------------------
  // Assignment and assignment arithmetic with rhs a MachineInt

  BigInt& BigInt::operator= (const MachineInt& rhs)
  {
    if (IsNegative(rhs))
      mpz_set_si(myRepr, AsSignedLong(rhs));
    else
      mpz_set_ui(myRepr, AsUnsignedLong(rhs));

    return *this;
  }


  BigInt& BigInt::operator+=(const MachineInt& rhs)
  {
    if (IsNegative(rhs))
      mpz_sub_ui(myRepr, myRepr, uabs(rhs));
    else
      mpz_add_ui(myRepr, myRepr, AsUnsignedLong(rhs));
    return *this;
  }


  BigInt& BigInt::operator-=(const MachineInt& rhs)
  {
    if (IsNegative(rhs))
      mpz_add_ui(myRepr, myRepr, uabs(rhs));
    else
      mpz_sub_ui(myRepr, myRepr, AsUnsignedLong(rhs));
    return *this;
  }


  BigInt& BigInt::operator*=(const MachineInt& rhs)
  {
    if (IsNegative(rhs))
      mpz_mul_si(myRepr, myRepr, AsSignedLong(rhs));
    else
      mpz_mul_ui(myRepr, myRepr, AsUnsignedLong(rhs));
    return *this;
  }


  BigInt& BigInt::operator/=(const MachineInt& rhs)
  {
    if (IsZero(rhs))
      CoCoA_THROW_ERROR1(ERR::DivByZero);
    if (IsNegative(rhs))
      mpz_neg(myRepr, myRepr);

    mpz_tdiv_q_ui(myRepr, myRepr, uabs(rhs));
    return *this;
  }


  BigInt& BigInt::operator%=(const MachineInt& rhs)
  {
    if (IsZero(rhs))
      CoCoA_THROW_ERROR1(ERR::ReqNonZeroModulus);
    mpz_tdiv_r_ui(myRepr, myRepr, uabs(rhs));
    return *this;
  }



  //---------------------------------------------------------------------------
  // increment and decrement

  BigInt& BigInt::operator++()
  {
    mpz_add_ui(myRepr, myRepr, 1);
    return *this;
  }


  BigInt& BigInt::operator--()
  {
    mpz_sub_ui(myRepr, myRepr, 1);
    return *this;
  }


  const BigInt BigInt::operator++(int)  // INEFFICIENT
  {
    BigInt ans(*this);
    ++*this;
    return ans;
  }


  const BigInt BigInt::operator--(int)  // INEFFICIENT
  {
    BigInt ans(*this);
    --*this;
    return ans;
  }


  //---------------------------------------------------------------------------

} // end of namespace CoCoA
