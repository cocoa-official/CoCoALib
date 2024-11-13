//   Copyright (c)  2011,2014  John Abbott and Anna M. Bigatti

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

#include "CoCoA/ToString.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRatOps.H"
#include "CoCoA/FloatApprox.H"
#include "CoCoA/error.H"
#include "CoCoA/utils.H"

#include <algorithm>
using std::fill;
using std::copy;
#include <iostream>
using std::ostream;
#include <limits>
using std::numeric_limits;
//#include <string>
using std::string;
#include <sstream>
using std::ostringstream;

namespace CoCoA
{

  namespace // anonymous for file local fns
  {

    std::string ScientificStr(const MantExp10& ME, long ndigits)
    {
      ostringstream MantissaStream;

      if (IsZero(ME.myMantissa)) MantissaStream << string(ndigits, '0');
      else MantissaStream << ME.myMantissa;
      const string& mantissa = MantissaStream.str();
//       string mantissa; mantissa.reserve(ndigits);
//       if (!IsZero(ME.myMantissa))
//         IsConvertible(mantissa, ME.myMantissa); // cannot fail, except std::bad_alloc
//       else
//         // fill string with ndigits of zeroes...
//         for (long i=0;i<ndigits;++i)mantissa+='0';///???fill(...);

      const int ExpDigits = numeric_limits<long>::digits10;
      string ans; ans.reserve(9+ndigits+ExpDigits); // 9 = size of "-.*10^(-)"

      if (ME.mySign < 0) ans += '-';
      ans += mantissa[0];
      ans += '.';
      //      copy(&mantissa[1], &mantissa[ndigits], back_inserter(ans));
      ans.insert(ans.end(), &mantissa[1], &mantissa[ndigits]);
      ans += "*10^";
      if (ME.myExponent < 0) ans += '(';
      ostringstream exp;
      exp << ME.myExponent;
      ans += exp.str();
      if (ME.myExponent < 0) ans += ')';
      return ans;
    }

    std::string FloatStr(const MantExp10& ME, long ndigits)
    {
      // See doc for info about magic number 8 in line below
      if (ME.myExponent >= ndigits || ME.myExponent < -8) return ScientificStr(ME, ndigits);

      ostringstream MantissaStream;
      if (IsZero(ME.myMantissa)) MantissaStream << string(ndigits, '0');
      else MantissaStream << ME.myMantissa;
      const string& mantissa = MantissaStream.str();
//       string mantissa; mantissa.reserve(ndigits);
//       if (!IsZero(ME.myMantissa))
//         IsConvertible(mantissa, ME.myMantissa); // cannot fail, except std::bad_alloc
//       else
//         // fill string with ndigits of zeroes...
//         for (long i=0;i<ndigits;++i)mantissa+='0';///???fill(...);

      string ans; ans.reserve(3+ndigits); // 3 = size of "-0."
      if (ME.mySign < 0) ans += '-';
      if (ME.myExponent < 0)
      {
        ans += "0.";
        for (long i=-1; i > ME.myExponent; --i)
          ans += '0';
        //        copy(&mantissa[0], &mantissa[ndigits], back_inserter(ans));
        ans.insert(ans.end(), &mantissa[0], &mantissa[ndigits]);
        return ans;
      }

      //   copy(&mantissa[0], &mantissa[ME.myExponent+1], back_inserter(ans));
      ans.insert(ans.end(), &mantissa[0], &mantissa[ME.myExponent+1]);
      ans += '.';
      // copy(&mantissa[ME.myExponent+1], &mantissa[ndigits], back_inserter(ans));
      ans.insert(ans.end(), &mantissa[ME.myExponent+1], &mantissa[ndigits]);
      return ans;
    }

  } // end of anonymous namespace



  std::string ScientificStr(const BigInt& N, const MachineInt& SigFig)
  {
    if (IsNegative(SigFig) || !IsSignedLong(SigFig) || IsZero(SigFig))
      CoCoA_THROW_ERROR(ERR::BadArg, "ScientificStr");
    return ScientificStr(MantissaAndExponent10(N, SigFig), AsSignedLong(SigFig));
  }

  std::string ScientificStr(const BigRat& q, const MachineInt& SigFig)
  {
    if (IsNegative(SigFig) || !IsSignedLong(SigFig) || IsZero(SigFig))
      CoCoA_THROW_ERROR(ERR::BadArg, "ScientificStr");
    return ScientificStr(MantissaAndExponent10(q, SigFig), AsSignedLong(SigFig));
  }

  std::string FloatStr(const BigInt& N, const MachineInt& SigFig)
  {
    if (IsNegative(SigFig) || !IsSignedLong(SigFig) || IsZero(SigFig))
      CoCoA_THROW_ERROR(ERR::BadArg, "FloatStr");
//???    if (abs(N) < power(10,SigFig)) return ToString(N);

    return FloatStr(MantissaAndExponent10(N, SigFig), AsSignedLong(SigFig));
  }

  std::string FloatStr(const BigRat& q, const MachineInt& SigFig)
  {
    if (IsNegative(SigFig) || !IsSignedLong(SigFig) || IsZero(SigFig))
      CoCoA_THROW_ERROR(ERR::BadArg, "FloatStr");
//    if (IsOneDen(q)) return FloatStr(num(q), SigFig);
    return FloatStr(MantissaAndExponent10(q, SigFig), AsSignedLong(SigFig));
  }


  //------------------------------------------------------------------
  // FixedStr

  std::string DecimalStr(const BigInt& N, const MachineInt& DecimalPlaces)
  {
    if (IsNegative(DecimalPlaces) || !IsSignedLong(DecimalPlaces) || IsZero(DecimalPlaces))
      CoCoA_THROW_ERROR(ERR::BadArg, "FixedStr");

    return ToString(N);
  }

  std::string DecimalStr(const BigRat& q, const MachineInt& DecimalPlaces)
  {
    if (IsNegative(DecimalPlaces) || !IsSignedLong(DecimalPlaces) || IsZero(DecimalPlaces))
      CoCoA_THROW_ERROR(ERR::BadArg, "FixedStr");
    if (IsOneDen(q)) return DecimalStr(num(q), DecimalPlaces);
    const long DigitsAfterPoint = AsSignedLong(DecimalPlaces);
    const BigInt N = RoundDiv(abs(num(q))*power(10,DigitsAfterPoint),den(q));
    string digits = ToString(N);
    if (len(digits) < 1+DigitsAfterPoint)
    {
      digits = string(DigitsAfterPoint+1-len(digits), '0') + digits;
    }
    string ans;
    if (q < 0) ans = '-';
    const long IntegerPart = len(digits) - DigitsAfterPoint;
    ans.insert(ans.end(), &digits[0], &digits[IntegerPart]);
    ans += '.';
    ans.insert(ans.end(), &digits[IntegerPart], &digits[IntegerPart+DigitsAfterPoint]);
    return ans;
  }

} // end of namespace CoCoA
