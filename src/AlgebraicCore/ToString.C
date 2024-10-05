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


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/ToString.C,v 1.6 2022/02/18 14:12:01 abbott Exp $
// $Log: ToString.C,v $
// Revision 1.6  2022/02/18 14:12:01  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.5  2021/01/07 15:16:53  abbott
// Summary: Corrected copyright
//
// Revision 1.4  2020/06/17 15:49:29  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.3  2018/05/22 14:16:40  abbott
// Summary: Split BigRat into BigRat (class defn + ctors) and BigRatOps
//
// Revision 1.2  2018/05/18 12:24:47  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.1  2014/04/11 15:05:49  abbott
// Summary: Fns for converting BigInt/BigRat to string (used to be in decimal)
// Author: JAA
//
// Revision 1.12  2014/04/10 15:31:36  abbott
// Summary: Replaced MantExp by MantExp10; new fns FixedStr and decimals (temporary)
// Author: JAA
//
// Revision 1.11  2013/03/26 14:56:06  abbott
// Updated the conversion fns (in ptic removed procedure "convert");
// numerous consequential changes.
//
// Revision 1.10  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.9  2012/02/08 17:07:33  bigatti
// -- definition of global in .C file
//
// Revision 1.8  2012/01/26 16:57:42  abbott
// Moved value of ourDefaultSigFig from H file to C file
// (this allows clean compilation on MSVC, also avoids
// excessive recompilation if we change the value).
//
// Revision 1.7  2012/01/26 16:47:00  bigatti
// -- changed back_inserter into insert
//
// Revision 1.6  2011/11/09 14:29:37  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.5  2011/08/24 10:32:04  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.4  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.3  2011/08/03 09:05:05  abbott
// Added new fn DecimalStr.
//
// Revision 1.2  2011/08/02 13:28:46  abbott
// Cleaned up impl.  Now MantissaAndExponent has default SigFig too.
// Used naming convention for members of MantExp.
//
// Revision 1.1  2011/03/01 15:21:22  abbott
// Added new functions: MantissaAndExponent, FloatStr (in CoCoALib).
//
//
