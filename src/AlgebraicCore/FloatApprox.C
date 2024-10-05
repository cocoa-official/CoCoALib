//   Copyright (c)  2014  John Abbott and Anna M. Bigatti

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

#include "CoCoA/FloatApprox.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRatOps.H"
#include "CoCoA/error.H"


#include <iostream>
using std::ostream;

namespace CoCoA
{

  const int MantExp2::ourDefaultMantBits = 53; // same as IEEE "double"

  MantExp2::MantExp2(int s, long e, const BigInt& m, long NumDigits):
      mySign(s),
      myExponent(e),
      myMantissa(m),
      myNumDigits(NumDigits)
  {
    CoCoA_ASSERT((s==0 && e == 0 && IsZero(m) && NumDigits==0) ||
                 (s*s==1 && m>0 && FloorLog2(m)==NumDigits-1));
  }


  MantExp2 MantissaAndExponent2(const MachineInt& n, const MachineInt& MantWidth)
  {
    return MantissaAndExponent2(BigRat(n,1), MantWidth);
  }


  // let rational version do the work so that halves are rounded consistently!
  MantExp2 MantissaAndExponent2(const BigInt& N, const MachineInt& MantWidth)
  {
    return MantissaAndExponent2(BigRat(N,1), MantWidth);
  }


  // Simple/compact rather than fast;  is speed so important here?
  MantExp2 MantissaAndExponent2(const BigRat& q, const MachineInt& MantWidth)
  {
    if (IsNegative(MantWidth) || !IsSignedLong(MantWidth) || AsSignedLong(MantWidth) < 2)
      CoCoA_THROW_ERROR(ERR::BadArg, "MantissaAndExponent2");
    if (IsZero(q)) return MantExp2(0,0,BigInt(0),0);
    const long MantBits = AsSignedLong(MantWidth);
    if (MantBits > 1000000000L) CoCoA_THROW_ERROR(ERR::ArgTooBig, "Precision too high");
    BigInt N = abs(num(q));
    BigInt D = den(q);
    const int SignQ = sign(q);
    const long LogQ = FloorLog2(q);
    const long exp = LogQ-MantBits+1;  // NB exp-1, exp+1, and -exp  will not overflow!
    if (exp <= 0)
      mpz_mul_2exp(mpzref(N), mpzref(N), -exp); // N *= 2^|exp|
    else
      mpz_mul_2exp(mpzref(D), mpzref(D), exp);  // D *= 2^exp

    N = RoundDiv(N,D);
    if (FloorLog2(N) == MantBits) // true iff mantissa has "overflowed"
      return MantExp2(SignQ, 1+LogQ, N/2, MantBits);
    return MantExp2(SignQ, LogQ, N, MantBits);
  }


  std::ostream& operator<<(std::ostream& out, const MantExp2& ME)
  {
    if (!out) return out;  // short-cut for bad ostreams
    out << "MantExp2(sign=" << ME.mySign << ", exp=" << ME.myExponent << ", mant=" << ME.myMantissa << ", NumDigits=" << ME.myNumDigits << ")";
    return out;
  }


  //------------------------------------------------------------------
  BigRat FloatApprox(const MachineInt& n, const MachineInt& MantBits)
  {
    return BigRat(MantissaAndExponent2(n, MantBits));
  }

  BigRat FloatApprox(const BigInt& N, const MachineInt& MantBits)
  {
    return BigRat(MantissaAndExponent2(N, MantBits));
  }

  BigRat FloatApprox(const BigRat& q, const MachineInt& MantBits)
  {
    return BigRat(MantissaAndExponent2(q, MantBits));
  }


  //------------------------------------------------------------------
  // Decimal "floating point" representation

  const int MantExp10::ourDefaultSigFig = 5;


  MantExp10::MantExp10(int s, long e, const BigInt& m, long NumDigits):
      mySign(s),
      myExponent(e),
      myMantissa(m),
      myNumDigits(NumDigits)
  {
    CoCoA_ASSERT((s==0 && e == 0 && IsZero(m) && NumDigits==0) ||
                 (s*s==1 && m>0 && FloorLog10(m)==NumDigits-1));
  }


  std::ostream& operator<<(std::ostream& out, const MantExp10& ME)
  {
    if (!out) return out;  // short-cut for bad ostreams
    out << "MantExp10(sign=" << ME.mySign << ", exp=" << ME.myExponent << ", mant=" << ME.myMantissa << ", NumDigits=" << ME.myNumDigits << ")";
    return out;
  }


  // Added flag to disable overflow checking when calling power: since the args to power
  // create a number small than one handed in as arg.  Concievbly there could be a
  // problem if SigFig is very large -- should it be limited? (maybe 10^9?)
  MantExp10 MantissaAndExponent10(const BigInt& N, const MachineInt& SigFig)
  {
    if (IsNegative(SigFig) || !IsSignedLong(SigFig) || IsZero(SigFig))
      CoCoA_THROW_ERROR(ERR::BadArg, "MantissaAndExponent10");
    if (IsZero(N)) return MantExp10(0,0,BigInt(0),0);
    const long ndigits = AsSignedLong(SigFig);
    const int SignN = sign(N);
    const long e = FloorLog10(N); // overflow???
    if (e < ndigits)
      return MantExp10(SignN, e, abs(N)*power(10, ndigits-e-1, BigIntPowerOverflowCheck::DISABLED), ndigits);
    const BigInt HalfULP = 5*power(10, e-ndigits, BigIntPowerOverflowCheck::DISABLED);
    const BigInt digits = (1+abs(N)/HalfULP)/2;
    // Must check whether digits has overflowed...
    if (abs(digits) == power(10, ndigits, BigIntPowerOverflowCheck::DISABLED))  // overflow check or not check??
      return MantExp10(SignN, e+1, digits/10, ndigits);
    return MantExp10(SignN, e, digits, ndigits);
  }


  MantExp10 MantissaAndExponent10(const BigRat& q, const MachineInt& SigFig)
  {
    if (IsNegative(SigFig) || !IsSignedLong(SigFig) || IsZero(SigFig))
      CoCoA_THROW_ERROR(ERR::BadArg, "MantissaAndExponent10");
    if (IsZero(q)) return MantExp10(0,0,BigInt(0),0);
    if (IsOneDen(q)) return MantissaAndExponent10(num(q), SigFig);
    const long ndigits = AsSignedLong(SigFig);
    const int signq = sign(q);
    const long e = FloorLog10(q); // overflow???
    BigInt digits;
    if (e < ndigits)
      digits = round(abs(q)*power(10,ndigits-e-1, BigIntPowerOverflowCheck::DISABLED));
    else
      digits = round(abs(q)/power(10,1+e-ndigits, BigIntPowerOverflowCheck::DISABLED)); 
    // Must check whether digits has overflowed...
    if (abs(digits) == power(10, ndigits, BigIntPowerOverflowCheck::DISABLED))  // overflow check or not check??
      return MantExp10(signq, e+1, digits/10, ndigits);
    return MantExp10(signq, e, digits, ndigits);
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/FloatApprox.C,v 1.17 2023/12/01 20:25:28 abbott Exp $
// $Log: FloatApprox.C,v $
// Revision 1.17  2023/12/01 20:25:28  abbott
// Summary: Added flag to disable overflow check in power (for integers)
//
// Revision 1.16  2022/02/18 14:11:54  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.15  2021/11/19 16:32:45  abbott
// Summary: Check for excessive precision
//
// Revision 1.14  2021/01/07 15:07:03  abbott
// Summary: Corrected copyright
//
// Revision 1.13  2020/06/17 15:49:23  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.12  2018/05/22 14:16:39  abbott
// Summary: Split BigRat into BigRat (class defn + ctors) and BigRatOps
//
// Revision 1.11  2018/05/18 12:13:36  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.10  2016/11/11 14:15:32  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.9  2016/03/25 20:03:56  abbott
// Summary: New simplified impls of FloatApprox (using new BigRat ctor from MantExp structures)
//
// Revision 1.8  2015/11/23 18:21:34  abbott
// Summary: Renamed ILogBase -> FloorLogBase; added FloorLog2, FloorLog10
//
// Revision 1.7  2015/10/09 18:26:36  abbott
// Summary: Corrected redmine reference
//
// Revision 1.6  2015/10/09 18:18:27  abbott
// Summary: Renamed "abs" to "uabs" for MachineInt; new fn "negate"; see redmine 783
//
// Revision 1.5  2014/05/14 13:18:11  abbott
// Summary: Updated impls of FloatApprox to follow new defn of MantissaAndExponent2
// Author: JAA
//
// Revision 1.4  2014/05/14 10:51:16  abbott
// Summary: Added new field myNumDigits to MantExp2 and MantExp10
// Author: JAA
//
// Revision 1.3  2014/05/13 11:13:02  abbott
// Summary: MantissaAndExponent2 now accepts NumBits from 2 onwards
// Author: JAA
//
// Revision 1.2  2014/04/11 13:33:03  abbott
// Summary: Added MantissaAndExponent2 and MantissaAndExponent10
// Author: JAA
//
// Revision 1.1  2014/04/10 15:32:10  abbott
// Summary: New fn FloatApprox
// Author: JAA
//
//
