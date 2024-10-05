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



// The next few lines contain RCS header/log information.
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/BigInt.H,v 1.24 2022/10/10 09:46:23 abbott Exp $
// $Log: BigInt.H,v $
// Revision 1.24  2022/10/10 09:46:23  abbott
// Summary: Added move-assignment op
//
// Revision 1.23  2022/02/18 14:11:44  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.22  2021/09/11 15:17:54  abbott
// Summary: Added comment
//
// Revision 1.21  2021/03/03 22:09:32  abbott
// Summary: New enum class (redmine 894)
//
// Revision 1.20  2021/01/31 10:07:41  abbott
// Summary: Added noexcept; removed IsDigitBase
//
// Revision 1.19  2021/01/07 14:09:37  abbott
// Summary: Corrected copyright
//
// Revision 1.18  2020/12/04 10:33:06  abbott
// Summary: Revised according to redmine 1529
//
// Revision 1.17  2020/07/28 08:00:51  abbott
// Summary: Added move ctor
//
// Revision 1.16  2018/05/17 15:26:55  bigatti
// -- removed doxygen initial comment
//
// Revision 1.15  2018/04/20 18:51:25  abbott
// Summary: Changed ctors for BigInt/BigRat from string or from MPZ/MPQ
//
// Revision 1.14  2016/10/08 19:45:04  abbott
// Summary: Exposed "new" fn ScanUnsignedIntegerLiteral (rtns a string)
//
// Revision 1.13  2016/07/21 14:13:49  abbott
// Summary: Added new fn ScanUnsignedIntegerLiteral
//
// Revision 1.12  2015/10/09 18:17:24  abbott
// Summary: Added commented out new ctors direct from machine integer types. (redmine 586)
//
// Revision 1.11  2014/07/15 12:21:11  abbott
// Summary: Removed some cruft
// Author: JAA
//
// Revision 1.10  2014/07/15 11:55:19  abbott
// Summary: Removed "explicit" from ctor for BigInt from MachineInt
// Author: JAA
//
// Revision 1.9  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.8  2012/05/25 13:01:23  abbott
// Added fn IsDivisible.
//
// Revision 1.7  2012/02/08 13:16:45  bigatti
// - fixed comment (... starting with -- ...)
//
// Revision 1.6  2012/02/02 15:55:19  abbott
// Replaced include of MachineInt by fwd decl; reactivated include of <iosfwd>
//
// Revision 1.5  2011/11/09 13:38:31  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.4  2011/08/23 14:29:38  bigatti
// -- added log for old name "ZZ.H"
//
// Revision 1.3  2011/08/14 15:52:18  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.2  2011/08/12 16:31:03  abbott
// COMMENTED OUT SOME FNS SO THAT BigInt CAN EXIST ALONGSIDE ZZ
// FOR THE TIME BEING.
//
// Revision 1.1  2011/08/12 15:21:26  abbott
// Added BigInt impl (derived from ZZ); not used by anyone yet.
//
// *************** log for old name: "ZZ.H" ***************
// Revision 1.23  2011/03/21 11:13:25  abbott
// Changed return type for operator%(ZZ,MachineInt) from
// unsigned long to long to fit in with the new coding conventions.
//
// Revision 1.22  2011/03/16 15:32:00  abbott
// Changed arg type "unsigned long" into "const MachineInt&"
// (in commented out code).
//
// Revision 1.21  2011/03/14 10:32:22  abbott
// Changed size_t into long (return types of NumDigits & exponent).
//
// Revision 1.20  2011/01/14 17:21:15  abbott
// Added isqrt, iroot, IsExactIroot, IsSquare, IsPower.
//
// Revision 1.19  2010/12/26 13:03:16  abbott
// Added ILogBase function (to ZZ & QQ).
//
// Revision 1.18  2010/03/22 11:50:31  abbott
// Added ctor from a string.
// Fixed stupid bug in operator-.
//
// Revision 1.17  2010/03/18 13:54:20  abbott
// Added openmath output fns (moved from OpenMath files).
//
// Revision 1.16  2010/03/05 21:35:48  abbott
// Removed unnecessary include of assert.H
//
// Revision 1.15  2010/02/16 10:22:21  abbott
// Changed include directive for gmp.h -- uses quotes instead of <...>
//
// Revision 1.14  2009/12/29 22:44:32  abbott
// Removed buggy proxy class ZZ::rtn.
// Consequent changes for function prototypes also in NumTheory.
// Corrected some minor buglets in NumTheory.
//
// Revision 1.13  2009/11/26 16:18:00  bigatti
// -- including string ZZ.C instead of ZZ.H
//
// Revision 1.12  2009/10/08 13:41:19  abbott
// Renamed "round" into "RoundDiv".
// Added two new signatures for "RoundDiv".
//
// Revision 1.11  2009/06/11 14:05:29  abbott
// CLEANING: Removed several functions which are now gathered in NumTheory.H/C
//           (for example: gcd, lcm, PowerMod, InvMod).
//
// Revision 1.10  2009/06/05 12:08:28  abbott
// Changed return type of operator%(ZZ,MachineInt); it is now unsigned long
// instead of ZZ.
//
// Revision 1.9  2008/12/16 21:14:47  abbott
// In functions taking a machine integer changed arg type from MachineInt
// to const-ref-MachineInt.
//
// Revision 1.8  2008/11/18 10:25:48  abbott
// Added function round.
//
// Revision 1.7  2008/04/22 13:09:16  abbott
// Removed IsPositive and IsNegative functions for ZZs.
// Comparison between RingElem and 0 is now handled specially (specially fast?).
//
// Revision 1.6  2007/10/30 17:14:11  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.5  2007/05/22 22:51:39  abbott
// Changed name of fn ndigits to NumDigits.
// Changed return type of exponent and NumDigits.
// Changed some exceptions from nonstandard to the appropriate standard one.
//
// Revision 1.4  2007/05/21 12:57:28  abbott
// New class for passing machine integers as args; includes some simple
// operations on machine integers (cmp, gcd, IsNegative,...).  Operations
// between ZZ and machine integers modified to use the new class.  Inexact
// integer division (of a ZZ) by a negative value now triggers an error;
// new error for reporting inexact integer division by a negative value.
//
// Revision 1.3  2007/03/23 18:38:43  abbott
// Separated the "convert" functions (and CheckedCast) into their own files.
// Many consequential changes.  Also corrected conversion to doubles.
//
// Revision 1.2  2007/03/16 17:43:05  abbott
// Added new convert function (from ZZ to double).
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.10  2007/03/08 18:42:05  cocoa
// Cleaned up whitespace.
//
// Revision 1.9  2007/01/13 14:14:34  cocoa
// Overhaul of RingHom code: it nows uses SmartPtrIRC, and printing is more logical.
// Have not yet updated the documentation.
//
// Revision 1.8  2006/12/06 17:17:44  cocoa
// -- removed #include "config.H"
//
// Revision 1.7  2006/11/29 11:59:35  cocoa
// -- fixed: convert(double& z, const ZZ& num, const ZZ& den) now returns
//    bool (was void) and does not throw errors
//
// Revision 1.6  2006/11/27 13:38:18  cocoa
// -- added #include <cstddef> //using std::size_t
//
// Revision 1.5  2006/11/27 13:06:23  cocoa
// Anna and Michael made me check without writing a proper message.
//
// Revision 1.4  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.3  2006/09/29 11:46:54  cocoa
// Corrected bug in convert(ZZ, ZZ, double) -- now correct and simpler.
// Previously went into infinite loop on negative doubles.
//
// Revision 1.2  2006/06/21 17:08:04  cocoa
// Removed a pointless const.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.4  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.3  2006/03/27 12:21:26  cocoa
// Minor silly changes to reduce number of complaints from some compiler or other.
//
// Revision 1.2  2006/03/14 15:01:49  cocoa
// Improved the implementation of ring member fns for computing powers.
// Should keep Intel C++ compiler quieter too.
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.2  2005/08/08 16:36:32  cocoa
// Just checking in before going on holiday.
// Don't really recall what changes have been made.
// Added IsIndet function for RingElem, PPMonoidElem,
// and a member function of OrdvArith.
// Improved the way failed assertions are handled.
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.2  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.7  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.6  2004/11/11 13:20:47  cocoa
// -- added \include *.txt for doxygen
//
// Revision 1.5  2004/11/04 18:47:42  cocoa
// (1) Ring member functions which previously expected mpz_t args
//     now expect ZZ args.  Numerous minor consequential changes.
// (2) Renamed function which gives access to the mpz_t value inside
//     a ZZ object: previously was raw(...), now is mpzref(...).
//     Plenty of calls had to be altered.
//
// Revision 1.4  2004/05/27 16:14:02  cocoa
// Minor revision for new coding conventions.
//
// Revision 1.3  2004/02/03 16:16:20  cocoa
// Removed pointless IamGCDDomain functions from several concrete rings.
// Added IamOrderedDomain functions where appropriate.
// Tidied ctors for the small finite fields.
//
// Revision 1.2  2003/10/17 10:51:06  cocoa
// Major cleaning, and new naming convention.
//
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//
// Revision 1.2  2003/06/23 16:12:44  abbott
// Minor cleaning prior to public release.
//
// Revision 1.1  2003/05/14 17:12:35  abbott
// Initial revision
//
//

#endif
