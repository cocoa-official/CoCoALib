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



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/BigRat.H,v 1.22 2022/10/10 09:53:34 abbott Exp $
// $Log: BigRat.H,v $
// Revision 1.22  2022/10/10 09:53:34  abbott
// Summary: Added move-assign op
//
// Revision 1.21  2022/02/18 14:11:44  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.20  2021/09/01 08:14:30  abbott
// Summary: Added comment
//
// Revision 1.19  2021/08/04 19:04:33  abbott
// Summary: Removed const (redmine 1606)
//
// Revision 1.18  2021/03/03 22:09:32  abbott
// Summary: New enum class (redmine 894)
//
// Revision 1.17  2021/01/15 16:59:33  abbott
// Summary: Redmine 1563: new ctor for BigRat direct from integer
//
// Revision 1.16  2021/01/07 14:10:03  abbott
// Summary: Corrected copyright
//
// Revision 1.15  2020/12/04 13:51:56  abbott
// Summary: Moved basic query fns to BigRat.H; updated BigRat.C (redmine 1529)
//
// Revision 1.14  2020/12/04 11:55:15  abbott
// Summary: Revised (redmine 1529)
//
// Revision 1.13  2020/10/05 19:23:05  abbott
// Summary: Change data member name to myRepr; added move ctor
//
// Revision 1.12  2018/05/25 09:26:46  abbott
// Summary: Minor cleaning
//
// Revision 1.11  2018/05/22 14:16:39  abbott
// Summary: Split BigRat into BigRat (class defn + ctors) and BigRatOps
//
// Revision 1.10  2018/04/20 18:51:25  abbott
// Summary: Changed ctors for BigInt/BigRat from string or from MPZ/MPQ
//
// Revision 1.9  2016/03/25 19:58:40  abbott
// Summary: Added BigRat ctors from MantExp2 and MantExp10 structures
//
// Revision 1.8  2015/11/23 18:19:56  abbott
// Summary: Renamed ILogBase -> FloorLogBase; added FloorLog2, FloorLog10
//
// Revision 1.7  2014/07/15 16:49:34  abbott
// Summary: Corrected include guard capitalization
// Author: JAA
//
// Revision 1.6  2014/06/14 19:25:08  abbott
// Summary: Added new fn CmpAbs (for BigRat)
// Author: JAA
//
// Revision 1.5  2014/05/02 13:55:29  abbott
// Summary: Replaced include by fwd decl
// Author: JAA
//
// Revision 1.4  2014/01/28 09:57:21  abbott
// Added option 2nd arg for conversion string -> BigRat to say whether the result
// should be "canonicalized".
//
// Revision 1.3  2013/05/20 15:50:20  abbott
// Added new ctor for BigRat from mpq_t.
//
// Revision 1.2  2011/11/09 13:38:31  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.1  2011/08/24 10:16:38  bigatti
// -- was file/class QQ
//
// ------- log for old file/class name: QQ.H --------------
// Revision 1.11  2011/08/14 15:52:18  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.10  2011/08/02 11:28:20  abbott
// Fixed typo in comment.
//
// Revision 1.9  2011/06/23 16:01:07  abbott
// Removed single arg ctor QQ(MachineInt), & consequential changes.
//
// Revision 1.8  2011/02/25 12:06:51  abbott
// Added new fn IsOneNum; also some minor code cleaning in QQ.C
//
// Revision 1.7  2010/12/26 13:03:16  abbott
// Added ILogBase function (to ZZ & QQ).
//
// Revision 1.6  2010/05/07 14:57:52  abbott
// Two main changes:
//   power(QQ,ZZ) now allows negative exponent
//   renamed QQ::AlreadyNormalized to QQ::AlreadyReduced
//           (and allowed denoms to be negative; the ctor then makes them positive).
//
// Revision 1.5  2010/03/22 11:49:28  abbott
// Added ctor from a string.
//
// Revision 1.4  2010/03/18 16:34:10  abbott
// Added new pseudo-ctors for QQ with optional flag to indicate that value is already normalized.
// Added OpenMath I/O operators.
//
// Revision 1.3  2009/11/26 16:12:54  bigatti
// -- commented out some includes
//
// Revision 1.2  2009/07/08 12:26:53  abbott
// Added floor and ceil functions for QQs.
// Added example program for QQs.
// Minor correction to convert.C; minor cleaning to ex-ZZ1.C
//
// Revision 1.1  2009/07/02 16:29:42  abbott
// Added new class QQ to represent rational numbers.
// Consequent change to the Makefile.
//
//

#endif
