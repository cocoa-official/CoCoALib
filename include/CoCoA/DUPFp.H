#ifndef CoCoA_DUPFp_H
#define CoCoA_DUPFp_H

//   Copyright (c)  2012  John Abbott,  Anna M. Bigatti

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

#include "CoCoA/SmallFpImpl.H"
#include "CoCoA/factorization.H"

#include <iosfwd>
//using std::ostream;
#include <vector>
//using std::vector;

namespace CoCoA
{

  struct DUPFp
  {
  public:
    DUPFp(long maxdeg, const SmallFpImpl& arith);//: myArith(arith) { myCoeffs.reserve(maxdeg); }
    DUPFp(const DUPFp&) = default; // to avoid compiler warning (g++ v.9)
    DUPFp& operator=(const DUPFp& rhs);
  public: // data members
    const SmallFpImpl& myArith;
    std::vector<SmallFpImpl::value> myCoeffs;  // convention last elem is non-zero (or vector is empty)
  };

  void AssignZero(DUPFp& f);
  void AssignOne(DUPFp& f);
  bool IsZero(const DUPFp& f);
  long deg(const DUPFp& f); // deg(0) = -1
  void FixDeg(DUPFp& f);
  SmallFpImpl::value LC(const DUPFp& f);
  void MakeMonic(DUPFp& f);
  DUPFp monic(const DUPFp& f);
  // Mult/div by a coeff:
  DUPFp operator*(const DUPFp& f, SmallFpImpl::value c);
  DUPFp& operator*=(DUPFp& f, SmallFpImpl::value c);
  DUPFp operator/(const DUPFp& f, SmallFpImpl::value c);
  DUPFp& operator/=(DUPFp& f, SmallFpImpl::value c);
  // Arith between polys:
  void add(DUPFp& lhs, const DUPFp& f, const DUPFp& g);
  DUPFp operator+(const DUPFp& f, const DUPFp& g);
  void sub(DUPFp& lhs, const DUPFp& f, const DUPFp& g);
  DUPFp operator-(const DUPFp& f, const DUPFp& g);
  void mul(DUPFp& lhs, const DUPFp& f, const DUPFp& g);
  DUPFp operator*(const DUPFp& f, const DUPFp& g);
  DUPFp operator/(const DUPFp& num, const DUPFp& den);
  DUPFp operator%(const DUPFp& num, const DUPFp& den);
  void ShiftAdd(DUPFp& f, const DUPFp& g, SmallFpImpl::value c, long exp);
  void QuoRem(DUPFp& q, DUPFp& r, const DUPFp& num, const DUPFp& den);

  DUPFp power(const DUPFp& base, long exp);
  DUPFp deriv(const DUPFp& f);
  DUPFp gcd(const DUPFp& f, const DUPFp& g);
  DUPFp exgcd(DUPFp& cf, DUPFp& cg, const DUPFp& f, const DUPFp& g);

  SmallFpImpl::value eval(const DUPFp& f, SmallFpImpl::value x);

  bool operator==(const DUPFp& f, const DUPFp& g);
  bool operator!=(const DUPFp& f, const DUPFp& g);
  std::ostream& operator<<(std::ostream& out, const DUPFp& f);

  DUPFp ConvertToDUPFp(const SmallFpImpl& ModP, ConstRefRingElem f);
  RingElem ConvertFromDUPFp(ConstRefRingElem x, const DUPFp& f);

  DUPFp MulMod(const DUPFp& f, const DUPFp& g, const DUPFp& m);
  DUPFp PowerMod(const DUPFp& f, long exp, const DUPFp& m);

  bool IsSqfr(const DUPFp& f);
  DUPFp PthRoot(const DUPFp& f);

  //-------------------------------------------------------
  // The following are needed to allow factorization objs -- see factorization.H
  template <>
  inline void factorization<DUPFp>::ourCheckNotField(const DUPFp& /*fac*/)
  { /* nothing to do */ }


  template <>
  inline void factorization<DUPFp>::ourCheckNotZeroDiv(const DUPFp& fac)
  {
    if (IsZero(fac))
      CoCoA_THROW_ERROR2(ERR::ReqNonZero, "factorization<DUPFp>");
  }

  template <>
  inline void factorization<DUPFp>::ourCheckNotUnit(const DUPFp& fac)
  {
    CoCoA_ASSERT(!IsZero(fac));
    if (deg(fac) == 0)
      CoCoA_THROW_ERROR2(ERR::BadArg,  "factorization<DUPFp> factor must not be invertible");
  }

  factorization<DUPFp> SqfrDecomp(DUPFp f);
  factorization<DUPFp> DistinctDegFactor(DUPFp f);

  factorization<DUPFp> factor(const DUPFp& f);

  inline DUPFp& DUPFp::operator=(const DUPFp& rhs)
  {
    CoCoA_ASSERT(myArith == rhs.myArith);
    myCoeffs = rhs.myCoeffs;
    return *this;
  }


  inline void swap(DUPFp& a, DUPFp& b) { std::swap(a.myCoeffs, b.myCoeffs); }

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/DUPFp.H,v 1.12 2024/07/31 07:52:46 abbott Exp $
// $Log: DUPFp.H,v $
// Revision 1.12  2024/07/31 07:52:46  abbott
// Summary: Switched to new CoCoA_THROW_ERROR macros (with suffix 1 or 2)
//
// Revision 1.11  2022/03/07 14:28:19  abbott
// Summary: Added explicit default copy ctor to avoid g++ v.9 warning
//
// Revision 1.10  2022/02/18 14:11:45  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.9  2021/01/07 14:11:05  abbott
// Summary: Corrected copyright
//
// Revision 1.8  2020/10/27 09:49:29  abbott
// Summary: Added check that factorization is not over a field (redmine 1499)
//
// Revision 1.7  2020/06/17 15:49:20  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.6  2015/11/30 21:58:13  abbott
// Summary: Minor cleaning; added 2 comments
//
// Revision 1.5  2015/11/04 12:16:27  abbott
// Summary: Consequential changes (after revising SmallFpImpl)
//
// Revision 1.4  2014/03/24 12:12:25  abbott
// Summary: Added template specializations for new factorization impl/interface
// Author: JAA
//
// Revision 1.3  2013/03/25 17:24:28  abbott
// Added several new fns (mostly just "copied" from old C4 code).
//
// Revision 1.2  2013/02/21 12:56:42  abbott
// Added new fns monic, deriv, gcd, and factor.
//
// Revision 1.1  2012/11/23 17:32:16  abbott
// Added new class DUPFp for dense univariate polys over a small finite field.
//
//

#endif
