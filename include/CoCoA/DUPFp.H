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

#endif
