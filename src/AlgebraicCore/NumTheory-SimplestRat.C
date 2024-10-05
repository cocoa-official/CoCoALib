//   Copyright (c)  1999,2009-2011  John Abbott  &  Anna M. Bigatti

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


#include "CoCoA/NumTheory-SimplestRat.H"
#include "CoCoA/NumTheory-ContFrac.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRatOps.H"
#include "CoCoA/assert.H"


namespace CoCoA
{

  BigRat SimplestBigRatBetween(const BigRat& lo, const BigRat& hi)
  {
    if (IsZero(lo) || IsZero(hi)) return BigRat(0);
    if (sign(lo) != sign(hi)) return BigRat(0);
    if (lo == hi) return lo;
    if (sign(lo) < 0) return -SimplestBigRatBetween(-hi, -lo);
    if (lo > hi) return SimplestBigRatBetween(hi,lo);
    // Now we have 0 < lo < hi.
    if (IsOneDen(lo) || IsOneDen(hi)) return BigRat(ceil(lo), 1);
    // Now lo & hi have the same integer part.
    ContFracIter L(lo);
    ContFracIter R(hi);
    ContFracApproximant ans;
    while (!IsEnded(L) && !IsEnded(R) && quot(L) == quot(R))
    {
      ans.myAppendQuot(quot(L));
      ++L;
      ++R;
    }
    if (IsEnded(L)) return lo;
    if (IsEnded(R)) return hi;
    if (quot(L) < quot(R))
    {
      if (IsFinal(L))
        ans.myAppendQuot(quot(L));
      else
        ans.myAppendQuot(quot(L)+1);
    }
    else
    {
      if (IsFinal(R))
        ans.myAppendQuot(quot(R));
      else
        ans.myAppendQuot(quot(R)+1);
    }
    return ans.myRational();
  }


  // This is really a bit messy/involved; isn't there a neater way?
  BigRat SimplestBinaryRatBetween(const BigRat& lo, const BigRat& hi)
  {
    if (IsZero(lo) || IsZero(hi)) return BigRat(0);
    if (sign(lo) != sign(hi)) return BigRat(0);
    if (lo == hi) return lo;
    if (sign(lo) < 0) return -SimplestBinaryRatBetween(-hi, -lo);
    if (lo > hi) return SimplestBinaryRatBetween(hi,lo);

    // Now we have 0 < lo < hi.
    // First check whether there is an integer between lo and hi.
    const BigInt L = ceil(lo);
    const BigInt H = floor(hi);
    if (H >= L)
    {
      // There is an integer in the interval.
      if (H == L) return BigRat(L,1);
      const long exp = FloorLog2(H-L);
      const BigInt pwr2 = power(2,exp);
      const BigInt q = H/pwr2; // integer division!
      // BigInt q,r;
      // quorem(q,r, H,pwr2);
      // if (IsZero(r)) // both L and H are divisible by pwr2
      // {
      //   if (IsEven(L/pwr2)) return BigRat(L,1); else return BigRat(H,1);
      // }
      if (IsEven(q) || (q-1)*pwr2 < L)
        return BigRat(q*pwr2,1);
      return BigRat((q-1)*pwr2,1);
    }

    // lo & hi have the same integer part.
    const BigRat recip = 1/(hi-lo);
    long exp = FloorLog2(recip);
    if (!IsOneDen(recip) || !IsPowerOf2(num(recip))) ++exp;
    const BigInt lwb = ceil(power(2,exp)*lo);
    const BigInt upb = floor(power(2,exp)*hi);
    if (lwb == upb) return BigRat(lwb,power(2,exp));
    CoCoA_ASSERT(lwb+1 == upb);
    if (IsEven(lwb)) return BigRat(lwb,power(2,exp));
    return BigRat(upb,power(2,exp));
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/NumTheory-SimplestRat.C,v 1.3 2022/02/18 14:11:55 abbott Exp $
// $Log: NumTheory-SimplestRat.C,v $
// Revision 1.3  2022/02/18 14:11:55  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.2  2021/01/15 16:59:33  abbott
// Summary: Redmine 1563: new ctor for BigRat direct from integer
//
// Revision 1.1  2020/01/26 14:14:31  abbott
// Summary: Finished splitting NumTheory into smaller pieces (redming 1161)
//
//
//
