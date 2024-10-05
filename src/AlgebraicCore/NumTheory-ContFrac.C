//   Copyright (c)  1999,2009-2011  John Abbott and Anna M. Bigatti

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


#include "CoCoA/NumTheory-ContFrac.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRatOps.H"
#include "CoCoA/assert.H"
#include "CoCoA/error.H"

#include <iostream>
using std::ostream;


namespace CoCoA
{

  const PlusOrMinusEpsilon PlusEpsilon = PlusOrMinusEpsilon::PlusEpsilon;
  const PlusOrMinusEpsilon MinusEpsilon = PlusOrMinusEpsilon::MinusEpsilon;


  // Some developer documentation:
  // ContFracIter has 3 data fields: myFrac, myQuot, myEpsilon
  // In "standard use" myEpsilon is zero, and stays zero at all times.
  // Let Q be the original input rational, and let the EPS be initial
  // input infinitesimal (its sign is +1, 0 or -1)
  // Let output "quotients" be q_0, q_1, q_2, etc.  These "quotients"
  // are the successive values placed in myQuot (see operator++).
  // q_0 is floor(Q), and is computed (in myQuot) in the ctor.
  // q_0 can be any integer; all the other q_k are positive integers,
  // *EXCEPT* perhaps the very last one which can be 0 (meaning +infinity)
  // if an infinitesimal was present in the original input.
  
  // The iterator has ended when myFirstQuot is false and myQuot=0 and myFrac=0.

  // Invariant:
  // (A) if myFrac = 0 then
  //     Q = [q_0, q_1, ..., q_k]   <--- RHS is a continued fraction
  // (B) otherwise myFrac > 0, and then
  //     Q+EPS = [q_0, q_1, q_k, myFrac+myEpsilon]
  //     where we consider only the signs of the infinitesimals.

  ContFracIter::ContFracIter(const BigRat& Q, PlusOrMinusEpsilon EPS)
  {
    switch (EPS)
    {
      case PlusOrMinusEpsilon::ZeroEpsilon: myEpsilon = 0; break;
      case PlusOrMinusEpsilon::PlusEpsilon: myEpsilon = 1; break;
      case PlusOrMinusEpsilon::MinusEpsilon: myEpsilon = -1; break;
    }
    myQuotIndex = 0;
    if (IsOneDen(Q) && myEpsilon >= 0)
    { myQuot = num(Q); /* myFrac = 0; */ myEpsilon = -myEpsilon; return; }
    myQuot = floor(Q);
    if (IsOneDen(Q)) myQuot -= 1;
    myNum = den(Q);
    myDen = num(Q) - myQuot*myNum;
    myEpsilon = -myEpsilon;
  }

  const BigInt& ContFracIter::operator*() const
  {
    if (IsEnded(*this)) CoCoA_THROW_ERROR(ERR::IterEnded, "ContFracIter::operator*");
    return myQuot;
  }

  ContFracIter& ContFracIter::operator++()
  {
    if (IsEnded(*this)) CoCoA_THROW_ERROR(ERR::IterEnded, "ContFracIter::operator++");
    ++myQuotIndex;
    if (IsZero(myNum)) { /*signal iter ended:*/ myQuot = 0; return *this; }
    if (!IsOne(myDen))
    { // the "usual" case
      quorem(myQuot, myNum, myNum, myDen);
      swap(myNum, myDen);
    }
    else
    { // myDen = 1, so almost ended.
      if (myEpsilon < 0) { myQuot = myNum-1;  myNum = 1; }
      else { myQuot = myNum; myNum = 0; }
    }
    myEpsilon = -myEpsilon;
    return *this;
  }

  ContFracIter ContFracIter::operator++(int)
  {
    ContFracIter prev = *this;
    operator++();
    return prev;
  }


  bool IsEnded(const ContFracIter& CFIter) noexcept
  {
    return (CFIter.myQuotIndex != 0) && IsZero(CFIter.myQuot) && IsZero(CFIter.myNum);
  }


  bool IsFinal(const ContFracIter& CFIter) noexcept
  {
    return IsZero(CFIter.myNum); // && !IsZero(CFIter.myQuot)    BUG BUG BUG ???
  }


  std::ostream& operator<<(std::ostream& out, const ContFracIter& CFIter)
  {
    if (!out) return out;  // short-cut for bad ostreams
    out << "ContFracIter(myNum = " << CFIter.myNum
        << ", myDen = " << CFIter.myDen
        << ", myQuotIndex = " << CFIter.myQuotIndex
        << ", myQuot = " << CFIter.myQuot;
    if (CFIter.myEpsilon != 0)
      out << ", myEpsilon = " << CFIter.myEpsilon;
    out << ")";
    return out;
  }


  //////////////////////////////////////////////////////////////////

  ContFracApproximant::ContFracApproximant():
      myCurr(BigRat::OneOverZero), // WARNING: anomalous value, 1/0
      myPrev(0,1)
  {}


  void ContFracApproximant::myAppendQuot(const MachineInt& q)
  {
    // Simple rather than fast.
    myAppendQuot(BigInt(q));
  }

  void ContFracApproximant::myAppendQuot(const BigInt& q)
  {
    // These 9 lines should avoid (all explicit) temporaries:
    // NB I have to use pointers to mpq_t because GMP's design won't let me use references.
    mpq_t* prev = &mpqref(myPrev);
    mpq_t* curr = &mpqref(myCurr);
    mpq_t* next = &mpqref(myNext);
    mpz_mul(mpq_numref(*next), mpq_numref(*curr), mpzref(q));
    mpz_add(mpq_numref(*next), mpq_numref(*next), mpq_numref(*prev));
    mpz_mul(mpq_denref(*next), mpq_denref(*curr), mpzref(q));
    mpz_add(mpq_denref(*next), mpq_denref(*next), mpq_denref(*prev));
    swap(myCurr, myPrev);
    swap(myNext, myCurr);
  }


  std::ostream& operator<<(std::ostream& out, const ContFracApproximant& CFConv)
  {
    if (!out) return out;  // short-cut for bad ostreams
    out << "ContFracApproximant(myCurr = " << CFConv.myCurr
        << ",  myPrev = " << CFConv.myPrev << ")";
    return out;
  }

  //////////////////////////////////////////////////////////////////


  // NB if (Q == 0) then myCFIter starts off "ended"
  CFApproximantsIter::CFApproximantsIter(const BigRat& Q):
      myCFIter(Q),
      myApproximant()
  {
    if (!IsEnded(myCFIter))
      myApproximant.myAppendQuot(quot(myCFIter));
  }

  CFApproximantsIter::CFApproximantsIter(const ContFracIter& CFIter):
      myCFIter(CFIter),
      myApproximant()
  {
    if (!IsEnded(myCFIter))
      myApproximant.myAppendQuot(quot(myCFIter));
  }


  CFApproximantsIter& CFApproximantsIter::operator++()
  {
    if (IsEnded(*this)) CoCoA_THROW_ERROR(ERR::IterEnded, "CFApproximantsIter::operator++");
    ++myCFIter;
    if (IsEnded(myCFIter)) return *this;
    myApproximant.myAppendQuot(quot(myCFIter));

    return *this;
  }

  CFApproximantsIter CFApproximantsIter::operator++(int)
  {
    CFApproximantsIter prev = *this;
    operator++();
    return prev;
  }


  std::ostream& operator<<(std::ostream& out, const CFApproximantsIter& CFAIter)
  {
    if (!out) return out;  // short-cut for bad ostreams
    out << "CFApproximantsIter(myApproximant = " << CFAIter.myApproximant
        << ",  myCFIter = " << CFAIter.myCFIter << ")";
    return out;
  }


  // Return first cont frac convergent having rel error at most MaxRelErr
  BigRat CFApprox(const BigRat& q, const BigRat& MaxRelErr)
  {
    // Simple rather than superfast.
    if (MaxRelErr < 0 || MaxRelErr > 1) CoCoA_THROW_ERROR(ERR::BadArg, "CFApprox: relative error must be between 0 and 1");
    if (IsZero(q) || IsZero(MaxRelErr)) return q;
    const BigRat MaxAbsErr = abs(q*MaxRelErr);
    if (MaxAbsErr >= 1) return BigRat(floor(q),1);
    CFApproximantsIter CFAIter(q);
    while (abs(q - *CFAIter) > MaxAbsErr)
    {
      ++CFAIter;
    }
    return *CFAIter;
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/NumTheory-ContFrac.C,v 1.7 2022/02/18 14:11:55 abbott Exp $
// $Log: NumTheory-ContFrac.C,v $
// Revision 1.7  2022/02/18 14:11:55  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.6  2021/02/10 19:40:00  abbott
// Summary: Added noexcept (sometimes instead of throw()) -- redmine 1572
//
// Revision 1.5  2021/01/07 15:16:51  abbott
// Summary: Corrected copyright
//
// Revision 1.4  2020/11/03 20:02:40  abbott
// Summary: Replaced BigRat computation by BigInt (redmine 897)
//
// Revision 1.3  2020/06/17 15:49:24  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.2  2020/02/11 16:12:18  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.1  2019/03/18 11:24:19  abbott
// Summary: Split NumTheory into several smaller files
//
//
