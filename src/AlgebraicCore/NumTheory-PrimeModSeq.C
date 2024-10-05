//   Copyright (c)  2022,2023  John Abbott and Anna M. Bigatti

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


#include "CoCoA/NumTheory-PrimeModSeq.H"

#include <limits>
#include <ostream>

namespace CoCoA
{

  /*static member*/ unsigned long PrimeSeq1ModN::ourComputeStepSize(long n)
  {
    const char* const FnName = "PrimeSeq1ModN ctor: arg must be >= 2";
    if (n < 2)
      CoCoA_THROW_ERROR(ERR::BadModulus, FnName);
    if (IsEven(n))
    {
      if (n > std::numeric_limits<long>::max()-1)  CoCoA_THROW_ERROR(ERR::ArgTooBig, FnName);
      return n;
    }
    // n is odd
    if (n > (std::numeric_limits<long>::max()-1)/2)  CoCoA_THROW_ERROR(ERR::ArgTooBig, FnName);
    return 2*n;
  }


  PrimeSeq1ModN::PrimeSeq1ModN(long n) : myModulus(ourComputeStepSize(n)), // will throw if n too small or too large
                                         myCurrPrime(0, ArgIsPrime), // this value is overwritten by call to myAdvance
                                         UPB(std::numeric_limits<long>::max()-myModulus)
  {
    // arg check already done by ourComputeStepSize
    myAdvance(1);
  }

  // advance to next prime greater than n, and 1 mod M
  void PrimeSeq1ModN::myAdvance(unsigned long n)  // arg type???  MachineInt?  Unsigned long??
  {
    do
    {
      if (n > UPB)
      {
        n = 0;
        break;
      }
      n += myModulus;
    } while (!IsPrime(n));
    myCurrPrime = SmallPrime(n, ArgIsPrime);
  }

  std::ostream &operator<<(std::ostream &out, const PrimeSeq1ModN &PSeq)
  {
    if (!out)
      return out; // short-cut for bad ostreams
    out << "PrimeSeq1ModN(curr=" << *PSeq << ")";
    return out;
  }

}


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/NumTheory-PrimeModSeq.C,v 1.5 2023/10/05 17:52:52 abbott Exp $
// $Log: NumTheory-PrimeModSeq.C,v $
// Revision 1.5  2023/10/05 17:52:52  abbott
// Summary: Removed arbitrary upper limit
//
// Revision 1.4  2023/02/22 20:20:47  abbott
// Summary: Moved impl from header file to source file
//
// Revision 1.3  2023/01/31 12:28:45  abbott
// Summary: Changed arg type of myAdvance to unsigned long
//
// Revision 1.2  2022/11/24 21:16:11  abbott
// Summary: Added copyright notice
//
//
