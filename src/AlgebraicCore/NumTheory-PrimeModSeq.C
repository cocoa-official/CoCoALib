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
    if (n < 2)
      CoCoA_THROW_ERROR1(ERR::BadModulus);
    if (IsEven(n))
    {
      if (n > std::numeric_limits<long>::max()-1)
        CoCoA_THROW_ERROR1(ERR::ArgTooBig);
      return n;
    }
    // n is odd
    if (n > (std::numeric_limits<long>::max()-1)/2)
      CoCoA_THROW_ERROR1(ERR::ArgTooBig);
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
