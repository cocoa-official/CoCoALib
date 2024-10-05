#ifndef CoCoA_PrimeModSeq
#define CoCoA_PrimeModSeq

//   Copyright (c)  2022  John Abbott,  Anna M. Bigatti

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


#include "CoCoA/BigIntOps.H"
#include "CoCoA/NumTheory-prime.H"

// Not needed?  #include <iosfwd>  // already incl in other header files

namespace CoCoA
{

  class PrimeSeq1ModN
  {
  public:
    explicit PrimeSeq1ModN(long n);
    SmallPrime operator*() const { return myCurrPrime; }
    PrimeSeq1ModN &operator++()
      {
        myAdvance(myCurrPrime);
        return *this;
      }
    PrimeSeq1ModN operator++(int)
      {
        PrimeSeq1ModN copy(*this);
        operator++();
        return copy;
      }
    void JumpTo(unsigned long n)
      { /*ARG CHECK*/
        myAdvance(1 + myModulus * (n / myModulus));
      }
    bool IamEnded() const { return (myCurrPrime == 0); }

  private: // data members
    const unsigned long myModulus;
    SmallPrime myCurrPrime;

  private: // impl details
    void myAdvance(unsigned long n);
    unsigned long UPB;
    static unsigned long ourComputeStepSize(long n);
  };


  inline SmallPrime CurrPrime(const PrimeSeq1ModN &PSeq)
  {
    return *PSeq;
  }


  inline SmallPrime NextPrime(PrimeSeq1ModN &PSeq)
  {
    return *++PSeq;
  }


  inline bool IsEnded(const PrimeSeq1ModN &PSeq)
  {
    return PSeq.IamEnded();
  }

  std::ostream &operator<<(std::ostream &out, const PrimeSeq1ModN &PSeq);


} // end of namespace CoCoA

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/NumTheory-PrimeModSeq.H,v 1.6 2024/01/21 19:54:27 abbott Exp $
// $Log: NumTheory-PrimeModSeq.H,v $
// Revision 1.6  2024/01/21 19:54:27  abbott
// Summary: JumpTo now takes UNSIGNED arg
//
// Revision 1.5  2023/02/22 20:20:46  abbott
// Summary: Moved impl from header file to source file
//
// Revision 1.4  2023/01/31 12:26:08  abbott
// Summary: Changed internal repr to unsigned long
//
// Revision 1.3  2022/11/24 21:17:05  abbott
// Summary: Added commented out include (why not needed?)
//
// Revision 1.2  2022/11/24 21:12:10  abbott
// Summary: Added copyright notice (& some spacing)
//


#endif
