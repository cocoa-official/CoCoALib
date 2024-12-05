#ifndef CoCoA_FloatApprox_H
#define CoCoA_FloatApprox_H

//   Copyright (c)  2014  John Abbott,  Anna M. Bigatti

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

#include "CoCoA/assert.H"
#include "CoCoA/BigInt.H"
#include "CoCoA/MachineInt.H"

#include <iosfwd>
//using std::ostream;


namespace CoCoA
{

  class BigRat; // fwd decl -- defined in BigRat.H
  class RingElem; // fwd decl -- defined in ring.H

  class MantExp2
  {
  public:
    MantExp2(): mySign(0), myExponent(0), myMantissa(/*0*/), myNumDigits(0) {}
    MantExp2(int s, long e, const BigInt& m, long NumDigits);
    // default copy ctor, assignment and dtor are all OK
  public: // data members
    int mySign;
    long myExponent;
    BigInt myMantissa;
    long myNumDigits;

    static const int ourDefaultMantBits; // actual value is in FloatApprox.C
  };

  std::ostream& operator<<(std::ostream& out, const MantExp2& ME);


  MantExp2 MantissaAndExponent2(const MachineInt& n, const MachineInt& MantBits = MantExp2::ourDefaultMantBits);
  MantExp2 MantissaAndExponent2(const BigInt& N, const MachineInt& MantBits = MantExp2::ourDefaultMantBits);
  MantExp2 MantissaAndExponent2(const BigRat& q, const MachineInt& MantBits = MantExp2::ourDefaultMantBits);
  MantExp2 MantissaAndExponent2(const RingElem& x); // only for x in RingTwinFloat, def'd in RingTwinFloat.C

  BigRat FloatApprox(const MachineInt& n, const MachineInt& MantBits = MantExp2::ourDefaultMantBits);
  BigRat FloatApprox(const BigInt& n, const MachineInt& MantBits = MantExp2::ourDefaultMantBits);
  BigRat FloatApprox(const BigRat& n, const MachineInt& MantBits = MantExp2::ourDefaultMantBits);

  //------------------------------------------------------------------
  // Decimal "floating point" representation

  class MantExp10
  {
  public:
    MantExp10(): mySign(0), myExponent(0), myMantissa(/*0*/), myNumDigits(0) {}
    MantExp10(int s, long e, const BigInt& m, long NumDigits);
    // default copy ctor, assignment and dtor are all OK
  public: // data members
    int mySign;
    long myExponent;
    BigInt myMantissa;
    long myNumDigits;

    static const int ourDefaultSigFig; // actual value is in decimal.C
  };

  std::ostream& operator<<(std::ostream& out, const MantExp10& ME);

  MantExp10 MantissaAndExponent10(const BigInt& N, const MachineInt& SigFig = MantExp10::ourDefaultSigFig);
  MantExp10 MantissaAndExponent10(const BigRat& N, const MachineInt& SigFig = MantExp10::ourDefaultSigFig);


} // end of namespace CoCoA

#endif
