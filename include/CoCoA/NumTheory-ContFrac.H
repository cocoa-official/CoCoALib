#ifndef CoCoA_NumTheory_ContFrac_H
#define CoCoA_NumTheory_ContFrac_H

//   Copyright (c)  2009,2012  John Abbott,  Anna M. Bigatti

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

#include "CoCoA/BigInt.H"
#include "CoCoA/BigRat.H"

#include <iosfwd>
// using std::ostream;


namespace CoCoA
{

  enum class PlusOrMinusEpsilon { MinusEpsilon, ZeroEpsilon, PlusEpsilon };
  extern const PlusOrMinusEpsilon PlusEpsilon;
  extern const PlusOrMinusEpsilon MinusEpsilon;
  
  class ContFracIter
  {
  public:
    explicit ContFracIter(const BigRat& Q, PlusOrMinusEpsilon EPS = PlusOrMinusEpsilon::ZeroEpsilon);
    const BigInt& operator*() const; // throws if IsEnded
    ContFracIter& operator++();
    ContFracIter operator++(int);
    friend bool IsEnded(const ContFracIter& CFIter) noexcept;
    friend bool IsFinal(const ContFracIter& CFIter) noexcept;
  private: // data members
//    BigRat myFrac;
    BigInt myNum, myDen;
    BigInt myQuot;
    int myEpsilon; // +1, 0 or -1
    long myQuotIndex;
//    bool myFirstQuot;
    friend std::ostream& operator<<(std::ostream& out, const ContFracIter& CFIter);
  };

  bool IsEnded(const ContFracIter& CFIter) noexcept;
  bool IsFinal(const ContFracIter& CFIter) noexcept;

  inline const BigInt& quot(const ContFracIter& CFIter)  // throws if IsEnded(CFIter)
  { return *CFIter; }

  //////////////////////////////////////////////////////////////////

  class ContFracApproximant
  {
  public:
    ContFracApproximant();
    void myAppendQuot(const MachineInt& q);
    void myAppendQuot(const BigInt& q);
    const BigRat& myRational() const noexcept { return myCurr; }
  private: // data members
    BigRat myCurr; // WARNING: initial value is anomalous, namely 1/0
    BigRat myPrev;
    BigRat myNext; // Just workspace; used only in operator++
    friend std::ostream& operator<<(std::ostream& out, const ContFracApproximant& CFConv);
  };

  std::ostream& operator<<(std::ostream& out, const ContFracApproximant& CFConv);

  //////////////////////////////////////////////////////////////////

  class CFApproximantsIter
  {
  public:
    explicit CFApproximantsIter(const BigRat& Q);
    explicit CFApproximantsIter(const ContFracIter& CFIter);
    const BigRat& operator*() const noexcept { return myApproximant.myRational(); }
    CFApproximantsIter& operator++();
    CFApproximantsIter operator++(int);
    friend bool IsEnded(const CFApproximantsIter& CFAIter) noexcept { return IsEnded(CFAIter.myCFIter); }
  private: // data members
    ContFracIter myCFIter;
    ContFracApproximant myApproximant;
    friend std::ostream& operator<<(std::ostream& out, const CFApproximantsIter& CFAIter);
  };

  std::ostream& operator<<(std::ostream& out, const ContFracIter& CFIter);
  std::ostream& operator<<(std::ostream& out, const CFApproximantsIter& CFAIter);

  BigRat CFApprox(const BigRat& q, const BigRat& MaxRelErr);


} // end of namespace CoCoA

#endif
