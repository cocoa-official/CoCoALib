#ifndef CoCoA_DenseUPolyClean_H
#define CoCoA_DenseUPolyClean_H

//   Copyright (c)  2007  Anna Bigatti

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


// Header file for the class DenseUPolyClean

#include "CoCoA/assert.H"
#include "CoCoA/ring.H"

#include <iosfwd>
// using std::ostream;
#include <vector>
// using std::vector;


namespace CoCoA
{


/*-- class DenseUPolyClean ----------------------------------------------*/
/**

  \brief Implementation of "clean" univariate polynomials

  \include DenseUPolyClean.txt
*/
/*-----------------------------------------------------------------*/
  class DenseUPolyClean
  {
  private:
    const ring myCoeffRingValue;   ///< the ring to which all the coeffs belong
    long myDegPlus1Value;   ///< always non-negative (zero only if poly is 0?)
    long mySizeValue;       ///< always non-negative
    std::vector<RingElem> myCoeffsValue;  ///< vector of coefficients

  public:
    DenseUPolyClean(const ring& CoeffRing, long MinCapacity);
    DenseUPolyClean(const DenseUPolyClean& copy, long MinCapacity);
    DenseUPolyClean(const DenseUPolyClean& copy) = delete; ///< no default copy ctor
    ~DenseUPolyClean();
    DenseUPolyClean& operator=(const DenseUPolyClean& rhs);
    DenseUPolyClean& operator=(const MachineInt& rhs);
    DenseUPolyClean& operator=(const BigInt& rhs);
    DenseUPolyClean& operator=(const BigRat& rhs);

    long myDegPlus1() const  { return myDegPlus1Value; }
    long mySize() const  { return mySizeValue; }

    static void ourSwap(DenseUPolyClean&, DenseUPolyClean&);
    void myAssignZero();
    void myAssignZeroCoeff(long d); ///< f_d = 0, no check on size
    void myAssignNonZeroCoeff(ConstRefRingElem c, long d); ///< f_d = c, no check on size
    void myResize(long NewSize);
    void myResetDeg(); ///< reset the correct value of deg (assumed to be less than the current value)
    ConstRefRingElem myCoeff(long d) const;

    friend ring CoeffRing(const DenseUPolyClean& f);
    friend long NumTerms(const DenseUPolyClean& f);
    friend ConstRefRingElem LC(const DenseUPolyClean& f);
    void myAddMulLM(ConstRefRingElem c, long d, const DenseUPolyClean& g);
    //    void myAddClear(DenseUPolyClean& g);
    void myMulByCoeff(ConstRefRingElem rawc);  // EXCEPTION SAFE?
    void myDivByCoeff(ConstRefRingElem rawc);  // EXCEPTION SAFE?
    void myMulByXExp(long n);
    void myMulByXExp(const BigInt& N);
    void myMulBy1MinusXExp(long n);
    void myMulBy1MinusXExp(const BigInt& N);
    void myNegate();

    friend void add(DenseUPolyClean& lhs, const DenseUPolyClean& g, const DenseUPolyClean& h);
    friend void sub(DenseUPolyClean& lhs, const DenseUPolyClean& g, const DenseUPolyClean& h);
    friend bool div(DenseUPolyClean& lhs, const DenseUPolyClean& g, const DenseUPolyClean& h);
    friend bool IsCompatible(const DenseUPolyClean& f, const DenseUPolyClean& g);
  };

  void deriv(DenseUPolyClean& lhs, const DenseUPolyClean&);


  // Friend functions - must also be declared outside the class (how very tedious).
  bool IsCompatible(const DenseUPolyClean& f, const DenseUPolyClean& g);
  long NumTerms(const DenseUPolyClean& f);
  ConstRefRingElem LC(const DenseUPolyClean& f);
  void add(DenseUPolyClean& lhs, const DenseUPolyClean& g, const DenseUPolyClean& h);
  void sub(DenseUPolyClean& lhs, const DenseUPolyClean& g, const DenseUPolyClean& h);
  bool div(DenseUPolyClean& lhs, const DenseUPolyClean& g, const DenseUPolyClean& h);
  bool IsZero(const DenseUPolyClean& f);
  bool IsEqual(const DenseUPolyClean& f, const DenseUPolyClean& g);
  void deriv(DenseUPolyClean& lhs, const DenseUPolyClean& f);

  // Inline functions.

  inline ConstRefRingElem DenseUPolyClean::myCoeff(long d) const
  {
    return myCoeffsValue[d];
  }


  inline void swap(DenseUPolyClean& f, DenseUPolyClean& g)
  {
    DenseUPolyClean::ourSwap(f, g);
  }

} // end of namespace CoCoA

#endif
