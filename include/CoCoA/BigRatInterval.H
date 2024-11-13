#ifndef CoCoA_BigRatInterval_H
#define CoCoA_BigRatInterval_H

//   Copyright (c)  2018  John Abbott,  Anna M. Bigatti

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


#include "CoCoA/BigRatOps.H"

#include <iosfwd>
//using std::ostream;

namespace CoCoA
{

  class BigRatInterval
  {
  public:
    BigRatInterval(const BigRat& lwb, const BigRat& upb);
    friend const BigRat& min(const BigRatInterval& I) { return I.myMin; }
    friend const BigRat& max(const BigRatInterval& I) { return I.myMax; }
  private:
    BigRat myMin;
    BigRat myMax;
  };

  BigRatInterval operator+(const BigRat& A, const BigRatInterval& B);
  BigRatInterval operator+(const BigRatInterval& A, const BigRat& B);
  BigRatInterval operator+(const BigRatInterval& A, const BigRatInterval& B);

  BigRatInterval operator-(const BigRat& A, const BigRatInterval& B);
  BigRatInterval operator-(const BigRatInterval& A, const BigRat& B);
  BigRatInterval operator-(const BigRatInterval& A, const BigRatInterval& B);

  BigRatInterval operator*(const BigRat& A, const BigRatInterval& B);
  BigRatInterval operator*(const BigRatInterval& A, const BigRat& B);
  BigRatInterval operator*(const BigRatInterval& A, const BigRatInterval& B);

  BigRatInterval operator/(const BigRatInterval& A, const BigRatInterval& B);
  bool IsZeroInside(const BigRatInterval& A);
  BigRat width(const BigRatInterval& A);
  BigRatInterval square(const BigRatInterval& A);
  BigRatInterval merge(const BigRatInterval& A, const BigRatInterval& B);
  BigRatInterval soften(const BigRatInterval& A); // slightly bigger interval with "nice" end points
  std::ostream& operator<<(std::ostream& out, const BigRatInterval& I);


  // below are inline impls

  inline BigRat width(const BigRatInterval& A)
  {
    return max(A) - min(A);
  }


  inline BigRatInterval operator+(const BigRat& A, const BigRatInterval& B)
  {
    return BigRatInterval(A+min(B), A+max(B));
  }

  inline BigRatInterval operator+(const BigRatInterval& A, const BigRat& B)
  {
    return BigRatInterval(min(A)+B, max(A)+B);
  }

  inline BigRatInterval operator+(const BigRatInterval& A, const BigRatInterval& B)
  {
    return BigRatInterval(min(A)+min(B), max(A)+max(B));
  }


  inline BigRatInterval operator-(const BigRat& A, const BigRatInterval& B)
  {
    return BigRatInterval(A-max(B), A-min(B));
  }
    
  inline BigRatInterval operator-(const BigRatInterval& A, const BigRat& B)
  {
    return BigRatInterval(min(A)-B, max(A)-B);
  }

  inline BigRatInterval operator-(const BigRatInterval& A, const BigRatInterval& B)
  {
    return BigRatInterval(min(A)-max(B), max(A)-min(B));
  }

} // end of namespace CoCoA

#endif
