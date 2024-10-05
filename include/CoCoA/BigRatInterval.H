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



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/BigRatInterval.H,v 1.5 2022/02/18 14:11:44 abbott Exp $
// $Log: BigRatInterval.H,v $
// Revision 1.5  2022/02/18 14:11:44  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.4  2018/05/22 14:16:39  abbott
// Summary: Split BigRat into BigRat (class defn + ctors) and BigRatOps
//
// Revision 1.3  2018/04/20 13:10:17  abbott
// Summary: Corrected fn name
//
// Revision 1.2  2018/04/20 12:55:30  abbott
// Summary: Added arith ops betw BigRatInterval and BigRat; added fn merge
//
// Revision 1.1  2018/04/18 14:15:21  abbott
// Summary: New files for BigRatInterval
//
//

#endif
