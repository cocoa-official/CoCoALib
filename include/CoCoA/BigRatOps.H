#ifndef CoCoA_BigRatOps_H
#define CoCoA_BigRatOps_H

//   Copyright (c)  2018  John Abbott,  Anna M Bigatti

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


#include "CoCoA/BigRat.H"
#include "CoCoA/BigIntOps.H"


#include <mutex>
#include <vector>
//using std::vector;  // only for CommonDenom

namespace CoCoA
{

//  class MantExp10;       // fwd decl -- defined in FloatApprox.H
//  class MantExp2;        // fwd decl -- defined in FloatApprox.H

  // STANDARD ARITHMETIC OPERATIONS

  void swap(BigRat& a, BigRat& b);

  BigRat abs(const BigRat& Q);
  BigRat operator-(const BigRat& Q);

  BigRat operator+(const BigRat& Q1, const BigRat& Q2);
  BigRat operator-(const BigRat& Q1, const BigRat& Q2);
  BigRat operator*(const BigRat& Q1, const BigRat& Q2);
  BigRat operator/(const BigRat& Q1, const BigRat& Q2);

  BigRat operator+(const BigRat& Q, const BigInt& N);
  BigRat operator-(const BigRat& Q, const BigInt& N);
  BigRat operator*(const BigRat& Q, const BigInt& N);
  BigRat operator/(const BigRat& Q, const BigInt& N);

  BigRat operator+(const BigInt& N, const BigRat& Q);
  BigRat operator-(const BigInt& N, const BigRat& Q);
  BigRat operator*(const BigInt& N, const BigRat& Q);
  BigRat operator/(const BigInt& N, const BigRat& Q);

  BigRat operator+(const BigRat& Q, const MachineInt& n);
  BigRat operator-(const BigRat& Q, const MachineInt& n);
  BigRat operator*(const BigRat& Q, const MachineInt& n);
  BigRat operator/(const BigRat& Q, const MachineInt& n);

  BigRat operator+(const MachineInt& n, const BigRat& Q);
  BigRat operator-(const MachineInt& n, const BigRat& Q);
  BigRat operator*(const MachineInt& n, const BigRat& Q);
  BigRat operator/(const MachineInt& n, const BigRat& Q);

  BigRat power(const BigRat& base, const BigInt& exponent);
  BigRat power(const BigRat& base, const MachineInt& exponent);


  class SumBigRat
  {
  public:
    SumBigRat() = default;
    SumBigRat(const SumBigRat&) = default;
    SumBigRat(SumBigRat&&) = default;
    SumBigRat& operator=(const SumBigRat&) = default;
    SumBigRat& operator=(SumBigRat&&) = default;
    ~SumBigRat() = default;
    SumBigRat& operator+=(long n);          // NOT THREADSAFE
    SumBigRat& operator+=(const BigInt& N); // NOT THREADSAFE
    SumBigRat& operator+=(const BigRat& Q); // NOT THREADSAFE
///    void myAdd(const BigInt& N); // threadsafe version of op+=
///    void myAdd(const BigRat& Q); // threadsafe version of op+=
    BigRat myTotal() const;      // compute & return the total (not threadsafe)
  private: // data members
///    std::mutex myMutex; // causes too much complication
    SumBigInt myBigIntPart;
    mutable std::vector<BigRat> myGeobucket;
  private:  // impl details
    static size_t ourIndex(const BigRat& Q) noexcept;
    void myCascade(size_t i);
    friend std::ostream& operator<<(std::ostream& out, const SumBigRat& S);
  };

  inline BigRat total(const SumBigRat& S)
  { return S.myTotal(); }

  std::ostream& operator<<(std::ostream& out, const SumBigRat& S);


  // COMPARISON FUNCTIONS

  int cmp(const BigRat& Q1, const BigRat& Q2);
  int cmp(const BigRat& Q, const BigInt& N);
  int cmp(const BigInt& N, const BigRat& Q);
  int cmp(const BigRat& Q, const MachineInt& n);
  int cmp(const MachineInt& n, const BigRat& Q);

  int CmpAbs(const BigRat& Q1, const BigRat& Q2);
  int CmpAbs(const BigRat& Q, const BigInt& N);
  int CmpAbs(const BigInt& N, const BigRat& Q);
  int CmpAbs(const BigRat& Q, const MachineInt& n);
  int CmpAbs(const MachineInt& n, const BigRat& Q);

  bool operator==(const BigRat& Q1, const BigRat& Q2);
  bool operator!=(const BigRat& Q1, const BigRat& Q2);
  bool operator> (const BigRat& Q1, const BigRat& Q2);
  bool operator>=(const BigRat& Q1, const BigRat& Q2);
  bool operator< (const BigRat& Q1, const BigRat& Q2);
  bool operator<=(const BigRat& Q1, const BigRat& Q2);
			
  bool operator==(const BigRat& Q, const BigInt& N);
  bool operator!=(const BigRat& Q, const BigInt& N);
  bool operator> (const BigRat& Q, const BigInt& N);
  bool operator>=(const BigRat& Q, const BigInt& N);
  bool operator< (const BigRat& Q, const BigInt& N);
  bool operator<=(const BigRat& Q, const BigInt& N);
			
  bool operator==(const BigInt& N, const BigRat& Q);
  bool operator!=(const BigInt& N, const BigRat& Q);
  bool operator> (const BigInt& N, const BigRat& Q);
  bool operator>=(const BigInt& N, const BigRat& Q);
  bool operator< (const BigInt& N, const BigRat& Q);
  bool operator<=(const BigInt& N, const BigRat& Q);
			
  bool operator==(const BigRat& Q, const MachineInt& n);
  bool operator!=(const BigRat& Q, const MachineInt& n);
  bool operator> (const BigRat& Q, const MachineInt& n);
  bool operator>=(const BigRat& Q, const MachineInt& n);
  bool operator< (const BigRat& Q, const MachineInt& n);
  bool operator<=(const BigRat& Q, const MachineInt& n);
      		
  bool operator==(const MachineInt& n, const BigRat& Q);
  bool operator!=(const MachineInt& n, const BigRat& Q);
  bool operator> (const MachineInt& n, const BigRat& Q);
  bool operator>=(const MachineInt& n, const BigRat& Q);
  bool operator< (const MachineInt& n, const BigRat& Q);
  bool operator<=(const MachineInt& n, const BigRat& Q);
			

  // MISCELLANEOUS FUNCTIONS

  double mantissa(long& exp, const BigRat& Q) noexcept;
//???  double mantissa(const BigRat& Q);
  long exponent(const BigRat& Q);
  double log(const BigRat& Q);
  double LogAbs(const BigRat& Q);
  long FloorLog2(const BigRat& Q);
  long FloorLog10(const BigRat& Q);
  long FloorLogBase(const BigRat& Q, const MachineInt& base);
  long FloorLogBase(const BigRat& Q, const BigInt& base);
  bool IsPowerOf2(const BigRat& Q);

  BigInt floor(const BigRat& Q);   ///< greatest integer <= Q
  BigInt ceil(const BigRat& Q);    ///< least integer >= Q
  BigInt round(const BigRat& Q);   ///< round rational Q to nearest integer (halves round toward +infinity)

  BigInt CommonDenom(const std::vector<BigRat>& L);


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/BigRatOps.H,v 1.15 2024/01/25 21:06:00 abbott Exp $
// $Log: BigRatOps.H,v $
// Revision 1.15  2024/01/25 21:06:00  abbott
// Summary: Removed thread-safety from SumBigInt/SumBigRat (costly/complicated to do it right)
//
// Revision 1.14  2023/12/10 15:31:18  abbott
// Summary: Cleaned and tidied SumBigInt & SumBigRat
//
// Revision 1.13  2023/12/09 21:37:32  abbott
// Summary: Tidied SumBigInt/SumBigRat; added doc
//
// Revision 1.12  2022/10/25 17:25:13  abbott
// Summary: Added mem fn myAdd, threadsafe version of op+=; needs mutex
//
// Revision 1.11  2022/10/11 18:06:36  abbott
// Summary: Added comment about swap, move ctor & mpve assign
//
// Revision 1.10  2022/10/10 09:58:11  abbott
// Summary: Added new class SumBigRat
//
// Revision 1.9  2022/02/18 14:11:44  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.8  2021/08/04 19:04:55  abbott
// Summary: Removed const (redmine 1606)
//
// Revision 1.7  2021/01/31 10:08:43  abbott
// Summary: Added LogAbs
//
// Revision 1.6  2020/12/04 13:51:57  abbott
// Summary: Moved basic query fns to BigRat.H; updated BigRat.C (redmine 1529)
//
// Revision 1.5  2020/10/23 07:53:40  abbott
// Summary: IsZero, IsOne, IsMinusOne are now inline
//
// Revision 1.4  2020/09/22 18:10:00  abbott
// Summary: Added new fn IsPowerOf2
//
// Revision 1.3  2020/02/28 08:53:48  abbott
// Summary: Added CommonDenom
//
// Revision 1.2  2020/01/26 14:15:57  abbott
// Summary: Changed copyright
//
// Revision 1.1  2018/05/22 14:15:14  abbott
// Summary: Split BigRat into BigRat (class defn) and BigRatOps
//
//

#endif
