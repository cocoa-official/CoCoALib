#ifndef CoCoA_NumTheory_CoprimeFactorBasis_H
#define CoCoA_NumTheory_CoprimeFactorBasis_H

//   Copyright (c)  2018  John Abbott, Anna M. Bigatti

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


#include <iosfwd>
//using std::ostream;
#include <vector>
// using std::vector;

namespace CoCoA
{
  
  class MachineInt; // fwd decl -- defined in MachineInt.H

  class CoprimeFactorBasis_BigInt
  {
  public:
    CoprimeFactorBasis_BigInt() {}
    void myAddInfo(const MachineInt& n);
    void myAddInfo(const BigInt& N);
    void myAddInfo(const std::vector<long>& v);
    void myAddInfo(const std::vector<BigInt>& v);
    friend const std::vector<BigInt>& FactorBase(const CoprimeFactorBasis_BigInt& CFB) /*NOT noexcept to allow buffering*/;
  private: // implementation detail
    struct LCR
    {
      LCR(const BigInt& L, const std::vector<BigInt>& C, const BigInt& R):
          myL(L), myR(R), myC(C) {}
    public: // data members
      BigInt myL;
      BigInt myR;
      std::vector<BigInt> myC;
    };
  private: // data members
    std::vector<BigInt> myCoprimeBasis;
    void myRefineBasis(BigInt N);
    struct LCR myLCR(BigInt A, BigInt B) const;
  };


  std::ostream& operator<<(std::ostream& out, const CoprimeFactorBasis_BigInt& CFB);


  inline const std::vector<BigInt>& FactorBase(const CoprimeFactorBasis_BigInt& CFB)  /*NOT noexcept to allow buffering*/
  { return CFB.myCoprimeBasis; }


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/NumTheory-CoprimeFactorBasis.H,v 1.4 2022/02/18 14:11:46 abbott Exp $
// $Log: NumTheory-CoprimeFactorBasis.H,v $
// Revision 1.4  2022/02/18 14:11:46  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.3  2021/02/10 19:32:00  abbott
// Summary: Added noexcept (redmine 1572)
//
// Revision 1.2  2019/09/16 17:28:53  abbott
// Summary: Changed arg names
//
// Revision 1.1  2019/03/27 14:14:27  bigatti
// (abbott) renamed GCDFreeBasis --> CoprimeFactorBasis
//
// Revision 1.1  2018/06/25 15:31:34  abbott
// Summary: GCDFreeBasis for integers
//
//

#endif
