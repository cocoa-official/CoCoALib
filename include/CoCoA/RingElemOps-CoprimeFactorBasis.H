#ifndef CoCoA_RingElemOps_CoprimeFactorBasis_H
#define CoCoA_RingElemOps_CoprimeFactorBasis_H

//   Copyright (c)  2016  John Abbott, and Anna M. Bigatti

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

#include "CoCoA/ring.H"

#include <vector>
//using std::vector;

namespace CoCoA
{

  class RingElem;  // fwd decl -- defined in ring.H

  class CoprimeFactorBasis_RingElem
  {
  public:
    CoprimeFactorBasis_RingElem() {}
    void myAddInfo(const RingElem& X);
    void myAddInfo(const std::vector<RingElem>& v);
    friend const std::vector<RingElem>& FactorBase(const CoprimeFactorBasis_RingElem& GFB);
  private: // implementation detail
    struct LCR
    {
      LCR(const RingElem& L, const std::vector<RingElem>& C, const RingElem& R):
          myL(L), myR(R), myC(C) {}
    public: // data members
      RingElem myL;
      RingElem myR;
      std::vector<RingElem> myC;
    };
  private: // data members
    std::vector<RingElem> myCoprimeBasis;
    void myRefineBasis(RingElem X);
    struct LCR myLCR(RingElem A, RingElem B) const;
  };


  std::ostream& operator<<(std::ostream& out, const CoprimeFactorBasis_RingElem& GFB);


  inline const std::vector<RingElem>& FactorBase(const CoprimeFactorBasis_RingElem& GFB)
  { return GFB.myCoprimeBasis; }



/////  std::vector<RingElem> GCDFreeBasis(const std::vector<RingElem>& L);
/////  RingElem RefineGCDFreeBasis(std::vector<RingElem>& basis, RingElem N);

} // end of namespace CoCoA

#endif
