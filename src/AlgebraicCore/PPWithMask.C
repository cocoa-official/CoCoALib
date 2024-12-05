//   Copyright (c)  2005  Anna Bigatti

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


#include "CoCoA/PPWithMask.H"
#include "CoCoA/assert.H"
#include "CoCoA/utils.H"

#include <iostream>
using std::ostream;

namespace CoCoA
{

//   PPWithMask::PPWithMask(const PPMonoid& PPM, const DivMaskRule& dm, const std::vector<long>& v):
//     myPP(PPM, v), myDivMaskRule(dm)
//   {
//     myUpdateDivMask();
//   }


  PPWithMask::PPWithMask(ConstRefPPMonoidElem pp, const DivMaskRule& dm):
    myPP(pp), myDivMaskRule(dm)
  {
    myUpdateDivMask();
  }


  PPWithMask& PPWithMask::operator=(const PPWithMask& pm)
  {
    CoCoA_ASSERT( myDivMaskRule == pm.myDivMaskRule );
    myPP = pm.myPP;
    myDivMask = pm.myDivMask;
    return *this;
  }


  PPWithMask& PPWithMask::operator=(ConstRefPPMonoidElem pp)
  {
    CoCoA_ASSERT( owner(myPP) == owner(pp) );  // or not?  (using expv)
    PPWithMask res(pp, myDivMaskRule);
    swap(myPP, res.myPP);
    std::swap(myDivMask, res.myDivMask);
    return *this;
  }


  PPWithMask& PPWithMask::operator=(const std::vector<long>& v)
  {
    CoCoA_ASSERT( NumIndets(owner(myPP)) == len(v) );
    PPWithMask res(PPMonoidElem(owner(myPP), v), myDivMaskRule);
    swap(myPP, res.myPP);
    std::swap(myDivMask, res.myDivMask);
    return *this;
  }


  void PPWithMask::myAssign(ConstRefPPMonoidElem pp)
  {
    myPP = pp;
    myUpdateDivMask();
  }


  void PPWithMask::mySwap(PPWithMask& pm)
  {
    swap(myPP, pm.myPP);
    std::swap(myDivMask, pm.myDivMask);
  }


  void PPWithMask::myUpdateDivMask()
  {
    owner(myPP)->myComputeDivMask(myDivMask, myDivMaskRule, raw(myPP));
  }


  std::ostream& operator<<(std::ostream& out, const PPWithMask& pm)
  {
    if (!out) return out;  // short-cut for bad ostreams
    out << "PPWithMask(" << pm.myPP << ", " << pm.myDivMask << ", " << pm.myDivMaskRule << ")";
    return out;
  }

}  // end of namespace CoCoA
