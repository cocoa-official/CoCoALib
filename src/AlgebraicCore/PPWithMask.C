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

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/PPWithMask.C,v 1.10 2022/02/18 14:11:56 abbott Exp $
// $Log: PPWithMask.C,v $
// Revision 1.10  2022/02/18 14:11:56  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.9  2019/12/11 14:51:18  abbott
// Summary: Removed 1 line of cruft
//
// Revision 1.8  2016/11/11 14:15:33  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.7  2011/03/10 17:28:56  bigatti
// -- using len instead of size in assert
//
// Revision 1.6  2010/12/26 13:04:37  abbott
// Changed "GlobalXXXput" into corresponding std C++ stream
// (even in commented out code).
//
// Revision 1.5  2010/07/14 08:59:05  bigatti
// -- added mySwap
//
// Revision 1.4  2008/09/16 15:02:52  bigatti
// -- added operator=(pp)
// -- added operator=(v)
//
// Revision 1.3  2008/05/27 16:20:10  bigatti
// *** empty log message ***
//
