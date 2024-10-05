//   Copyright (c)  2006  John Abbott and Anna M. Bigatti

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

#include "CoCoA/QBGenerator.H"

#include "CoCoA/VectorOps.H"   // for template fn to output vector/list.
#include "CoCoA/assert.H"
#include "CoCoA/error.H"

//#include <list>
using std::list;
#include <iostream>
using std::ostream;
#include <algorithm>
using std::find_if;

namespace CoCoA
{

  namespace // anonymous namespace for file local definitions
  {

    class IsFactorOf
    {
    public:
      explicit IsFactorOf(ConstRefPPMonoidElem pp);
      bool operator()(ConstRefPPMonoidElem candidate);
    private:
      PPMonoidElem myPP;
    };

    inline IsFactorOf::IsFactorOf(ConstRefPPMonoidElem pp):
        myPP(pp)
    {}

    inline bool IsFactorOf::operator()(ConstRefPPMonoidElem candidate)
    {
      return IsDivisible(myPP, candidate);
    }

  } // end of anonymous namespace


  QBGenerator::QBGenerator(const PPMonoid& PPM):
      myPPMValue(PPM),
      myCornerList(),
      myNewCornerList(),
      myAvoidList(),
      myQBList()
  {
    myCornerList.push_back(one(PPM));
    myNewCornerList.push_back(one(PPM));
  }


  void QBGenerator::myCornerPPIntoQB(PPMonoidElem pp) // NB I do want a private copy of the arg!!!
  {
    CoCoA_ASSERT(owner(pp) == myPPMValue);
    const list<PPMonoidElem>::iterator posn = find(myCornerList.begin(), myCornerList.end(), pp);
    if (posn == myCornerList.end())
      CoCoA_THROW_ERROR("PP does not belong to the corner set", "myCornerPPIntoQB");

    // To ensure EXCEPTION SAFETY, we create a list of all possible new corners
    // which will later be filtered to remove the bad candidates.
    list<PPMonoidElem> NewCornerList;
    for (long xi=0; xi < NumIndets(myPPMValue); ++xi)
    {
      NewCornerList.push_back(pp * indet(myPPMValue, xi));
    }
    list<PPMonoidElem> NewCornerListCopy(NewCornerList);

    myQBList.push_back(pp);
    myCornerList.erase(posn); // remove pp from CornerList --> DESTROYS ORIGINAL COPY OF pp!

    // Now find out which multiples of pp (if any) are to be put into the corner lists.
    // We consider pp*x[i] in turn for each indeterminate x[i].
    list<PPMonoidElem>::iterator it = NewCornerList.begin();
    list<PPMonoidElem>::iterator it2 = NewCornerListCopy.begin();
    while (it != NewCornerList.end())
    {
      if (find_if(myAvoidList.begin(), myAvoidList.end(), IsFactorOf(*it)) == myAvoidList.end() &&
          find_if(myCornerList.begin(), myCornerList.end(), IsFactorOf(*it)) == myCornerList.end())
      { ++it; ++it2; continue; }
      it = NewCornerList.erase(it);
      it2 = NewCornerListCopy.erase(it2);
    }

    // Sort the new corner elements,
    // and merge a copy of them into the full set of corner elements.
    NewCornerList.sort();
    NewCornerListCopy.sort(); //??? useful ???
    myCornerList.merge(NewCornerListCopy); // merge moves the PPs into CornerList
    myNewCornerList.swap(NewCornerList);
  }


  void QBGenerator::myCornerPPIntoAvoidSet(ConstRefPPMonoidElem pp)
  {
    CoCoA_ASSERT(owner(pp) == myPPMValue);
    const list<PPMonoidElem>::iterator it = find(myCornerList.begin(), myCornerList.end(), pp);
    if (it == myCornerList.end())
      CoCoA_THROW_ERROR("PP does not belong to the corner set", "myCornerPPIntoAvoidSet");
    myAvoidList.push_back(pp); // Do this *before* deleting pp from myCornerList as pp may alias the element we delete!
    myCornerList.erase(it);
    myNewCornerList.clear();
  }


  const std::list<PPMonoidElem>& QBGenerator::myNewCorners() const
  {
    return myNewCornerList;
  }


  const std::list<PPMonoidElem>& QBGenerator::myCorners() const
  {
    return myCornerList;
  }


  const std::vector<PPMonoidElem>& QBGenerator::myQB() const
  {
    return myQBList;
  }


  const PPMonoid& QBGenerator::myPPM() const
  {
    return myPPMValue;
  }


  void QBGenerator::myOutputSelf(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "QBGenerator("
        << "QB=" << myQBList
        << ", corners=" << myCornerList
        << ", NewCorners=" << myNewCornerList
        << ", avoid=" << myAvoidList
        << ")";
  }


  const PPMonoid& PPM(const QBGenerator& QBG)
  {
    return QBG.myPPM();
  }


  std::ostream& operator<<(std::ostream& out, const QBGenerator& QBG)
  {
    if (!out) return out;  // short-cut for bad ostreams
    QBG.myOutputSelf(out);
    return out;
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/QBGenerator.C,v 1.14 2022/02/18 14:11:56 abbott Exp $
// $Log: QBGenerator.C,v $
// Revision 1.14  2022/02/18 14:11:56  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.13  2021/01/07 15:16:52  abbott
// Summary: Corrected copyright
//
// Revision 1.12  2020/06/17 15:49:25  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.11  2020/02/11 16:56:41  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.10  2020/02/11 16:12:18  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.9  2018/05/17 15:39:44  bigatti
// -- renamed VectorOperations --> VectorOps
//
// Revision 1.8  2016/11/11 14:15:33  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.7  2015/04/16 16:39:26  abbott
// Summary: Rewrote myCornerPPIntoQB: wastes fewer copies, still exception safe
// Author: JAA
//
// Revision 1.6  2015/04/15 13:52:34  abbott
// Changed impl of myCornerPPIntoQB; avoids copying lists of PPs.
//
// Revision 1.5  2014/07/31 14:45:18  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.4  2013/01/21 13:29:53  abbott
// Changed impl so that it is exception safe.
//
// Revision 1.3  2011/03/10 16:39:34  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.2  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.6  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.5  2007/01/09 15:52:08  cocoa
// Changed QBGenerator to use std::vector instead of std::list for the result.
// Minor mod to configure script.
//
// Revision 1.4  2006/11/24 17:06:10  cocoa
// -- reorganized includes of header files
//
// Revision 1.3  2006/10/06 14:04:14  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.2  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.2  2006/04/27 13:45:30  cocoa
// Changed name of NewIdentityRingHom to NewIdentityHom.
// Changed name of member functions which print out their own object
// into myOutputSelf (to distinguish from "transitive" myOutput fns).
//
// Revision 1.1  2006/04/21 15:03:23  cocoa
// New code for Buchberger-Moeller and variants.
//
//
