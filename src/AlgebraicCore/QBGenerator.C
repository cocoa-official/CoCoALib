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
