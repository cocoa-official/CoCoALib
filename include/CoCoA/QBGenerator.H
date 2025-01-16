#ifndef CoCoA_QBGenerator_H
#define CoCoA_QBGenerator_H

//   Copyright (c)  2006  John Abbott,  Anna M. Bigatti

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


#include "CoCoA/PPMonoid.H"

#include <iosfwd>
//using std::ostream;
#include <list>
//using std::list;
#include <vector>
//using std::vector;

namespace CoCoA
{

  class QBGenerator
  {
  public:
    explicit QBGenerator(const PPMonoid& PPM);
    // default copy ctor and assignment are OK
    // default dtor is OK
    void myCornerPPIntoQB(PPMonoidElem pp);
    void myCornerPPIntoAvoidSet(ConstRefPPMonoidElem pp);
    const std::list<PPMonoidElem>& myNewCorners() const;
    const std::list<PPMonoidElem>& myCorners() const;
    const std::vector<PPMonoidElem>& myQB() const;
    const PPMonoid& myPPM() const;
    void myOutputSelf(std::ostream& out) const;
  private: // data members
    PPMonoid myPPMValue;
    std::list<PPMonoidElem> myCornerList;
    std::list<PPMonoidElem> myNewCornerList;
    std::vector<PPMonoidElem> myAvoidList;
    std::vector<PPMonoidElem> myQBList;
  };

  const PPMonoid& PPM(const QBGenerator& QBG);
  std::ostream& operator<<(std::ostream& out, const QBGenerator& QBG);

} // end of namespace CoCoA

#endif
