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



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/QBGenerator.H,v 1.4 2022/02/18 14:11:47 abbott Exp $
// $Log: QBGenerator.H,v $
// Revision 1.4  2022/02/18 14:11:47  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.3  2021/01/07 14:21:25  abbott
// Summary: Corrected copyright
//
// Revision 1.2  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.6  2007/02/22 17:27:02  bigatti
// -- removed useless comment
//
// Revision 1.5  2007/01/09 15:52:08  cocoa
// Changed QBGenerator to use std::vector instead of std::list for the result.
// Minor mod to configure script.
//
// Revision 1.4  2006/12/06 17:15:26  cocoa
// -- removed #include "config.H"
//
// Revision 1.3  2006/10/06 14:04:15  cocoa
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

#endif
