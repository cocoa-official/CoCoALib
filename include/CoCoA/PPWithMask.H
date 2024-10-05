#ifndef CoCoA_PPWithMask_H
#define CoCoA_PPWithMask_H

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


#include "CoCoA/PPMonoid.H"
#include "CoCoA/DivMask.H"

#include <iosfwd>
//using std::ostream;

namespace CoCoA
{
  class PPWithMask
  {
  public:
    PPWithMask(const PPMonoid& PPM, const DivMaskRule& DMRule);
    //PPWithMask(const PPMonoid& PPM, const DivMaskRule& DMRule, const std::vector<long>& v);
    PPWithMask(ConstRefPPMonoidElem pp, const DivMaskRule& DMRule);
    PPWithMask(const PPWithMask&) = default; // to avoid compiler warnings
    ~PPWithMask() {};
    void myAssign(ConstRefPPMonoidElem pp);
    void mySwap(PPWithMask& pm);
    friend RefPPMonoidElem PP(PPWithMask& pm) { return pm.myPP; }
    friend ConstRefPPMonoidElem PP(const PPWithMask& pm) { return pm.myPP; }
    friend const DivMaskRule& DivMaskImpl(const PPWithMask& pm) { return pm.myDivMaskRule; }
    friend bool IsEqualFast(const PPWithMask& pm1, const PPWithMask& pm2); ///< suffix "Fast" because it does not check compatibility
    friend bool IsDivisibleFast(const PPWithMask& pm1, const PPWithMask& pm2);  ///< suffix "Fast" because it does not check compatibility
    PPWithMask& operator=(const PPWithMask& pm);
    PPWithMask& operator=(ConstRefPPMonoidElem pp);
    PPWithMask& operator=(const std::vector<long>& v);
    bool operator==(const PPWithMask& pm) const; ///< WARNING: it does not check compatibility
    bool operator!=(const PPWithMask& pm) const; ///< WARNING: it does not check compatibility
    friend std::ostream& operator<<(std::ostream& out, const PPWithMask& pm);
//???    friend OpenMath::OutputChannel& operator<<(OpenMath::OutputChannel& OMOut, const PPWithMask& pm);

  private:
    void myUpdateDivMask(); ///< assign the correct DivMask for myPP

    /**@name member fields of PPWithMask  */
    //@{
  private: // data members
    DivMask myDivMask;  ///< the DivMask associated to PP by myDivMaskRule
    PPMonoidElem myPP;  ///< the PP itself
    const DivMaskRule myDivMaskRule; ///< the DivMask "style": this might disappear
    //@}
  };


  std::ostream& operator<<(std::ostream& out, const PPWithMask& pm);

  //----------------------------------------------------------------------
  // Inline functions

  inline PPWithMask::PPWithMask(const PPMonoid& PPM, const DivMaskRule& DMRule):
    myPP(PPM), myDivMaskRule(DMRule)
  {}


  inline bool IsEqualFast(const PPWithMask& pm1, const PPWithMask& pm2)
  {
    return (pm1.myDivMask == pm2.myDivMask) && (pm1.myPP == pm2.myPP);
  }


  inline bool IsDivisibleFast(const PPWithMask& pm1, const PPWithMask& pm2)
  {
    // chiamare myIsDivisible per PP?
    // assert su compatibilita' DivMaskRule
    return IsSubset(pm2.myDivMask, pm1.myDivMask) && IsDivisible(pm1.myPP, pm2.myPP);
  }


  inline bool PPWithMask::operator==(const PPWithMask& pm) const
  {
    // == dovrebbe garantire la cmpatibilita'!!!!!!!
    return IsEqualFast(*this, pm);
  }


  inline bool PPWithMask::operator!=(const PPWithMask& pm) const
  {
    return !IsEqualFast(*this, pm);
  }


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/PPWithMask.H,v 1.8 2022/03/07 14:28:48 abbott Exp $
// $Log: PPWithMask.H,v $
// Revision 1.8  2022/03/07 14:28:48  abbott
// Summary: Added explicit default copy ctor to avoid g++ v.9 warning
//
// Revision 1.7  2022/02/18 14:11:47  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.6  2014/03/05 10:01:57  bigatti
// - -some improvements
//
// Revision 1.5  2010/07/14 08:59:05  bigatti
// -- added mySwap
//
// Revision 1.4  2010/07/13 12:44:12  bigatti
// -- moved documentation to doc/txt files
//
// Revision 1.3  2008/09/16 15:02:52  bigatti
// -- added operator=(pp)
// -- added operator=(v)
//
// Revision 1.2  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.5  2007/03/08 18:42:06  cocoa
// Cleaned up whitespace.
//
// Revision 1.4  2006/11/27 13:45:56  cocoa
// -- removed #include "config.H"
//
// Revision 1.3  2006/11/24 17:40:38  cocoa
// -- minor reorganization
//
// Revision 1.2  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.7  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
// Revision 1.6  2006/02/03 14:26:12  cocoa
// -- changed "DivMaskImpl" --> "DivMaskRule"
//
// Revision 1.5  2006/02/01 16:56:13  cocoa
// Added some missing assignment operators for (Ref)PPMonoidElems.
//
// Revision 1.4  2006/01/20 15:43:30  cocoa
// -- fixed: use of RefPPMonoidElem and ConstRefPPMonoidElem
//
// Revision 1.3  2006/01/18 16:15:16  cocoa
// Cleaned up DivMask considerably; everything still works,
// so I'm checking in (and then going home).
//
// Revision 1.2  2006/01/17 10:23:08  cocoa
// Updated DivMask; many consequential changes.
// A few other minor fixes.
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.2  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.2  2004/11/02 15:09:00  cocoa
// -- only minor changes to documentation
//
// Revision 1.1  2004/10/29 15:33:23  cocoa
// -- new class pairing a PPMonoidElem with a DivMask for fast
//    divisibility tests
//

#endif
