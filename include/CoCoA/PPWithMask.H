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

#endif
