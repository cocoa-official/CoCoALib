#ifndef CoCoA_DivMask_H
#define CoCoA_DivMask_H

//   Copyright (c)  2005-2007  John Abbott,  Anna M. Bigatti

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


#include "CoCoA/config.H"
#include "CoCoA/SmartPtrIRC.H"

#include <bitset>
//using std::bitset;
#include <iosfwd>
//using std::ostream;

namespace CoCoA
{

  //-- class DivMask -----------------------------------------
  /**

  \brief Implementation of DivMask for fast divisibility test

  DivMask implements a generalization of the present/absent bitset of a PP
  for fast divisibility tests.
  It is represented as a bitset whose size must be defined at compile-time
  for efficiency (its main reason for existence!) and it must satisfy:

      pp1 | pp2  ==> dm1 | dm2      ( lcm(pp1,pp2)=pp2  ==> dm1 && dm2 = dm2 )

  See DivMaskRuleBase.
  */
  //-----------------------------------------------------------------

  class DivMaskRuleBase; // forward declaration for use in friend declaration

  class DivMask
  {
  public:
    static constexpr long ourMaskWidth = 64;  // ourMaskWidth must be a power of 2.
    typedef std::bitset<ourMaskWidth> mask_t;
  public:
    // Default ctor, copy ctor, assignment and dtor are OK.
    friend const mask_t bits(const DivMask& dm) noexcept;
  private: // data members
    friend class DivMaskRuleBase; // see mem.fn. DivMaskRuleBase::myBits
    mask_t myBits;
  };

  std::ostream& operator<<(std::ostream& out, const DivMask& dm);


  class DivMaskRule
  {
  public:
    explicit DivMaskRule(const DivMaskRuleBase* DMRPtr): mySmartPtr(DMRPtr) {}
    // assignment disabled because SmartPtrIRC has no assignment
    const DivMaskRuleBase* operator->() const noexcept { return mySmartPtr.operator->(); }  ///< Allow const member fns to be called.
    bool operator==(const DivMaskRule& DMR) const noexcept { return mySmartPtr==DMR.mySmartPtr; }
  private: // data members
    SmartPtrIRC<const DivMaskRuleBase> mySmartPtr;
  };


  class DivMaskRuleBase: protected IntrusiveReferenceCount
  {
    friend class SmartPtrIRC<const DivMaskRuleBase>; // morally "friend DivMaskRule", so it can alter reference count
  protected:
    DivMaskRuleBase(): IntrusiveReferenceCount() {};
    // Copy ctor and assignment disabled (see below).
    virtual ~DivMaskRuleBase() {};
  public:
    virtual void myAssignFromExpv(DivMask& dm, const SmallExponent_t* expv, long NumIndets) const = 0;
//???    virtual void myAssignFromExpv(DivMask& dm, const std::vector<long>& expv) const = 0;
//???    virtual void myAssignFromExpv(DivMask& dm, const std::vector<BigInt>& EXPV) const = 0;
    virtual void myOutputSelf(std::ostream& out) const = 0;
//????    virtual void myOutputSelf_OM(OpenMath::OutputChannel& OMOut) const = 0;
  public: // Disable copy ctor and assignment.
    explicit DivMaskRuleBase(const DivMaskRuleBase&) = delete;
    DivMaskRuleBase& operator=(const DivMaskRuleBase&) = delete;
  protected:
    DivMask::mask_t& myBits(DivMask& dm) const noexcept; ///< WRITE ACCESS to the bitset for derived concrete classes
  };


  std::ostream& operator<<(std::ostream& out, const DivMaskRule& DMR);

  DivMaskRule NewDivMaskNull();
  DivMaskRule NewDivMaskSingleBit();
  DivMaskRule NewDivMaskSingleBitWrap();
  DivMaskRule NewDivMaskEvenPowers();
  DivMaskRule NewDivMaskHashing();


  //----------------------------------------------------------------------
  // Inline functions -- as always, be wary of changing the order!

  inline const DivMask::mask_t bits(const DivMask& dm) noexcept
  {
    return dm.myBits;
  }

  inline bool operator==(const DivMask& dm1, const DivMask& dm2) noexcept
  {
    return bits(dm1) == bits(dm2);
  }

  inline bool operator!=(const DivMask& dm1, const DivMask& dm2) noexcept
  {
    return !(dm1 == dm2);
  }


  inline bool IsSubset(const DivMask& dm1, const DivMask& dm2) noexcept // is dm1 subset of dm2?
  {
    return (bits(dm1) & bits(dm2)) == bits(dm1);
  }


  // This is like bits(dm) but gives write access to the bitset.
  inline DivMask::mask_t& DivMaskRuleBase::myBits(DivMask& dm) const noexcept
  {
    return dm.myBits;
  }


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/DivMask.H,v 1.13 2023/07/07 09:24:57 abbott Exp $
// $Log: DivMask.H,v $
// Revision 1.13  2023/07/07 09:24:57  abbott
// Summary: Maskwidth reset to 64
//
// Revision 1.12  2023/05/11 19:33:17  abbott
// Summary: INcr MaskWidth to 128 (why not 64?)
//
// Revision 1.11  2022/02/18 14:11:45  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.10  2022/02/08 20:18:53  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.9  2021/10/30 17:13:07  abbott
// Summary: Used keyword delete (redmine 1627)
//
// Revision 1.8  2021/02/10 19:32:00  abbott
// Summary: Added noexcept (redmine 1572)
//
// Revision 1.7  2021/01/07 14:12:10  abbott
// Summary: Corrected copyright
//
// Revision 1.6  2012/02/06 18:29:28  bigatti
// -- changed MaskWidth --> ourMaskWidth
// -- removed doxygen include
//
// Revision 1.5  2011/08/14 15:52:18  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.4  2011/03/14 10:24:03  abbott
// Changed size_t into long.
//
// Revision 1.3  2009/09/24 14:10:04  abbott
// Moved a forward declaration.  Added "explicit" keyword.
//
// Revision 1.2  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.8  2007/03/08 18:42:06  cocoa
// Cleaned up whitespace.
//
// Revision 1.7  2006/11/29 17:26:58  cocoa
// -- just shortened the log
//
// Revision 1.6  2006/11/23 17:06:55  cocoa
// -- changed: DivMaskRule is now a class (instead of typedef)
//
// Revision 1.5  2006/11/16 16:08:21  cocoa
// -- changed: DivMaskRuleBase ctor and dtor are now protected
//
// Revision 1.4  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.3  2006/10/06 14:04:16  cocoa
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
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//

#endif
