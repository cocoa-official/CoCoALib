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

#endif
