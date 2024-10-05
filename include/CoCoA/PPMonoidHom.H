#ifndef CoCoA_PPMonoidHom_H
#define CoCoA_PPMonoidHom_H

//   Copyright (c)  2010  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/SmartPtrIRC.H"

#include <iosfwd>
//using std::ostream;


namespace CoCoA
{

  /*-----------------------------------------------------------------*/
  /** \include PPMonoidHom.txt  */
  /*-----------------------------------------------------------------*/

  class PPMonoidHom; // fwd decl (see below)

  class PPMonoidHomBase: protected IntrusiveReferenceCount
  {
    friend class SmartPtrIRC<const PPMonoidHomBase>; // Morally "friend PPMonoidHom", so it can alter reference count.
  protected:
    PPMonoidHomBase(const PPMonoid& domain, const PPMonoid& codomain): myDomain(domain), myCodomain(codomain) {}
    virtual ~PPMonoidHomBase() {}
  public: // Disable copy ctor and assignment.
    PPMonoidHomBase(const PPMonoidHomBase&) = delete;
    PPMonoidHomBase& operator=(const PPMonoidHomBase&) = delete;
  public:
    friend const PPMonoid& domain(const PPMonoidHom& phi);
    friend const PPMonoid& codomain(const PPMonoidHom& phi);
    virtual void myApply(PPMonoidElemRawPtr image, PPMonoidElemConstRawPtr arg) const = 0;
    void myOutputSelf(std::ostream& out) const;
//???    virtual void myOutputSelf_OM(OpenMath::OutputChannel& OMOut) const = 0;
  private: // Implementation details.
    virtual void myOutputSelfDetails(std::ostream& out) const; ///< called only by myOutputSelf

  protected: // data members
    const PPMonoid myDomain;
    const PPMonoid myCodomain;
  };



  //---------------------------------------------------------------------------

  class PPMonoidHom
  {
  public:
    explicit PPMonoidHom(const PPMonoidHomBase* phi): mySmartPtr(phi) {}
    // Default copy ctor works fine.
    // Assignment disabled because SmartPtrIRC has no assignment.
    // Default dtor works fine.
    PPMonoidElem operator()(ConstRefPPMonoidElem x) const;      ///< apply hom to an element
    PPMonoidHom operator()(const PPMonoidHom&) const;           ///< composition of two PPMonoidHoms
    const PPMonoidHomBase* operator->() const;
//     const PPMonoidHomBase* myHomPtr() const;
  private: // data member
    const SmartPtrIRC<const PPMonoidHomBase> mySmartPtr;
  };


  const PPMonoid& domain(const PPMonoidHom& phi);    ///< inline accessor fn
  const PPMonoid& codomain(const PPMonoidHom& phi);  ///< inline accessor fn
  std::ostream& operator<<(std::ostream& out, const PPMonoidHom& phi);

  PPMonoidHom IdentityHom(const PPMonoid& PPM);
  PPMonoidHom GeneralHom(const PPMonoid& PPM, const std::vector<PPMonoidElem>& images);
///???  PPMonoidHom PermutationHom(const PPMonoid& domain, const PPMonoid& codomain, const std::vector<int>& );
  PPMonoidHom RestrictionHom(const PPMonoid& PPM, const std::vector<long>& IndetIndices);

  //----------------------------------------------------------------------
  // inline functions

  inline const PPMonoidHomBase* PPMonoidHom::operator->() const
  { return mySmartPtr.operator->(); }

  inline const PPMonoid& domain(const PPMonoidHom& phi)
  { return phi->myDomain; }

  inline const PPMonoid& codomain(const PPMonoidHom& phi)
  { return phi->myCodomain; }


}  // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/PPMonoidHom.H,v 1.7 2022/02/18 14:11:47 abbott Exp $
// $Log: PPMonoidHom.H,v $
// Revision 1.7  2022/02/18 14:11:47  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.6  2022/02/08 20:18:53  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.5  2021/10/30 17:11:20  abbott
// Summary: Used keyword delete (redmine 1627)
//
// Revision 1.4  2021/01/07 14:20:46  abbott
// Summary: Corrected copyright
//
// Revision 1.3  2012/02/10 17:08:06  abbott
// Added new pseudo-ctor for RestrictionHom (& related impl class).
//
// Revision 1.2  2011/03/10 16:39:35  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.1  2010/07/09 17:03:42  abbott
// First simple implementation of PPMonoid homs.
//
//

#endif
