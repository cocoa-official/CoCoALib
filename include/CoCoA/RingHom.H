#ifndef CoCoA_RingHom_H
#define CoCoA_RingHom_H

//   Copyright (c)  2005  John Abbott, and Anna M. Bigatti

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


#include "CoCoA/ring.H"
#include "CoCoA/SmartPtrIRC.H"

#include<vector>
//using std::vector;
#include <iosfwd>
//using std::ostream;


namespace CoCoA
{

  class ideal; // Fwd decl for IsInKer, defined in ideal.H.
  class ConstMatrixView; // Fwd decl for phi(M), defined in MatrixView.H.
  class matrix; // Fwd decl for phi(M), defined in matrix.H.
  class PartialRingHom; // Fwd decl

  class RingHomBase: protected IntrusiveReferenceCount
  {
    friend class SmartPtrIRC<const RingHomBase>; // Morally "friend RingHom", so it can alter reference count.
  protected:
    RingHomBase(const ring& domain, const ring& codomain): myDomain(domain), myCodomain(codomain) {}
    virtual ~RingHomBase() {}
  public: // Disable copy ctor and assignment.
    RingHomBase(const RingHomBase&) = delete;
    RingHomBase& operator=(const RingHomBase&) = delete;
  public:
    friend const ring& domain(const RingHom& phi);
    friend const ring& codomain(const RingHom& phi);
    friend const ring& domain(const PartialRingHom& phi);
    friend const ring& codomain(const PartialRingHom& phi);
    virtual void myApply(RingElemRawPtr image, RingElemConstRawPtr arg) const = 0;
    virtual bool IamPartial() const = 0; ///??? { return true; } ///??? default defn or pure virtual????
    void myOutputSelf(std::ostream& out) const;
//???    virtual void myOutputSelf_OM(OpenMath::OutputChannel& OMOut) const = 0;
///???  private: // Implementation details.
    virtual void myOutputSelfDetails(std::ostream& out) const; ///< called only by myOutputSelf

  protected: // const data members
    const ring myDomain;
    const ring myCodomain;
  };



  //---------------------------------------------------------------------------

  class PartialRingHom
  {
  public:
    explicit PartialRingHom(const RingHomBase* phi): mySmartPtr(phi) {}
    // Default copy ctor works fine.
    // Default assignment works fine.
    // Default dtor works fine.
    RingElem operator()(ConstRefRingElem x) const;      ///< apply hom to an element
////????    RingElem operator()(const RingElem& x) const;       ///< apply hom to an element
    RingElem operator()(const MachineInt& n) const;     ///< map n into domain then apply hom
    RingElem operator()(const BigInt& N) const;         ///< map N into domain then apply hom
    RingElem operator()(const BigRat& q) const;         ///< map q into domain then apply hom
    // apply phi to other types:
    //  matrix apply(RingHom phi, ConstMatrixView M);  // result is phi(M)
    matrix operator()(ConstMatrixView M) const;   // result is phi(M)
    //-------------------------
    PartialRingHom operator()(const PartialRingHom&) const;    ///< composition of two RingHoms
    PartialRingHom operator()(const RingHom&) const;           ///< composition of two RingHoms
    const RingHomBase* operator->() const;
//     const RingHomBase* myHomPtr() const;
  protected: // data member
    SmartPtrIRC<const RingHomBase> mySmartPtr;
  };


  class RingHom: public PartialRingHom
  {
  public:
    explicit RingHom(const RingHomBase* phi): PartialRingHom(phi) { /*******CoCoA_ASSERT(jaa() && !phi->IamPartial());********/ }  ///// EMERGENCY COMMENTED OUT CODE
    // Default copy ctor works fine.
    // Default assignment works fine.
    // Default dtor works fine.
    RingElem operator()(ConstRefRingElem x) const;      ///< apply hom to an element
    RingElem operator()(const MachineInt& n) const;     ///< map n into domain then apply hom
    RingElem operator()(const BigInt& N) const;         ///< map N into domain then apply hom
    RingElem operator()(const BigRat& q) const;         ///< map q into domain then apply hom
    // // apply phi to other types:
    // //  matrix apply(RingHom phi, ConstMatrixView M);  // result is phi(M)
    matrix operator()(const ConstMatrixView M) const;   // result is phi(M)
    std::vector<RingElem> operator()(const std::vector<RingElem>& v) const;  ///< map each elem of v to produce new vector
    // //-------------------------
    RingHom operator()(const RingHom&) const;           ///< composition of two RingHoms
    PartialRingHom operator()(const PartialRingHom&) const; ///< composition of two RingHoms
//     const RingHomBase* operator->() const;
//     const RingHomBase* myHomPtr() const;
//   private: // data member
//     SmartPtrIRC<const RingHomBase> mySmartPtr;
  };


  const ring& domain(const RingHom& phi);    ///< inline accessor fn
  const ring& codomain(const RingHom& phi);  ///< inline accessor fn
  const ring& domain(const PartialRingHom& phi);    ///< inline accessor fn
  const ring& codomain(const PartialRingHom& phi);  ///< inline accessor fn
  std::ostream& operator<<(std::ostream& out, const RingHom& phi);
  bool IsInKer(const ideal& I, const RingHom& phi); ///< is I in ker(phi)?
  RingHom IdentityHom(const ring& R);
  PartialRingHom sequential(const PartialRingHom& phi, const PartialRingHom& theta);
  bool ImageLiesInSubfield(const RingHom& phi);


  //---------------------------------------------------------------------------
  // Two further abstract base classes; the only change wrt RingHomBase is
  // the message printed out in myOutputSelfDetails.

  class RingHomEmbeddingBase: public RingHomBase
  {
  protected:
    RingHomEmbeddingBase(const ring& domain, const ring& codomain): RingHomBase(domain, codomain) {}
    // Copy ctor and assignment disabled in RingHomBase
    // Default dtor is fine.
  private: // Implementation details.
    virtual void myOutputSelfDetails(std::ostream& out) const; ///< called only by myOutputSelf
    virtual bool IamPartial() const { return false; }
  };


  class RingHomInducedBase: public RingHomBase
  {
  protected:
    RingHomInducedBase(const ring& NewDomain, const RingHom& InducingHom);
    // Copy ctor and assignment disabled in RingHomBase
    // Default dtor is fine.
  private: // Implementation details.
    virtual void myOutputSelfDetails(std::ostream& out) const; ///< called only by myOutputSelf
  protected: // data member
    const RingHom myInducingHom;
  };


  //----------------------------------------------------------------------
  // inline functions

  inline const RingHomBase* PartialRingHom::operator->() const
  { return mySmartPtr.operator->(); }

  inline const ring& domain(const RingHom& phi)
  { return phi->myDomain; }

  inline const ring& codomain(const RingHom& phi)
  { return phi->myCodomain; }

  inline const ring& domain(const PartialRingHom& phi)
  { return phi->myDomain; }

  inline const ring& codomain(const PartialRingHom& phi)
  { return phi->myCodomain; }

  inline bool IsPartial(const RingHom& phi)
  { CoCoA_ASSERT(!phi->IamPartial()); (void)phi/*avoid unused var warning*/; return false; }

  inline bool IsPartial(const PartialRingHom& phi)
  { return phi->IamPartial(); }


}  // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/RingHom.H,v 1.19 2022/03/07 14:29:26 abbott Exp $
// $Log: RingHom.H,v $
// Revision 1.19  2022/03/07 14:29:26  abbott
// Summary: Corrected two false comments (about assignment being disabled)
//
// Revision 1.18  2022/02/18 14:11:48  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.17  2022/02/08 20:18:53  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.16  2021/10/30 17:15:37  abbott
// Summary: Used keyword delete (redmine 1627)
//
// Revision 1.15  2021/08/02 07:55:19  abbott
// Summary: Added prototype phi(vec); redmine 1467
//
// Revision 1.14  2021/07/30 15:28:37  bigatti
// Summary: added phi(MAT)  (redmine #1598)
//
// Revision 1.13  2021/01/07 14:44:46  abbott
// Summary: Corrected copyright
//
// Revision 1.12  2018/06/25 15:32:02  abbott
// Summary: Removed debugging kludge
//
// Revision 1.11  2013/02/21 14:14:42  abbott
// First attempt at implementing PartialRingHom -- some problems remain!!
//
// Revision 1.10  2012/06/19 15:44:01  abbott
// Realigned comments.
//
// Revision 1.9  2011/11/09 13:50:01  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.8  2011/08/24 10:21:09  bigatti
// -- renamed QQ --> BigRat
//
// Revision 1.7  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.6  2011/03/14 10:27:02  abbott
// Removed superfluous #include directive.
//
// Revision 1.5  2011/02/23 15:00:03  bigatti
// -- RingHom can now be applied also to BigRat
//
// Revision 1.4  2010/10/06 14:09:08  abbott
// Removed const qualifier so that RingHom values can be assigned.
//
// Revision 1.3  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.2  2007/10/30 17:14:11  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.9  2007/03/05 21:06:07  cocoa
// New names for homomorphism pseudo-ctors: removed the "New" prefix.
//
// Revision 1.8  2007/03/03 15:24:14  cocoa
// Improved a comment.
//
// Revision 1.7  2007/01/13 14:14:34  cocoa
// Overhaul of RingHom code: it nows uses SmartPtrIRC, and printing is more logical.
// Have not yet updated the documentation.
//
// Revision 1.6  2007/01/11 14:06:15  cocoa
// -- added prefix "raw" to RawPtr arguments names
//
// Revision 1.5  2006/12/06 17:21:55  cocoa
// -- removed #include "config.H"
//
// Revision 1.4  2006/11/29 16:52:07  cocoa
// -- doxygen style comments
//
// Revision 1.3  2006/11/08 16:21:59  cocoa
// Structural cleaning of RingHom; many consequential changes.
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
// Revision 1.4  2006/04/27 13:45:30  cocoa
// Changed name of NewIdentityRingHom to NewIdentityHom.
// Changed name of member functions which print out their own object
// into myOutputSelf (to distinguish from "transitive" myOutput fns).
//
// Revision 1.3  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.2  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.4  2005/04/20 15:40:48  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.3  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.2  2005/02/11 14:15:20  cocoa
// New style ring elements and references to ring elements;
// I hope I have finally got it right!
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.7  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.6  2004/11/11 13:04:30  cocoa
// -- added \include *.txt for doxygen
//
// Revision 1.5  2004/07/20 15:04:06  cocoa
// The next step in the new "ring element" conversion process:
// handling the case of creating a "const RefRingElem" object
// (since C++ refuses to do this properly itself).
//
// Revision 1.4  2004/05/27 16:14:02  cocoa
// Minor revision for new coding conventions.
//
// Revision 1.3  2004/05/24 15:52:13  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.2  2003/10/09 12:16:38  cocoa
// New coding convention for rings.
//
// Revision 1.2  2003/06/23 16:59:19  abbott
// Minor cleaning prior to public release.
//
// Revision 1.1  2003/05/14 17:16:35  abbott
// Initial revision
//

#endif
