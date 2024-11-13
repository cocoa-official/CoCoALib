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

#endif
