#ifndef CoCoA_ReductionCog_H
#define CoCoA_ReductionCog_H

//   Copyright (c)  2006  Anna M. Bigatti

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
#include "CoCoA/ring.H"
#include "CoCoA/SmartPtrIRC.H"

#include <iosfwd>
//using std::ostream;

namespace CoCoA
{

  //JAA  class RingElem;         // forward declaration -- defined in ring.H
  class SparsePolyRing;   // forward declaration -- defined in SparsePolyRing.H

  class ReductionCogBase; // forward declaration for SmartPtrIRC
  class ReductionCog
  {
  public:
    explicit ReductionCog(ReductionCogBase* ptr): mySmartPtr(ptr) {}
    // assignment disabled because SmartPtrIRC has no assignment
    ReductionCogBase* operator->() { return mySmartPtr.operator->(); }  ///< Allow member fns to be called.
    const ReductionCogBase* operator->() const { return mySmartPtr.operator->(); }  ///< Allow const member fns to be called.
    bool operator==(const ReductionCog& f) const { return mySmartPtr==f.mySmartPtr; }
  private: // data members
    SmartPtrIRC<ReductionCogBase> mySmartPtr;
  };


  ConstRefPPMonoidElem ActiveLPP(const ReductionCog& F);
  bool IsActiveZero(const ReductionCog& F);
  std::ostream& operator<<(std::ostream& out, const ReductionCog& F);

  class ReductionCogBase: protected IntrusiveReferenceCount
  {
    friend class SmartPtrIRC<ReductionCogBase>; // morally "friend ReductionCog", so it can alter reference count
  protected:
    ReductionCogBase(): IntrusiveReferenceCount() {}
  public:
    virtual ~ReductionCogBase() {};
    virtual void myMoveToNextLM() =0;
    virtual void myReduce(ConstRefRingElem reducer, long RedLen=0) =0;
    virtual std::ostream& myOutput(std::ostream& out) const =0;
    virtual void myAssignReset(RingElem& f) =0; // almost no-throw guarantee (?)
    virtual void myAssignReset(RingElem& f, long fLen) =0; // almost no-throw guarantee (?)
    virtual void myRelease(RingElem& f) =0; // almost no-throw guarantee (?)

    friend ConstRefPPMonoidElem ActiveLPP(const ReductionCog& F);
    friend bool IsActiveZero(const ReductionCog& F);

  protected:
    virtual bool IamActiveZero() const =0;
    virtual ConstRefPPMonoidElem myActiveLPP() const =0;
  };


  ReductionCog NewRedCogPolyField(const SparsePolyRing& P);
  ReductionCog NewRedCogPolyGCD(const SparsePolyRing& P);
  ReductionCog NewRedCogGeobucketField(const SparsePolyRing& P);
  ReductionCog NewRedCogGeobucketGCD(const SparsePolyRing& P);

  //------- inline functions ----------------------------

  inline ConstRefPPMonoidElem ActiveLPP(const ReductionCog& F)
  { return F->myActiveLPP(); }

  inline bool IsActiveZero(const ReductionCog& F)
  { return F->IamActiveZero(); }


} // end of namespace CoCoA

#endif
