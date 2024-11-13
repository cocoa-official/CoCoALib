#ifndef CoCoA_SugarDegree_H
#define CoCoA_SugarDegree_H

//   Copyright (c)  2009  Anna M. Bigatti

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


// #include "CoCoA/SmartPtrIRC.H"

#include <iosfwd>
//using std::ostream;
#include <cstddef>
//using std::size_t;
#include <memory>
// using std::unique_ptr;

namespace CoCoA
{

  class ConstRefPPMonoidElem; // fwd decl -- defined in PPMonoid.H
  class RingElemAlias;        // fwd decl -- defined in ring.H
  typedef const RingElemAlias& ConstRefRingElem;     // fwd decl -- defined in ring.H
  class GPoly;                // fwd decl -- defined in TmpGPoly.H
  class ReductionCog;         // fwd decl -- defined in ReductionCog.H
  class degree;               // fwd decl -- defined in degree.H


  enum UninitializedMarker { uninitialized }; // just used for sugar ctor

  class SugarDegreeBase //: protected IntrusiveReferenceCount
  {
  public:
    virtual ~SugarDegreeBase() {};
    virtual SugarDegreeBase* myClone() const =0;
    ///@name Warning! Throw an error if the wrong type of value is asked!
    //@{
    virtual const degree& myWSugar() const =0; ///< only if impl stores this value
    virtual long myStdSugar() const =0; ///< only if impl stores this value
    //@}
    virtual void myMul(ConstRefPPMonoidElem pp) =0; ///< sugar after multiplying by pp
    virtual void myUpdate(ReductionCog F, const GPoly& g); ///< sugar after reducing F by g
    virtual void myUpdate(ConstRefPPMonoidElem CofactorPP, const GPoly& g) =0; ///< sugar after adding pp*g
    virtual int myCmp(const SugarDegreeBase& s) const =0;  // this <=> s ? <0,=0,>0
    virtual void myOutput(std::ostream& out) const =0;
  };


  class SugarDegree
  {
  public:
    SugarDegree(UninitializedMarker) {} // default: myPtr(nullptr)
    explicit SugarDegree(SugarDegreeBase* ptr): myPtr(ptr) {}
    SugarDegree(const SugarDegree& sd); ///< unique_ptr: default copy is wrong!
    SugarDegree& operator=(const SugarDegree& rhs); ///< unique_ptr: makes a copy!
    SugarDegreeBase* operator->() { return myPtr.operator->(); }  ///< Allow non-const member fns to be called.
    const SugarDegreeBase* operator->() const { return myPtr.operator->(); }  ///< Allow const member fns to be called.
    bool operator!=(const SugarDegree& s) const { return myPtr->myCmp(*(s.myPtr.get())) != 0; }
    bool operator<(const SugarDegree& s) const { return myPtr->myCmp(*(s.myPtr.get())) < 0; }
    friend bool IsInitialized(const SugarDegree& sd) { return sd.myPtr.get()!=nullptr; }
    friend int cmp(const SugarDegree& s1, const SugarDegree& s2) { return s1->myCmp(*(s2.myPtr.get())); }
  private: // data member
    std::unique_ptr<SugarDegreeBase> myPtr;
  };


  std::ostream& operator<<(std::ostream& out, const SugarDegree& s);

  ///@name pseudo-constructors
  //@{
  SugarDegree NewStdSugar(ConstRefRingElem f);
  SugarDegree NewStdSugarNoIdx(ConstRefRingElem f, long PosIndet);
  SugarDegree NewStdSugarSat(ConstRefRingElem f);
  SugarDegree NewStdSugarNoIdxSat(ConstRefRingElem f, long PosIndet);

  SugarDegree NewWSugar(ConstRefRingElem f);
  SugarDegree NewWDeg1CompTmp(ConstRefRingElem f); ///< temporary: only for testing
  SugarDegree NewWSugarConst(ConstRefRingElem f);  ///< stays constant in myUpdate
  SugarDegree NewWSugarSat(ConstRefRingElem f);
  //@}


  //------- inline functions ----------------------------

  // none

} // end of namespace CoCoA

#endif
