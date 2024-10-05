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


//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/SugarDegree.H,v 1.19 2022/02/18 14:11:50 abbott Exp $
// $Log: SugarDegree.H,v $
// Revision 1.19  2022/02/18 14:11:50  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.18  2020/02/11 17:16:17  abbott
// Summary: All mem fns like myOutput now return void (redmine 969)
//
// Revision 1.17  2019/03/18 17:56:49  abbott
// Summary: Changed 0 into nullptr where appropriate
//
// Revision 1.16  2019/03/04 10:27:28  abbott
// Summary: Changed auto_ptr into unqiue_ptr
//
// Revision 1.15  2012/10/24 12:07:14  abbott
// Revised fwd decl of ConstRefRingElem.
//
// Revision 1.14  2011/08/27 21:46:56  abbott
// Added friend fn "cmp"
//
// Revision 1.13  2010/06/10 08:14:07  bigatti
// -- commented out include for SmartPtrIRC.H
//
// Revision 1.12  2010/05/21 16:11:37  bigatti
// -- removed CheckSugar
//
// Revision 1.11  2010/05/14 09:53:09  bigatti
// -- removed empty ctor for SugarDegree
// -- added marker for SugarDegree(uninitialized)
// -- SugarDegree for GBasis input is initialized by myPrepareGBasis
//
// Revision 1.10  2010/03/30 16:05:33  bigatti
// -- added comments
//
// Revision 1.9  2009/09/25 12:36:48  bigatti
// -- added CheckSugar (temporary)
//
// Revision 1.8  2009/07/20 14:28:45  bigatti
// -- modified interface for pseudo constructors (now with RingElem argument)
//
// Revision 1.7  2009/04/24 16:34:11  bigatti
// -- added some myDeg and pseudo-ctor
//
// Revision 1.6  2009/03/20 14:01:14  bigatti
// -- minor changes (names)
//
// Revision 1.5  2009/03/18 16:37:13  bigatti
// -- "almost-final" cleanup, tested against some current sugar definitions
//
// Revision 1.4  2009/03/18 15:11:34  abbott
// Minor cosmetic changes.
//
// Revision 1.3  2009/03/16 16:34:42  bigatti
// -- changed pointers: now using auto_ptr instead of SmartPtrIRC
// -- added myCmp, myClone, operator < and !=
//
// Revision 1.2  2009/03/16 07:27:20  bigatti
// -- added necessary "const"
// -- added WDeg1CompTmp (temporary, for testing)
//
// Revision 1.1  2009/02/20 13:27:15  bigatti
// -- renamed from "sugar.[CH]"
//
// Revision 1.3  2009/02/20 11:01:15  bigatti
// -- added NewHomogWSugar for graded=homoeneous case (just constant)
//
// Revision 1.2  2009/02/20 09:51:54  bigatti
// -- changed name: sweetener --> SugarDegree
//
// Revision 1.1  2009/02/09 13:57:05  bigatti
// -- first import
//

#endif
