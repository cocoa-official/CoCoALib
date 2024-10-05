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


//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/ReductionCog.H,v 1.10 2022/02/18 14:11:48 abbott Exp $
// $Log: ReductionCog.H,v $
// Revision 1.10  2022/02/18 14:11:48  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.9  2018/04/09 16:32:54  bigatti
// -- back to original name
//
// Revision 1.1  2018/04/06 15:34:53  bigatti
// -- renaming ReductionCog.H
//
// Revision 1.7  2014/04/30 16:01:54  abbott
// Summary: Replaced size_t by long
// Author: JAA
//
// Revision 1.6  2012/10/24 12:06:43  abbott
// Revised fwd decl of ConstRefRingElem.
//
// Revision 1.5  2012/10/16 09:48:31  abbott
// Replaced  RefRingElem  by  RingElem&  (several times)
//
// Revision 1.4  2008/09/16 15:03:43  bigatti
// -- added LPPForDiv
// -- changed LPP into LPPForOrd
//
// Revision 1.3  2008/04/21 12:52:34  abbott
// Changed size_t into std::size_t.
//
// Revision 1.2  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.8  2007/03/08 18:42:05  cocoa
// Cleaned up whitespace.
//
// Revision 1.7  2007/02/22 17:27:43  bigatti
// -- removed useless include <memory> (auto_ptr)
//
// Revision 1.6  2006/12/06 17:15:26  cocoa
// -- removed #include "config.H"
//
// Revision 1.5  2006/11/27 14:17:41  cocoa
// -- just a comment
//
// Revision 1.4  2006/11/24 17:53:14  cocoa
// -- using SmartPtrIRC for ref counting (was auto_ptr)
//
// Revision 1.3  2006/10/06 15:17:08  cocoa
// -- myReduce now returns void instead of RingElem
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
// Revision 1.3  2006/04/11 14:16:29  cocoa
// -- reorganization of fns between reduce,SparsePolyRing,GPoly
// -- added optional "len" argument to myAssignReset in ReductionCog
//
// Revision 1.2  2006/04/10 13:38:06  cocoa
// -- added myTmpLM field to avoid repeated creations/deletions
//
// Revision 1.1  2006/03/16 14:28:29  cocoa
// -- first import
//

#endif
