//   Copyright (c)  2005  John Abbott and Anna M. Bigatti

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


// Implementation file for the class submodule

#include "CoCoA/QuotientModule.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/FGModule.H"
#include "CoCoA/FreeModule.H"
#include "CoCoA/MemPool.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/ring.H"

#include <vector>
using std::vector;
#include <iostream>
using std::ostream;


namespace CoCoA
{

  class QuotientModuleImpl: public FGModuleBase
  {
    // Two typedefs to save typing.
    typedef ModuleBase::RawPtr RawPtr;
    typedef const ModuleBase::RawPtr& ConstRawPtr;

  public:
    QuotientModuleImpl(const FGModule& Mnumer, const FGModule& Mdenom);
    long myNumCompts() const override;
    const ring& myRing() const override;
    const FreeModule& myAmbientFreeModule() const override;
    const std::vector<ModuleElem>& myGens() const override;
    const std::vector<ModuleElem>& myMinGens(const CpuTimeLimit& CheckForTimeout) const override;
    const std::vector<ModuleElem>& myTidyGens(const CpuTimeLimit& CheckForTimeout) const override;

    const ModuleElem& myZero() const override;
    void myNew(RawPtr& rawv) const override;
    void myNew(RawPtr& rawv, ConstRawPtr rawt) const override;
    void myDelete(RawPtr& rawv) const override;                                       // destroys v (incl all resources)
    void mySwap(RawPtr& rawv, RawPtr& raww) const override;                           // swap(v, w);
    void myAssign(RawPtr& rawlhs, ConstRawPtr rawv) const override;                   // lhs = v;
    ConstRefRingElem myCompt(const RawPtr& rawv, long pos) const override;            ///< v[pos]
    void myNegate(RawPtr& rawlhs, ConstRawPtr rawv) const override;                   // lhs = -v
    void myAdd(RawPtr& rawlhs, ConstRawPtr rawv, ConstRawPtr raww) const override;    // lhs = v+w;
    void mySub(RawPtr& rawlhs, ConstRawPtr rawv, ConstRawPtr raww) const override;    // lhs = v-w;

    void myMul(RawPtr& rawlhs, RingElemConstRawPtr rawx, ConstRawPtr rawv) const override; // lhs = r*v;
    void myDiv(RawPtr& rawlhs, RingElemConstRawPtr rawx, ConstRawPtr rawv) const override; // lhs = (1/r)*v;
    void myOutput(std::ostream& out, ConstRawPtr rawv) const override;                // out << v
    void myOutputSelf(std::ostream& out) const override;                              // out << M
    void myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawv) const override;            // OMOut << v
    void myOutputSelf_OM(OpenMathOutput& OMOut) const override;                          // OMOut << M
    bool myIsZero(ConstRawPtr rawv) const override;                                   // v == 0
//???    bool IsZeroAddMul(RawPtr& rawlhs, RingElemConstRawPtr rawy, ConstRawPtr rawz) const;  // lhs += y*z, result says whether lhs == 0.
    bool myIsEqual(ConstRawPtr, ConstRawPtr) const override;

  private: // data members
    const FreeModule myM;
    const FGModule myNumer;
    const FGModule myDenom;
    std::vector<ModuleElem> myGensArray;
    mutable bool myTidyGensIsValid;
    mutable std::vector<ModuleElem> myTidyGensArray;
//???    std::vector<ModuleElem>& ComputeTidyGens() const;
  };



  QuotientModuleImpl::QuotientModuleImpl(const FGModule& Mnumer, const FGModule& Mdenom):
      myM(AmbientFreeModule(Mnumer)),
      myNumer(Mnumer),
      myDenom(Mdenom)
  {
    if (AmbientFreeModule(Mnumer) != AmbientFreeModule(Mdenom))
      CoCoA_THROW_ERROR("Modules reside in different ambient free modules", "QuotientModule ctor");
    //??????? MISSING CODE  inherit gens of Mnumer directly??
    myRefCountZero();
  }


  long QuotientModuleImpl::myNumCompts() const
  {
    return NumCompts(myM);
  }


  const ring& QuotientModuleImpl::myRing() const
  {
    return RingOf(myM);
  }


  const FreeModule& QuotientModuleImpl::myAmbientFreeModule() const
  {
    return myM;
  }


  const vector<ModuleElem>& QuotientModuleImpl::myGens() const
  {
    return myGensArray;
  }


  const vector<ModuleElem>& QuotientModuleImpl::myMinGens(const CpuTimeLimit& /*CheckForTimeout*/) const
  {
//    if (!IsSparsePolyRing(myRing()))
//      CoCoA_THROW_ERROR(ERR::NotSparsePolyRing, "QuotientModuleImpl::myMinGens()");
    CoCoA_THROW_ERROR(ERR::NYI, "QuotientModuleImpl::myMinGens");
    return myGensArray;
  }


  const vector<ModuleElem>& QuotientModuleImpl::myTidyGens(const CpuTimeLimit& /*CheckForTimeout*/) const
  {
    if (!myTidyGensIsValid)
    {
      CoCoA_THROW_ERROR(ERR::NYI, "QuotientModuleImpl::myTidyGens");
    }
    return myGensArray;
  }


  const ModuleElem& QuotientModuleImpl::myZero() const
  {
    return zero(myM);
  }


  void QuotientModuleImpl::myNew(RawPtr& rawv) const
  {
    myM->myNew(rawv);
  }


  void QuotientModuleImpl::myNew(RawPtr& rawv, ConstRawPtr rawcopy) const
  {
    myM->myNew(rawv, rawcopy);
  }


  void QuotientModuleImpl::myDelete(RawPtr& rawv) const
  {
    myM->myDelete(rawv);
  }


  void QuotientModuleImpl::mySwap(RawPtr& rawv, RawPtr& raww) const
  {
    myM->mySwap(rawv, raww);
  }


  void QuotientModuleImpl::myAssign(RawPtr& rawlhs, ConstRawPtr rawv) const
  {
    myM->myAssign(rawlhs, rawv);
  }


  ConstRefRingElem QuotientModuleImpl::myCompt(const RawPtr& rawv, long pos) const
  {
    CoCoA_ASSERT(0 <= pos && pos < myNumCompts());
    return myM->myCompt(rawv, pos);
  }

  void QuotientModuleImpl::myNegate(RawPtr& rawlhs, ConstRawPtr rawv) const
  {
    myM->myNegate(rawlhs, rawv);
  }


  void QuotientModuleImpl::myAdd(RawPtr& rawlhs, ConstRawPtr rawv, ConstRawPtr raww) const
  {
    myM->myAdd(rawlhs, rawv, raww);
    //??? reduce module Mdenom
  }


  void QuotientModuleImpl::mySub(RawPtr& rawlhs, ConstRawPtr rawv, ConstRawPtr raww) const
  {
    myM->mySub(rawlhs, rawv, raww);
    //??? reduce module Mdenom
  }


  void QuotientModuleImpl::myMul(RawPtr& rawlhs, RingElemConstRawPtr rawx, ConstRawPtr rawv) const
  {
    myM->myMul(rawlhs, rawx, rawv);
    //??? reduce module Mdenom
  }


  void QuotientModuleImpl::myDiv(RawPtr& rawlhs, RingElemConstRawPtr rawx, ConstRawPtr rawv) const
  {
    myM->myDiv(rawlhs, rawx, rawv);
    //??? reduce module Mdenom
  }


  void QuotientModuleImpl::myOutput(ostream& out, ConstRawPtr rawv) const
  {
    if (!out) return;  // short-cut for bad ostreams
    myM->myOutput(out, rawv);
  }


  void QuotientModuleImpl::myOutputSelf(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "QuotientModule(" << myNumer << ", " << myDenom << ")";
  }


  void QuotientModuleImpl::myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawv) const
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("???", "ModuleElement"); // BUG: what should this OMSymbol be???
    OMOut << myM;
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("???", "list"); // BUG: what should this OMSymbol be???
    OMOut << myNumCompts();
    myM->myOutput_OM(OMOut, rawv); // BUG: this should be a "naked" output???
    OMOut->mySendApplyEnd();
    OMOut->mySendApplyEnd();
  }


  void QuotientModuleImpl::myOutputSelf_OM(OpenMathOutput& OMOut) const
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("???", "QuotientModule"); // BUG: what should this OMSymbol be???
    OMOut << myNumer << myDenom;
    OMOut->mySendApplyEnd();
  }


  bool QuotientModuleImpl::myIsZero(ConstRawPtr rawv) const
  {
    return myM->myIsZero(rawv);
  }


//???    bool IsZeroAddMul(RawPtr& rawlhs, RingElemConstRawPtr rawy, ConstRawPtr rawz) const;  // lhs += y*z, result says whether lhs == 0.


  bool QuotientModuleImpl::myIsEqual(ConstRawPtr rawv, ConstRawPtr raww) const
  {
    return myM->myIsEqual(rawv, raww);
  }


  FGModule NewQuotientModule(const FGModule& Mnumer, const FGModule& Mdenom)
  {
    return FGModule(new QuotientModuleImpl(Mnumer, Mdenom));
  }


}  // end of namespace CoCoA

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/QuotientModule.C,v 1.21 2023/07/04 09:20:12 abbott Exp $
// $Log: QuotientModule.C,v $
// Revision 1.21  2023/07/04 09:20:12  abbott
// Summary: Changed TimeOut to Timeout
//
// Revision 1.20  2023/06/27 17:44:28  abbott
// Summary: Added commented out arg check
//
// Revision 1.19  2023/03/27 14:00:25  abbott
// Summary: Commented out unused CheckForTimeout arg (in NYU fns)
//
// Revision 1.18  2022/02/18 14:11:56  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.17  2022/02/08 20:18:54  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.16  2021/10/30 16:49:49  abbott
// Summary: Used keywords override & delete (redmine 1625 & 1627)
//
// Revision 1.15  2021/01/07 15:16:52  abbott
// Summary: Corrected copyright
//
// Revision 1.14  2020/06/17 15:49:25  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.13  2020/02/11 16:56:41  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.12  2020/02/11 16:12:19  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.11  2018/05/25 09:24:46  abbott
// Summary: Major redesign of CpuTimeLimit (many consequences)
//
// Revision 1.10  2018/05/18 12:15:04  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.9  2014/07/30 14:07:57  abbott
// Summary: Changed BaseRing into RingOf; myBaseRing --> myRing
// Author: JAA
//
// Revision 1.8  2014/03/26 15:23:44  bigatti
// -- added MinGens for submodules
//
// Revision 1.7  2012/10/12 12:38:18  abbott
// Removed element accessor (via operator[]) and non-const mem fn  ModuleBase::myCompt.
//
// Revision 1.6  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.5  2011/03/10 16:39:34  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.4  2009/12/03 17:26:34  abbott
// Renamed EFGModule to FGModule.
// Renamed ModuleBase member fns  myInit -> myNew, myKill -> myDelete.
// Removed some cruft (old code that was not used by anyone).
//
// Revision 1.3  2008/05/30 12:50:04  abbott
// Aligned some comments.
//
// Revision 1.2  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.6  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.5  2007/01/15 13:39:54  cocoa
// -- added prefix "raw" to RawPtr arguments names
//
// Revision 1.4  2006/11/24 17:06:10  cocoa
// -- reorganized includes of header files
//
// Revision 1.3  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.2  2006/10/06 10:15:52  cocoa
// In response to Susan's bug: a fiasco when compiling with CoCoA_MEMPOOL_DEBUG
// set wrongly.  Moved several implementation classes out of their header files
// into the implementation files.  Several functions had to be uninlined.
// Also corrected position of #include, etc.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.7  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.6  2006/04/21 14:56:33  cocoa
// Changed return type of myCompt member function: now it returns a
// ConstRefRingElem instead of a RingElem (i.e. a copy).
//
// Revision 1.5  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.4  2006/03/12 21:28:33  cocoa
// Major check in after many changes
//
// Revision 1.3  2005/11/29 13:04:47  cocoa
// -- added "const" to myCompt argument
//
// Revision 1.2  2005/11/24 16:09:38  cocoa
// -- added operator[] for ModuleElem
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
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
// Revision 1.3  2004/11/11 14:19:20  cocoa
// -- moved CVS log to the bottom
// -- minor changes for doxygen
//
// Revision 1.2  2004/05/24 15:52:14  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.1  2004/01/28 15:54:09  cocoa
// Sundry additions.
//
//
