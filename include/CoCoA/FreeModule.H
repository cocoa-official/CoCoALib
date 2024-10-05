#ifndef CoCoA_FreeModule_H
#define CoCoA_FreeModule_H

//   Copyright (c)  2005  John Abbott,  Anna M. Bigatti

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

#include "CoCoA/FGModule.H"
#include "CoCoA/degree.H"

namespace CoCoA
{

  class module;               // fwd decl -- defined in module.H
  class ring;                 // fwd decl -- defined in ring.H
  class RingElem;             // fwd decl -- defined in ring.H
  class ConstRefPPMonoidElem; // fwd decl -- defined in PPMonoid.H
  class ModuleOrdering;       // fwd decl -- defined in ModuleOrdering.H
  class ModuleOrderingCtor;   // fwd decl -- defined in ModuleOrdering.H
  class SparsePolyRing;       // fwd decl -- defined in SparsePolyRing.H

  class FreeModuleBase: public FGModuleBase
  {};

  class FreeModule: public FGModule
  {
  public:
    // Default copy ctor is OK.
    /*implicit*/ FreeModule(const module& M);
    FreeModule(const module& M, const ErrorContext& ErrCtx);
    explicit FreeModule(const FreeModuleBase* MPtr);
    // Assignment is disabled because it is disabled for module & FGModule
///???    ~FreeModule();

    const FreeModuleBase* operator->() const;  // allow const member fns to be called
  };


  FreeModule NewFreeModule(const ring& R, long NumCompts);
  bool IsFreeModule(const module& M);
  const std::vector<ModuleElem>& CanonicalBasis(const FreeModule& F);
  const FreeModuleBase* FreeModulePtr(const module& M);
  const FreeModuleBase* FreeModulePtr(const module& M, const ErrorContext& ErrCtx);
//  const FreeModuleBase* FreeModulePtr(const module& M, const char* const FnName);
  RingElem FirstNonZero(const ModuleElem& v); // copies nz element (returning a ref make the compiler moan)
  long FirstNonZeroPosn(const ModuleElem& v);


  //---------------------------------------------------------------------------
  // Inline functions  -- order is important

  inline const FreeModuleBase* FreeModule::operator->() const
  {
    return static_cast<const FreeModuleBase*>(myRawPtr());
  }


  inline bool IsFreeModule(const module& M)
  {
    return FreeModulePtr(M) != nullptr;
  }


  inline FreeModule::FreeModule(const module& M):
      FGModule(FreeModulePtr(M,CoCoA_ERROR_CONTEXT))
  {}


  inline FreeModule::FreeModule(const module& M, const ErrorContext& ErrCtx):
      FGModule(FreeModulePtr(M,ErrCtx))
  {}


  inline FreeModule::FreeModule(const FreeModuleBase* MPtr):
      FGModule(MPtr)
  {}

//   class GradedFreeModule: public FreeModule  // copied from SparsePolyRing
//   {
//   public:
//     explicit GradedFreeModule(const GradedFreeModuleImpl* MPtr);
//     // Default copy ctor works fine.
//     // Default dtor works fine.
//   private: // disable assignment
//     GradedFreeModule& operator=(const GradedFreeModule& rhs) = delete;
//   public:
//     const GradedFreeModuleImpl* operator->() const; // allow member fns to be called
//   };


  const std::vector<degree>& shifts(const FreeModule& M); ///< defined only if GradedFreeModuleImpl
  const ModuleOrdering& ordering(const FreeModule& M); ///< defined only if GradedFreeModuleImpl

  FreeModule NewFreeModule(const ring& P, long NumCompts, const ModuleOrderingCtor& O); // inline
  FreeModule NewFreeModule(const ring& P, const std::vector<degree>& shifts, const ModuleOrderingCtor& O);
  FreeModule NewFreeModule(const ring& P, const std::vector<degree>& shifts);
  //  bool IsGradedFreeModule(const module& M);

  long LPosn(const ModuleElem& v);
  RingElemAlias LC(const ModuleElem& v);
  ConstRefPPMonoidElem LPP(const ModuleElem& v);
  degree wdeg(const ModuleElem& v);
  long StdDeg(const ModuleElem& v);  
  long deg(const ModuleElem& v);  
  int CmpWDeg(const ModuleElem& v1, const ModuleElem& v2);
  bool IsHomog(const ModuleElem& v);
  ModuleElem homog(const ModuleElem& v, ConstRefRingElem h);

  FreeModule NewFreeModuleForSyz(const std::vector<RingElem>& L, const ErrorContext& ErrCtx);
  FreeModule NewFreeModuleForSyz(const std::vector<RingElem>& L);
  FreeModule NewFreeModuleForSyz(const std::vector<ModuleElem>& L, const ErrorContext& ErrCtx);
  FreeModule NewFreeModuleForSyz(const std::vector<ModuleElem>& L);

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/FreeModule.H,v 1.22 2024/03/26 10:51:13 bigatti Exp $
// $Log: FreeModule.H,v $
// Revision 1.22  2024/03/26 10:51:13  bigatti
// Summary: added NewFreeModuleForSyz with ErrorContext
//
// Revision 1.21  2024/03/15 19:04:45  abbott
// Summary: Arg checking fns now expect a ErrorContext object (instead of string; related to redmine 1726)
//
// Revision 1.20  2023/03/09 22:32:15  abbott
// Summary: Added comment about FirstNonZero
//
// Revision 1.19  2022/02/18 14:11:45  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.18  2022/02/05 19:57:32  abbott
// Summary: Added CanonicalBasis (redmine 1235)
//
// Revision 1.17  2021/10/30 17:12:00  abbott
// Summary: Used keyword delete (redmine 1627)
//
// Revision 1.16  2021/01/07 14:14:21  abbott
// Summary: Corrected copyright
//
// Revision 1.15  2019/03/18 17:13:24  abbott
// Summary: Used nullptr instead of 0
//
// Revision 1.14  2014/07/09 14:28:46  abbott
// Summary: Removed AsFreeModule and AsFGModule
// Author: JAA
//
// Revision 1.13  2013/08/02 14:40:54  bigatti
// -- added LC
// -- changed LPos --> LPosn
//
// Revision 1.12  2013/05/27 16:36:53  bigatti
// -- new constructor for FreeModule (with new implementation of orderings)
//
// Revision 1.11  2013/02/21 17:07:02  bigatti
// -- added NewFreeModuleForSyz
//
// Revision 1.10  2013/02/13 09:08:06  bigatti
// -- added homog(ModuleElem)
//
// Revision 1.9  2013/02/12 16:22:46  bigatti
// -- added FirstNonZero, FirstNonZeroPos
// -- added deg, StdDeg
//
// Revision 1.8  2013/01/23 14:02:31  bigatti
// -- removed doxygen include
//
// Revision 1.7  2011/11/09 13:38:31  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.6  2011/03/10 16:39:35  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.5  2009/12/03 17:26:35  abbott
// Renamed EFGModule to FGModule.
// Renamed ModuleBase member fns  myInit -> myNew, myKill -> myDelete.
// Removed some cruft (old code that was not used by anyone).
//
// Revision 1.4  2008/04/21 12:32:54  abbott
// Corrected size_t into std::size_t in several header files; in some cases,
// replaced size_t with MachineInt (with consequent changes to impl files).
//
// Revision 1.3  2008/03/12 16:35:18  bigatti
// -- changed: IsHomogeneous --> IsHomog
// -- changed: ERR:ZeroPoly --> ERR::ZeroRingElem
//
// Revision 1.2  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.6  2006/12/06 17:11:16  cocoa
// -- removed #include "config.H"
//
// Revision 1.5  2006/11/24 17:17:58  cocoa
// -- reorganized includes of header files
//
// Revision 1.4  2006/11/14 17:46:20  cocoa
// -- changed: reference counting in modules now uses SmartPtrIRC
//
// Revision 1.3  2006/10/06 09:53:06  cocoa
// Moved definition of FreeModuleImpl out of the header file.
// Merged GradedFreeModule code into here so it can access FreeModuleImpl
// (which is no longer visible in the header file).
//
// Revision 1.2  2006/07/18 10:51:25  cocoa
// -- added: FirstNonZeroPos
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.10  2006/05/12 16:10:59  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.9  2006/05/04 14:25:16  cocoa
// -- major cleaning of FreeModule: created GradedFreeModule and moved
//    some code around
//
// Revision 1.8  2006/04/21 14:56:33  cocoa
// Changed return type of myCompt member function: now it returns a
// ConstRefRingElem instead of a RingElem (i.e. a copy).
//
// Revision 1.7  2006/03/27 12:21:26  cocoa
// Minor silly changes to reduce number of complaints from some compiler or other.
//
// Revision 1.6  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.5  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
// Revision 1.4  2005/12/31 12:22:18  cocoa
// Several minor tweaks to silence the Microsoft compiler:
//  - added some missing #includes and using directives
//  - moved some function defns into the right namespace
//  - etc.
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
// Revision 1.3  2005/09/28 11:50:35  cocoa
// -- new code for graded modules
//
// Revision 1.2  2005/06/22 14:42:16  cocoa
// Renamed MemPool data member to myMemMgr
// (seems more sensible than myMemory).
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
// Revision 1.4  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.3  2004/11/08 14:03:32  cocoa
// -- small changes for doxygen support
//
// Revision 1.2  2004/01/28 15:56:49  cocoa
// "Old style" code, brought in alignment with the new coding conventions.
//
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//

#endif
