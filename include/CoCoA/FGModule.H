#ifndef CoCoA_FGModule_H
#define CoCoA_FGModule_H

//   Copyright (c)  2004-2009  John Abbott,  Anna M. Bigatti

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

#include "CoCoA/CpuTimeLimit.H"
#include "CoCoA/error.H"
#include "CoCoA/module.H"
#include "CoCoA/utils.H" // for len

#include <vector>
//using std::vector;


// Header file for the classes FGModule and FGModuleBase
// FG = Finitely Generated (i.e. submodules and quotients of free modules).

namespace CoCoA
{

  class FreeModule;   // fwd decl -- defined in FreeModule.H
  class matrix;       // fwd decl -- defined in matrix.H

  class FGModuleBase; // forward declaration
  class FGModule: public module
  {
  public:
    // Default copy ctor is OK.
    /*implicit*/ FGModule(const module& M);
    FGModule(const module& M, const ErrorContext& ErrCtx);
    explicit FGModule(const FGModuleBase* MPtr);   //??? public???
    // Assignment effectively disabled (because CoCoA::module has no assignment)
    // default dtor is fine

    const FGModuleBase* myModulePtr() const;
    const FGModuleBase* operator->() const;  // allow const member fns to be called
  };



  class FGModuleBase: public ModuleBase
  {
  public:
    typedef ModuleRawPtr RawPtr;
    typedef const ModuleRawPtr& ConstRawPtr;

    virtual long myNumCompts() const = 0;
    bool IamZero() const override;
    virtual ConstRefRingElem myCompt(const RawPtr& v, long pos) const = 0;
    virtual const FreeModule& myAmbientFreeModule() const = 0;
    virtual const std::vector<ModuleElem>& myGens() const = 0;
    virtual const std::vector<ModuleElem>& myMinGens(const CpuTimeLimit& CheckForTimeout) const = 0;
    virtual const std::vector<ModuleElem>& myTidyGens(const CpuTimeLimit& CheckForTimeout) const = 0;
  };


  long NumCompts(const module& M); // inline
  long NumCompts(const ModuleElem& v); // inline
  const std::vector<ModuleElem>& gens(const FGModule& M); // inline
  const std::vector<ModuleElem>& MinGens(const FGModule& M, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit()); // inline
  const std::vector<ModuleElem>& TidyGens(const FGModule& M, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit()); // inline
  bool IsFGModule(const module& M);
  const FGModuleBase* FGModulePtr(const module& M);
  const FGModuleBase* FGModulePtr(const module& M, const ErrorContext& ErrCtx);

  matrix GensAsRows(const FGModule& M);
  matrix GensAsCols(const FGModule& M);


  /////////////////////////////////////////////////////////////////////////////
  // Inline functions


  inline bool IsFGModule(const module& M)
  { return FGModulePtr(M) != nullptr; }


  inline FGModule::FGModule(const module& M):
      module(FGModulePtr(M,CoCoA_ERROR_CONTEXT))
  {}

  inline FGModule::FGModule(const module& M, const ErrorContext& ErrCtx):
      module(FGModulePtr(M,ErrCtx))
  {}


  inline FGModule::FGModule(const FGModuleBase* MPtr):
      module(MPtr)
  {}


  inline const FGModuleBase* FGModule::myModulePtr() const
  { return static_cast<const FGModuleBase*>(myRawPtr()); }


  inline const FGModuleBase* FGModule::operator->() const
  { return static_cast<const FGModuleBase*>(myRawPtr()); }


  inline long NumCompts(const module& M)
  { return FGModulePtr(M,CoCoA_ERROR_CONTEXT)->myNumCompts(); }


  inline long NumCompts(const ModuleElem& v)
  { return NumCompts(owner(v)); }


  inline const FreeModule& AmbientFreeModule(const FGModule& M)
  { return M->myAmbientFreeModule(); }


  inline const std::vector<ModuleElem>& gens(const FGModule& M)
  { return M->myGens(); }


  inline long NumGens(const FGModule& M)
  { return len(M->myGens()); }


  inline const std::vector<ModuleElem>& MinGens(const FGModule& M, const CpuTimeLimit& CheckForTimeout)
  { return M->myMinGens(CheckForTimeout); }


  inline const std::vector<ModuleElem>& TidyGens(const FGModule& M, const CpuTimeLimit& CheckForTimeout)
  { return M->myTidyGens(CheckForTimeout); }

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/FGModule.H,v 1.18 2024/09/04 15:12:51 bigatti Exp $
// $Log: FGModule.H,v $
// Revision 1.18  2024/09/04 15:12:51  bigatti
// Summary: just removed a comment
//
// Revision 1.17  2024/03/15 19:04:45  abbott
// Summary: Arg checking fns now expect a ErrorContext object (instead of string; related to redmine 1726)
//
// Revision 1.16  2023/07/04 09:15:35  abbott
// Summary: Changed TimeOut to Timeout
//
// Revision 1.15  2022/02/18 14:11:45  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.14  2021/10/30 17:12:53  abbott
// Summary: Used keyword override (redmine 1625)
//
// Revision 1.13  2021/01/07 14:12:29  abbott
// Summary: Corrected copyright
//
// Revision 1.12  2019/03/18 17:13:24  abbott
// Summary: Used nullptr instead of 0
//
// Revision 1.11  2019/03/15 16:35:52  bigatti
// -- added NumGens
//
// Revision 1.10  2018/05/25 09:24:46  abbott
// Summary: Major redesign of CpuTimeLimit (many consequences)
//
// Revision 1.9  2014/07/09 14:28:46  abbott
// Summary: Removed AsFreeModule and AsFGModule
// Author: JAA
//
// Revision 1.8  2014/03/26 15:23:44  bigatti
// -- added MinGens for submodules
//
// Revision 1.7  2013/07/31 14:46:31  bigatti
// -- removed empty lines
//
// Revision 1.6  2013/07/31 09:49:07  bigatti
// -- added IsZero(module)
// -- deleted empty lines in inline functions
//
// Revision 1.5  2013/02/21 16:55:18  bigatti
// -- added GensAsRows, GensAsCols
//
// Revision 1.4  2012/10/12 12:38:18  abbott
// Removed element accessor (via operator[]) and non-const mem fn  ModuleBase::myCompt.
//
// Revision 1.3  2011/11/09 13:38:31  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.2  2011/03/10 16:39:35  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.1  2009/12/03 17:26:35  abbott
// Renamed EFGModule to FGModule.
// Renamed ModuleBase member fns  myInit -> myNew, myKill -> myDelete.
// Removed some cruft (old code that was not used by anyone).
//
// Revision 1.6  2008/05/30 12:47:48  abbott
// Removed a blank line.
//
// Revision 1.5  2008/05/29 15:42:34  bigatti
// -- added ugly fix for myCompt (should be pure virtual)
//
// Revision 1.4  2008/04/21 12:32:54  abbott
// Corrected size_t into std::size_t in several header files; in some cases,
// replaced size_t with MachineInt (with consequent changes to impl files).
//
// Revision 1.3  2008/02/11 15:11:00  abbott
// Added new NumCompts function (for ModuleElems).
//
// Revision 1.2  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.4  2006/12/06 17:11:16  cocoa
// -- removed #include "config.H"
//
// Revision 1.3  2006/11/14 17:46:20  cocoa
// -- changed: reference counting in modules now uses SmartPtrIRC
//
// Revision 1.2  2006/10/06 14:04:16  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.6  2006/04/21 14:56:33  cocoa
// Changed return type of myCompt member function: now it returns a
// ConstRefRingElem instead of a RingElem (i.e. a copy).
//
// Revision 1.5  2006/03/27 12:21:26  cocoa
// Minor silly changes to reduce number of complaints from some compiler or other.
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
// Revision 1.2  2005/09/28 11:50:35  cocoa
// -- new code for graded modules
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
// Revision 1.3  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.2  2004/06/29 17:10:22  cocoa
// Partially tidied use of "protected" and "private" in various
// base classes.  Checking in at the end of the day -- it works,
// and I wouldn't want it to be lost next time point's disk
// misbehaves.
//
// Revision 1.1  2004/01/28 15:54:09  cocoa
// Sundry additions.
//
//

#endif
