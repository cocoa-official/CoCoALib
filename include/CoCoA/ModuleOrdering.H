#ifndef CoCoA_ModuleOrdering_H
#define CoCoA_ModuleOrdering_H

//   Copyright (c)  2005-2007,2011,2013 John Abbott, Anna Bigatti e Massimo Caboara

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


#include "CoCoA/PPOrdering.H"
#include "CoCoA/SmartPtrIRC.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/degree.H"
#include "CoCoA/ring.H"

#include <iosfwd>
//using std::ostream;
#include <vector>
//using std::vector;


namespace CoCoA
{

  class OpenMathOutput; // fwd decl -- defined in OpenMath.H

  class ModuleOrderingBase; // fwd decl -- defined later in this file

  class ModuleOrdering
  {
  public:
    explicit ModuleOrdering(const ModuleOrderingBase* ptr): mySmartPtr(ptr) {}
    // assignment disabled because SmartPtrIRC has no assignment
    const ModuleOrderingBase* operator->() const { return mySmartPtr.operator->(); }  ///< Allow const member fns to be called.
    const ModuleOrderingBase* myRawPtr() const { return mySmartPtr.myRawPtr(); } ///< Used by "downcasting" functions Is(..), As(...) etc.
    bool operator==(const ModuleOrdering& MOrd) const { return mySmartPtr==MOrd.mySmartPtr; }
  private: // data members
    SmartPtrIRC<const ModuleOrderingBase> mySmartPtr;
  };


  std::ostream& operator<<(std::ostream& out, const ModuleOrdering& MOrd);
  OpenMathOutput& operator<<(OpenMathOutput& OMOut, const ModuleOrdering& MOrd);

  //---------------------------------------------------------------------------
  class ModuleOrderingBase: private IntrusiveReferenceCount
  {
    friend class SmartPtrIRC<const ModuleOrderingBase>; // morally "friend ModuleOrdering", so it can alter ref count.

  protected:
    ModuleOrderingBase(const PPOrdering& PPO, const std::vector<degree>& shifts);
    ModuleOrderingBase(const PPOrdering& PPO, const std::vector<degree>& shifts, const std::vector<long>& perm);
    virtual ~ModuleOrderingBase();
  public: // disable copy ctor & assignment
    ModuleOrderingBase(const ModuleOrderingBase& copy) = delete;
    ModuleOrderingBase& operator=(const ModuleOrderingBase&) = delete;
  public:
    virtual void myOutputSelf(std::ostream& out) const;
    virtual void myOutputSelf_OM(OpenMathOutput& OMOut) const;
    virtual void myOutputName(std::ostream& out) const = 0;
    virtual void myOutputName_OM(OpenMathOutput& OMOut) const = 0;
    //    virtual void myOrdMat(matrix& M) const = 0;
    const std::vector<degree>& myShifts() const;
    const std::vector<long>& myPerm() const;
    const PPOrdering& myPPOrdering() const;

  protected:
    PPOrdering myPPO;                       ///< the PPOrdering (on the base ring)
    std::vector<degree> myShiftsValue;      ///< the components shifts (degrees)
    std::vector<long> myPermutationValue;   ///< the components priorities
  };

  const std::vector<degree>& shifts(const ModuleOrdering& O);
  long NumComponents(const ModuleOrdering& MOrd);
  long GradingDim(const ModuleOrdering& MOrd);
  const PPOrdering& ModPPOrdering(const ModuleOrdering& MOrd);

  //---------------------------------------------------------------------------
  // Pseudo-constructors for ModuleOrderings
  ModuleOrdering NewOrdPosn(const PPOrdering& PPO, long NumComponents);
  ModuleOrdering NewPosnOrd(const PPOrdering& PPO, long NumComponents);
  ModuleOrdering NewWDegPosnOrd(const PPOrdering& PPO, long NumComponents);

  ModuleOrdering NewOrdPosn(const PPOrdering& PPO, const std::vector<degree>& shifts);
  ModuleOrdering NewWDegPosnOrd(const PPOrdering& PPO, const std::vector<degree>& shifts);
  ModuleOrdering NewPosnOrd(const PPOrdering& PPO, const std::vector<degree>& shifts);

  ModuleOrdering NewOrdPosn(const PPOrdering& PPO, const std::vector<long>& perm);
  ModuleOrdering NewWDegPosnOrd(const PPOrdering& PPO, const std::vector<long>& perm);
  ModuleOrdering NewOrdPosn(const PPOrdering& PPO, const std::vector<degree>& shifts, const std::vector<long>& perm);

  ModuleOrdering NewWDegPosnOrd(const PPOrdering& PPO, const std::vector<degree>& shifts, const std::vector<long>& perm);

  //  ModuleOrdering NewMatrixOrdering(ConstMatrix OrderMatrix, ??);

  // ??? Compiler bug: the friend declarations below need explicit namespace
  //     scoping with g++ 3.2.2, otherwise g++ thinks New... etc are
  //     in the CoCoA::MOrd namespace.  Very odd!

  // ANNA: philosophical problem: WDeg must include shifts; what about Ord?
  bool IsOrdPosn(const ModuleOrdering& MOrd);
  bool IsPosnOrd(const ModuleOrdering& MOrd);
  bool IsWDegPosnOrd(const ModuleOrdering& MOrd);

  //----------------------------------------------------------------------
  class ModuleOrderingCtor
  {
  public:
    virtual ~ModuleOrderingCtor() {}
    virtual ModuleOrdering operator()(const PPOrdering& PPO, long NumCompts) const = 0;
    virtual ModuleOrdering operator()(const PPOrdering& PPO, const std::vector<degree>& shifts) const = 0;
  };

  struct PosnOrdCtor: public ModuleOrderingCtor
  {
    ModuleOrdering operator()(const PPOrdering& PPO, long NumCompts) const
    { return NewPosnOrd(PPO, NumCompts); }
    ModuleOrdering operator()(const PPOrdering& PPO, const std::vector<degree>& shifts) const; // NYI!!!
  };

  struct OrdPosnCtor: public ModuleOrderingCtor
  {
    ModuleOrdering operator()(const PPOrdering& PPO, long NumCompts) const
    { return NewOrdPosn(PPO, NumCompts); }
    ModuleOrdering operator()(const PPOrdering& PPO, const std::vector<degree>& shifts) const
    { return NewOrdPosn(PPO, shifts); }
  };

  struct WDegPosnOrdCtor: public ModuleOrderingCtor
  {
    ModuleOrdering operator()(const PPOrdering& PPO, long NumCompts) const
    { return NewWDegPosnOrd(PPO, NumCompts); }
    ModuleOrdering operator()(const PPOrdering& PPO, const std::vector<degree>& shifts) const
    { return NewWDegPosnOrd(PPO, shifts); }
  };

  extern OrdPosnCtor OrdPosn;
  extern PosnOrdCtor PosnOrd;
  extern WDegPosnOrdCtor WDegPosnOrd;

  //---------------------------------------------------------------------

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/ModuleOrdering.H,v 1.6 2022/02/18 14:11:46 abbott Exp $
// $Log: ModuleOrdering.H,v $
// Revision 1.6  2022/02/18 14:11:46  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.5  2022/02/08 20:18:53  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.4  2021/10/30 17:07:48  abbott
// Summary: Used keyword delete (redmine 1627)
//
// Revision 1.3  2015/12/01 13:11:00  abbott
// Summary: Changed mem fn PPOrderingCtor::myCtor into operator(); also for ModuleOrderingCtor; see issue 829
//
// Revision 1.2  2014/07/31 13:10:46  bigatti
// -- GetMatrix(PPO) --> OrdMat(PPO)
// -- added OrdMat and GradingMat to PPOrdering, PPMonoid, SparsePolyRing
//
// Revision 1.1  2013/06/03 08:51:58  bigatti
// -- was ModuleTermOrdering
//
// Revision 1.5  2013/05/31 10:31:09  abbott
// Moved NYI impl body of PosnOrdCtor::myCtor to *.C file to avoid multiple compiler
// warnings about unused parameter.
//
// Revision 1.4  2013/05/27 16:35:05  bigatti
// -- major reorganisation of the implementation, changed names
// ---- WDegPosTO --> WDegPosnOrd,  WDegTOPos --> OrdPosn,  PosWDegTO --> PosnOrd
//
// Revision 1.3  2011/03/10 16:39:35  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.2  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.7  2006/12/06 17:14:12  cocoa
// -- removed #include "config.H"
//
// Revision 1.6  2006/11/24 17:14:10  cocoa
// -- reorganized includes of header files
// -- SmartPtrIRC for reference counting
//
// Revision 1.5  2006/11/10 13:30:57  cocoa
// -- fixed: "const &" to PPOrdering arguments
// -- some more documentation
//
// Revision 1.4  2006/11/10 13:06:03  cocoa
// -- cleaned code
//
// Revision 1.3  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
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
// Revision 1.7  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.6  2006/05/11 15:59:23  cocoa
// -- changed reference count is done using SmartPtrIRC
//
// Revision 1.5  2006/05/04 14:19:02  cocoa
// -- moved some code from .H to .C
//
// Revision 1.4  2006/04/28 11:32:17  cocoa
// -- moved concrete class definition from .H to .C
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
// Revision 1.1  2005/09/28 11:50:35  cocoa
// -- new code for graded modules
//

#endif
