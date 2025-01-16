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

#endif
