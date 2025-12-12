#ifndef CoCoA_FinGenModule_H
#define CoCoA_FinGenModule_H

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


// Header file for the classes FinGenModule and FinGenModuleBase
// FG = Finitely Generated (i.e. submodules and quotients of free modules).

namespace CoCoA
{

  class FreeModule;   // fwd decl -- defined in FreeModule.H
  class matrix;       // fwd decl -- defined in matrix.H

  class FinGenModuleBase; // forward declaration
  class FinGenModule: public module
  {
  public:
    // Default copy ctor is OK.
    /*implicit*/ FinGenModule(const module& M);
    FinGenModule(const module& M, const ErrorContext& ErrCtx);
    explicit FinGenModule(const FinGenModuleBase* MPtr);   //??? public???
    // Assignment effectively disabled (because CoCoA::module has no assignment)
    // default dtor is fine

    const FinGenModuleBase* myModulePtr() const;
    const FinGenModuleBase* operator->() const;  // allow const member fns to be called
  };



  class FinGenModuleBase: public ModuleBase
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
  const std::vector<ModuleElem>& gens(const FinGenModule& M); // inline
  const std::vector<ModuleElem>& MinGens(const FinGenModule& M, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit()); // inline
  const std::vector<ModuleElem>& TidyGens(const FinGenModule& M, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit()); // inline
  bool IsFinGenModule(const module& M);
  const FinGenModuleBase* FinGenModulePtr(const module& M);
  const FinGenModuleBase* FinGenModulePtr(const module& M, const ErrorContext& ErrCtx);

  matrix GensAsRows(const FinGenModule& M);
  matrix GensAsCols(const FinGenModule& M);


  /////////////////////////////////////////////////////////////////////////////
  // Inline functions


  inline bool IsFinGenModule(const module& M)
  { return FinGenModulePtr(M) != nullptr; }


  inline FinGenModule::FinGenModule(const module& M):
      module(FinGenModulePtr(M,CoCoA_ERROR_CONTEXT))
  {}

  inline FinGenModule::FinGenModule(const module& M, const ErrorContext& ErrCtx):
      module(FinGenModulePtr(M,ErrCtx))
  {}


  inline FinGenModule::FinGenModule(const FinGenModuleBase* MPtr):
      module(MPtr)
  {}


  inline const FinGenModuleBase* FinGenModule::myModulePtr() const
  { return static_cast<const FinGenModuleBase*>(myRawPtr()); }


  inline const FinGenModuleBase* FinGenModule::operator->() const
  { return static_cast<const FinGenModuleBase*>(myRawPtr()); }


  inline long NumCompts(const module& M)
  { return FinGenModulePtr(M,CoCoA_ERROR_CONTEXT)->myNumCompts(); }


  inline long NumCompts(const ModuleElem& v)
  { return NumCompts(owner(v)); }


  inline const FreeModule& AmbientFreeModule(const FinGenModule& M)
  { return M->myAmbientFreeModule(); }


  inline const std::vector<ModuleElem>& gens(const FinGenModule& M)
  { return M->myGens(); }


  inline long NumGens(const FinGenModule& M)
  { return len(M->myGens()); }


  inline const std::vector<ModuleElem>& MinGens(const FinGenModule& M, const CpuTimeLimit& CheckForTimeout)
  { return M->myMinGens(CheckForTimeout); }


  inline const std::vector<ModuleElem>& TidyGens(const FinGenModule& M, const CpuTimeLimit& CheckForTimeout)
  { return M->myTidyGens(CheckForTimeout); }

} // end of namespace CoCoA

#endif
