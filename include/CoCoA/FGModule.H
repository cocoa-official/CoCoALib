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

#endif
