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

#include "CoCoA/FinGenModule.H"
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

  class FreeModuleBase: public FinGenModuleBase
  {};

  class FreeModule: public FinGenModule
  {
  public:
    // Default copy ctor is OK.
    /*implicit*/ FreeModule(const module& M);
    FreeModule(const module& M, const ErrorContext& ErrCtx);
    explicit FreeModule(const FreeModuleBase* MPtr);
    // Assignment is disabled because it is disabled for module & FinGenModule
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
      FinGenModule(FreeModulePtr(M,CoCoA_ERROR_CONTEXT))
  {}


  inline FreeModule::FreeModule(const module& M, const ErrorContext& ErrCtx):
      FinGenModule(FreeModulePtr(M,ErrCtx))
  {}


  inline FreeModule::FreeModule(const FreeModuleBase* MPtr):
      FinGenModule(MPtr)
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

#endif
