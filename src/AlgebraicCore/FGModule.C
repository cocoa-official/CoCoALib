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


// Source code for classes FGModule and FGModuleBase

#include "CoCoA/FGModule.H"
#include "CoCoA/utils.H" // for len


namespace CoCoA
{

  const FGModuleBase* FGModulePtr(const module& M)
  { return dynamic_cast<const FGModuleBase*>(M.myRawPtr()); }

  
  const FGModuleBase* FGModulePtr(const module& M, const ErrorContext& ErrCtx)
  {
    const FGModuleBase* ptr = FGModulePtr(M);
    if (ptr == nullptr)
      CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::ReqFGModule, ErrCtx);
    return ptr;
  }


  bool FGModuleBase::IamZero() const
  {
    const std::vector<ModuleElem>& g = myGens();
    for (long i=0; i<len(g); ++i)
      if (!IsZero(g[i])) return false;
    return true;
  }


} // end of namespace CoCoA
