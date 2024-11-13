#ifndef CoCoA_submodule_H
#define CoCoA_submodule_H

//   Copyright (c)  2005,2006  John Abbott, and Anna M. Bigatti

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


// Header file for the class submodule


#include <vector>
// using vector

namespace CoCoA
{

  class ConstMatrixView;  // forward declaration -- defined in MatrixView.H
  class FGModule;    // forward declaration -- defined in FGModule.H
  class FreeModule;    // forward declaration -- defined in FGModule.H
  class ModuleElem;  // forward declaration -- defined in module.H
  class RingElem;  // forward declaration -- defined in ring.H
  class ideal;  // forward declaration -- defined in ideal.H
  class module;  // forward declaration -- defined in module.H

  FGModule submodule(const FGModule& M, const std::vector<ModuleElem>& gens);
  FGModule submodule(const std::vector<ModuleElem>& gens);
  FGModule submodule(const ModuleElem& v1);
  FGModule submodule(const ModuleElem& v1, const ModuleElem& v2);
  FGModule submodule(const ModuleElem& v1, const ModuleElem& v2, const ModuleElem& v3);
  FGModule submodule(const ModuleElem& v1, const ModuleElem& v2, const ModuleElem& v3, const ModuleElem& v4);
  FGModule SubmoduleCols(const FGModule& F, ConstMatrixView M);
  FGModule SubmoduleRows(const FGModule& F, ConstMatrixView M);
  FGModule SubmoduleOfMinGens(const FGModule& F);
  FGModule syz(const std::vector<RingElem>& g);
  FGModule syz(const FreeModule& F, const std::vector<RingElem>& g);
  FGModule SyzOfGens(const FreeModule& F, const ideal& I);
  FGModule SyzOfGens(const FreeModule& F, const FGModule& N);
  FGModule SyzOfGens(const ideal& I);
  FGModule SyzOfGens(const FGModule& N);
  bool IsElem(const ModuleElem& v, const module& M);
  bool IsContained(const module& M, const module& N);
  bool IsHomog(const module& M);
  FGModule LT(const module& M);

} // end of namespace CoCoA

#endif
