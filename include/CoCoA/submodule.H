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


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/submodule.H,v 1.21 2022/03/07 11:01:21 abbott Exp $
// $Log: submodule.H,v $
// Revision 1.21  2022/03/07 11:01:21  abbott
// Summary: Added new fn syntax syz(L)
//
// Revision 1.20  2022/03/07 10:19:19  bigatti
// Summary: added fwd decl RingElem
//
// Revision 1.19  2022/03/07 09:47:10  bigatti
// Summary: added func syz(const FreeModule& F, const std::vector<RingElem>& g);
//
// Revision 1.18  2022/02/18 14:11:52  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.17  2021/01/07 14:48:33  abbott
// Summary: Corrected copyright
//
// Revision 1.16  2017/11/20 20:10:26  bigatti
// -- mimimalized --> SubmoduleOfMinGens
//
// Revision 1.15  2014/04/17 13:37:49  bigatti
// -- MatrixViews --> MatrixView
//
// Revision 1.14  2014/04/10 13:04:27  bigatti
// -- minimalized(FGModule)
//
// Revision 1.13  2014/04/09 13:11:06  bigatti
// -- added submodule(vector<ModuleElem>)
//
// Revision 1.12  2013/07/31 14:46:53  bigatti
// -- added LT(module)
//
// Revision 1.11  2013/07/31 09:49:36  bigatti
// -- added IsHomog
//
// Revision 1.10  2013/06/10 10:15:18  bigatti
// -- renamed args for IsContained
//
// Revision 1.9  2013/06/06 05:47:34  bigatti
// -- pseudoctor is now called "submodule" instead of "NewSubmodule"
// -- added pseudoctor with 1,2,3,4 generators
//
// Revision 1.8  2013/06/03 10:05:59  bigatti
// -- added IsElem, IsContained for modules
//
// Revision 1.7  2013/03/25 17:29:59  abbott
// Changed formal parameter names (M for matrix, N for submodule).
//
// Revision 1.6  2013/03/22 10:26:48  bigatti
// -- added SyzOfGens for modules
//
// Revision 1.5  2013/02/21 17:32:54  bigatti
// -- added SyzOfGens
//
// Revision 1.4  2013/01/30 15:45:19  bigatti
// -- added NewSubmoduleCols/Rows
//
// Revision 1.3  2009/12/03 17:26:34  abbott
// Renamed EFGModule to FGModule.
// Renamed ModuleBase member fns  myInit -> myNew, myKill -> myDelete.
// Removed some cruft (old code that was not used by anyone).
//
// Revision 1.2  2007/10/30 17:14:11  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.3  2006/11/24 17:41:51  cocoa
// -- reorganized includes of header files
//
// Revision 1.2  2006/10/06 10:15:53  cocoa
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
// Revision 1.4  2006/03/12 21:28:34  cocoa
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
// Revision 1.5  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.4  2004/11/11 13:05:49  cocoa
// -- added \include *.txt for doxygen
//
// Revision 1.3  2004/11/09 15:57:01  cocoa
// -- minor changes for doxygen
//
// Revision 1.2  2004/01/28 15:37:58  cocoa
// Fairly major update: resuscitated "old style" code which didn't compile
// under the current organization.
//
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//

#endif
