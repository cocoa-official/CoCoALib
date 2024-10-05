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


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/FGModule.C,v 1.11 2024/08/02 16:00:31 bigatti Exp $
// $Log: FGModule.C,v $
// Revision 1.11  2024/08/02 16:00:31  bigatti
// Summary: just spacing
//
// Revision 1.10  2024/07/22 15:43:57  abbott
// Summary: Renamed macros CoCoA_THROW_ERROR_WITH_CONTEXT so that last char is num of args
//
// Revision 1.9  2024/07/03 14:46:46  bigatti
// Summary: error codes: Not... changed into ReqFGModule, ReqFreeModule, ReqModuleSpPR
//
// Revision 1.8  2024/03/15 18:57:32  abbott
// Summary: Improved use of PolyRing inside fns; employed CoCoA_ERROR_CONTEXT
//
// Revision 1.7  2022/02/18 14:11:54  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.6  2021/01/07 15:07:02  abbott
// Summary: Corrected copyright
//
// Revision 1.5  2020/06/17 15:49:23  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.4  2019/03/19 11:07:07  abbott
// Summary: Replaced 0 by nullptr where appropriate
//
// Revision 1.3  2014/07/09 14:27:53  abbott
// Summary: Removed AsFreeModule and AsFGModule
// Author: JAA
//
// Revision 1.2  2013/07/31 09:50:43  bigatti
// -- added IsZero(module)
//
// Revision 1.1  2009/12/03 17:26:34  abbott
// Renamed EFGModule to FGModule.
// Renamed ModuleBase member fns  myInit -> myNew, myKill -> myDelete.
// Removed some cruft (old code that was not used by anyone).
//
// Revision 1.5  2008/06/30 17:14:32  abbott
// Commented out unused parameter names (to avoid annoying compiler warnings)
//
// Revision 1.4  2008/05/30 12:49:18  abbott
// Added SERIOUS errors to the two stopgap myCompt default implementations
// (which are just a nasty hack to let things compile)
//
// Revision 1.3  2008/05/29 15:42:33  bigatti
// -- added ugly fix for myCompt (should be pure virtual)
//
// Revision 1.2  2007/10/30 17:14:08  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.2  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.2  2006/03/27 12:21:25  cocoa
// Minor silly changes to reduce number of complaints from some compiler or other.
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.2  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.2  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.1  2004/01/28 15:54:09  cocoa
// Sundry additions.
//
//
