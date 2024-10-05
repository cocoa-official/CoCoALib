//   Copyright (c)  2007  Anna Bigatti,  John Abbott

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


namespace CoCoA
{

  namespace CoCoAServerOperationsFromCoCoALib
  {
    bool RegisterOpsOnce();
    bool GlobalDummyVar = RegisterOpsOnce();
  }

  namespace CoCoAServerOperationsFromUser
  {
    bool RegisterOpsOnce();
    bool GlobalDummyVar = RegisterOpsOnce();
  }

#ifdef CoCoA_WITH_FROBBY
  namespace CoCoAServerOperationsFromFrobby
  {
    bool RegisterOpsOnce();
    bool GlobalDummyVar = RegisterOpsOnce();
  }
#endif

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/server/RegisterServerOps.H,v 1.4 2022/02/22 20:40:04 abbott Exp $
// $Log: RegisterServerOps.H,v $
// Revision 1.4  2022/02/22 20:40:04  abbott
// Summary: Updated copyright message (redmine 1555)
//
// Revision 1.3  2010/02/03 14:19:00  bigatti
// -- removed RegisterServerOpsFrobby.H (declaration moved into RegisterServerOps.H)
// -- added RegisterServerOpsUser declaration in RegisterServerOps.H
//
// Revision 1.2  2007/10/30 17:14:06  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
//
