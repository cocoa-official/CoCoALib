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
