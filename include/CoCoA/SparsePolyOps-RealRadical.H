#ifndef CoCoA_SparsePolyOps_RealRadical_H
#define CoCoA_SparsePolyOps_RealRadical_H

//   Copyright (c)  2017  John Abbott,  Anna M. Bigatti

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

  class bool3;     // fwd decl -- defined in bool3.H
  class RingElem;  // fwd decl -- defined in ring.H
  
  bool3 HasRealRoot3(const RingElem& f);
  RingElem RealRadical(const RingElem& f);


} // end of namespace CoCoA

#endif
