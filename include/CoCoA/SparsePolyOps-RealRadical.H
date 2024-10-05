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



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/SparsePolyOps-RealRadical.H,v 1.2 2022/02/18 14:11:49 abbott Exp $
// $Log: SparsePolyOps-RealRadical.H,v $
// Revision 1.2  2022/02/18 14:11:49  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.1  2017/12/18 13:10:34  abbott
// Summary: Renamed files (underscore becomes minus)
//
// Revision 1.1  2017/12/01 21:35:41  abbott
// Summary: (approx) RealRadical
//
//

#endif
