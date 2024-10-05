#ifndef CoCoA_HomomorphismOps_H
#define CoCoA_HomomorphismOps_H

//   Copyright (c)  2017  John Abbott, Anna Bigatti

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

#include "CoCoA/ring.H"


namespace CoCoA
{

  class RingHom;  // fwd decl -- defined in RingHom.H
  class ideal;    // fwd decl -- defined in ideal.H

  ideal ker(const RingHom& phi);
  ideal ker_H(const RingHom& phi);
  bool IsInjective(const RingHom& phi);
  bool IsSurjective(const RingHom& phi);
  bool IsInImage(const RingHom& phi, ConstRefRingElem y);

  RingElem preimage(const RingHom& phi, ConstRefRingElem y); // a preimage, throws exception if none exists
  RingElem preimage0(const RingHom& phi, ConstRefRingElem y); // a preimage, or 0 if none exists

  
} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/HomomorphismOps.H,v 1.3 2022/02/18 14:11:46 abbott Exp $
// $Log: HomomorphismOps.H,v $
// Revision 1.3  2022/02/18 14:11:46  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.2  2021/06/21 11:16:48  abbott
// Summary: Fixed include guard
//
// Revision 1.1  2017/12/18 13:08:38  bigatti
// -- renamed Ops from Fns
//
// Revision 1.4  2017/12/07 15:41:15  bigatti
// -- added ker_H (ker via homogenous computation: to be completed)
//
// Revision 1.3  2017/08/08 13:49:28  abbott
// Summary: Removed cruft
//
// Revision 1.2  2017/07/24 12:06:37  abbott
// Summary: Added preimage0 (and some cruft)
//
// Revision 1.1  2017/07/22 16:12:51  abbott
// Summary: Addin to allow compilation; these files are NOT YET READY TO BE USED
//
//

#endif
