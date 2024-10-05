#ifndef CoCoA_RootBound_H
#define CoCoA_RootBound_H

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

#include "CoCoA/ring.H"

namespace CoCoA
{

  class BigRat; // fwd decl -- defined in BigRat.H

  BigRat RootBound(ConstRefRingElem f, long NumIters = -1);
  BigRat RootBound2(ConstRefRingElem f, long NumIters = -1);

  BigRat RootBound_Cauchy(ConstRefRingElem f);
  BigRat RootBound_Lagrange(ConstRefRingElem f);
  BigRat RootBound_Birkhoff(ConstRefRingElem f);
  BigRat RootBound_LMS(ConstRefRingElem f);
  BigRat RootBound_simple(ConstRefRingElem f);

  RingElem RootBoundTransform(ConstRefRingElem f);
//  RingElem RootBound_preprocess(ConstRefRingElem f);  // Should be public???

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/RootBound.H,v 1.4 2022/02/18 14:11:48 abbott Exp $
// $Log: RootBound.H,v $
// Revision 1.4  2022/02/18 14:11:48  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.3  2020/07/28 08:01:48  abbott
// Summary: Added RootBoundTransform
//
// Revision 1.2  2017/12/12 14:17:24  abbott
// Summary: Added new fns RootBound_Cauchy and RootBound_LMS; commented out preprocess
//
// Revision 1.1  2017/09/14 15:54:54  abbott
// Summary: Added RootBound
//
//

#endif
