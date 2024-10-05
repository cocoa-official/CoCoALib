#ifndef CoCoA_SparsePolyOps_IDEAL_FGLM_H
#define CoCoA_SparsePolyOps_IDEAL_FGLM_H

//   Copyright (c)  2024  John Abbott and Anna M. Bigatti

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


// Header file for functions FGLM etc.
//#include "CoCoA/ring.H"
//#include "CoCoA/VerificationLevel.H"

#include <vector>
// using std::vector;

namespace CoCoA
{

  class RingElem;         // fwd decl -- defined in ring.H
  class RingElemAlias;    // fwd decl -- defined in ring.H
  class ideal;            // fwd decl -- defined in ideal.H
  class PPMonoidElem;     // fwd decl -- defined in PPMonoid.H
  class SparsePolyRing;         // fwd decl -- defined in PolyRing.H
  //  class VerificationLevel;// fwd decl -- defined in VerificationLevel.H

  std::vector<RingElem> FGLM(const SparsePolyRing& P_new, const ideal& I);
  

  void FGLM(std::vector<PPMonoidElem>& QB,
            std::vector<RingElem>& GBasis,
            const SparsePolyRing& P_new,
            const ideal& I);


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/SparsePolyOps-ideal-FGLM.H,v 1.1 2024/05/13 13:12:28 bigatti Exp $
// $Log: SparsePolyOps-ideal-FGLM.H,v $
// Revision 1.1  2024/05/13 13:12:28  bigatti
// Summary: first import: added FGLM
//

#endif
