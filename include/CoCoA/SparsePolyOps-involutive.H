#ifndef CoCoA_SparsePolyOps_involutive_H
#define CoCoA_SparsePolyOps_involutive_H

//   Copyright (c)  2017  John Abbott,  Anna M. Bigatti
//   Author: Mario Albert

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

#include "CoCoA/PolyRing.H"
#include "CoCoA/FGModule.H"

#include <map>
#include <vector>

namespace CoCoA
{

  class ideal;     // fwd decl -- defined in ideal.H
  class RingElem;  // fwd decl -- defined in ring.H
  class PPMonoidElem;  // fwd decl -- defined in PPMonoid.H

  //----- Involutive -------------------------------------------------
  namespace Involutive {
    class UniversalInvolutiveBasisContainer; // fwd decl -- defined in TmpUniversalInvolutiveBasisContainer.H

    const std::vector<RingElem>& JanetBasis(const ideal& I);
    bool IsDeltaRegular(const ideal& I);
    bool IsMonomial(const ideal& I);
    bool IsHomogeneous(const ideal& I);
    bool IsCohenMacaulay(const ideal& I);
    std::map<PPMonoidElem, std::vector<bool> > MultVars(const ideal& I);
    std::map<PPMonoidElem, std::vector<bool> > NonMultVars(const ideal& I);
    RingElem HilbertPol(const ideal& I, ConstRefRingElem var);
    RingElem HilbertSeries(const ideal& I, ConstRefRingElem var);
    FGModule FirstSyzygy(const ideal& I);
    long Dimension(const ideal& I);
    std::vector<std::pair<PPMonoidElem, std::vector<bool> > > ComplementaryDecomposition(const ideal& I);
    long Depth(const ideal& I);
    long ProjDim(const ideal& I);
    std::vector<RingElem> Socle(const ideal& I);
    std::map<std::pair<long, long>, long> ExtremalBettiNumbers(const ideal& I);
    std::vector<RingElem> RegularSequence(const ideal& I);
    std::vector<RingElem> MaximalStronglyIndependentSet(const ideal& I);
    long Regularity(const ideal& I);
    long Satiety(const ideal& I);
    std::vector<RingElem> Saturation(const ideal& I);
  }

  

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/SparsePolyOps-involutive.H,v 1.3 2024/03/18 20:57:28 abbott Exp $
// $Log: SparsePolyOps-involutive.H,v $
// Revision 1.3  2024/03/18 20:57:28  abbott
// Summary: Removed const from return value type (redmine 1606)
//
// Revision 1.2  2022/02/18 14:11:49  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.1  2018/03/15 14:56:54  bigatti
// -- from SparsePolyRing.H
//

#endif
