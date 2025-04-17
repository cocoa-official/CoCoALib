#ifndef CoCoA_GOperations_H
#define CoCoA_GOperations_H

//   Copyright (c)  2005  John Abbott,  Anna M. Bigatti
//   Original author: 2005  Massimo Caboara

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

#include "CoCoA/TmpGTypes.H"
#include "CoCoA/TmpGReductor.H"

//#include <vector>
// using vector;

// Operations depending on Groebner Bases


namespace CoCoA
{

  class ConstRefPPMonoidElem; // forward declaration -- defined in PPMonoid.H
  class RingElemAlias;        // forward declaration -- defined in ring.H
  typedef const RingElemAlias& ConstRefRingElem;     // forward declaration -- defined in ring.H
  class ModuleElem;           // forward declaration -- defined in PPMonoid.H
  class RingElem;             // forward declaration -- defined in ring.H

  void ComputeGBasis(PolyList& outGB, PolyList& outMinGens, const PolyList& , const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
  void ComputeGBasisTrunc(PolyList& outGB, PolyList& outMinGens, long& TruncDeg, const PolyList& , const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
  void ComputeGBasisSelfSatCore(PolyList& outGB, const PolyList& inGens, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
  void ComputeGBasisRealSolve(PolyList& outGB, const PolyList& inGens, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
  //  void ComputeSATMixGBasis(PolyList& outGB, const PolyList& inGens);

  void ComputeGBasis(VectorList& outGB, VectorList& outMinGens, const VectorList& inGens, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
  void ComputeGBasisFrameWork(PolyList&, const PolyList&);

// The second argument is the list of indices of the
// indeterminates that will be eliminated
void ComputeElim(VectorList&,const VectorList&,ConstRefPPMonoidElem, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());///< Not yet implemented
void ComputeElim(PolyList&,const PolyList&,ConstRefPPMonoidElem, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());

  void ComputeSyz(VectorList&, const FreeModule&, const VectorList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
  void ComputeSyz(VectorList&, const FreeModule&, const PolyList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());

void ComputeIntersection(VectorList&,const VectorList&,const VectorList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
void ComputeIntersection(PolyList&,const PolyList&,const PolyList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());

void ComputeColonByPrincipal(PolyList&, const VectorList&,const VectorList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
void ComputeColonByPrincipal(PolyList&, const PolyList&, ConstRefRingElem f, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
void ComputeColonByPrincipal(VectorList&, const VectorList&,const PolyList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());

void ComputeColon(PolyList&,const VectorList&,const VectorList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
void ComputeColon(PolyList&,const PolyList&,const PolyList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit() );
void ComputeColon(VectorList&,const VectorList&,const PolyList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());

void ComputeSaturationByPrincipal(VectorList&,const VectorList&,const PolyList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
void ComputeSaturationByPrincipal(PolyList&, const PolyList&, ConstRefRingElem, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());

void ComputeSaturation(PolyList&,const PolyList&,const PolyList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
void ComputeSaturation(VectorList&,const VectorList&,const PolyList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());

void ComputeHomogenization(VectorList&,const VectorList&,const PolyList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());///< Not yet implemented
void ComputeHomogenization(PolyList&,const PolyList&,const PolyList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());///< Not yet implemented

bool RadicalMembership(const PolyList&,ConstRefRingElem, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());

void ComputeLT(VectorList&,const VectorList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());///< Not yet implemented
void ComputeLT(PolyList&,const PolyList&, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());

  RingElem ComputeElimFirst(const PolyList& inPL, ConstRefPPMonoidElem inds, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());

}// end namespace cocoa

#endif
