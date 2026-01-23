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

  // "2" because they "return" 2 arguments GB_out and MinGens_out
  void ComputeGBasis2(PolyList& GB_out, PolyList& MinGens_out,
                      const PolyList& , const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
  void ComputeGBasisTrunc2(PolyList& GB_out, PolyList& MinGens_out,
                           long& TruncDeg, const PolyList& , const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
  void ComputeGBasis2(VectorList& GB_out, VectorList& MinGens_out,
                      const VectorList& G_in, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());

  PolyList ComputeGBasisSelfSatCore(const PolyList& G, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
  PolyList ComputeGBasisRealSolve(const PolyList& G, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());

  // second arg is product of the indets to be eliminated
  VectorList ComputeElim(const VectorList& G, ConstRefPPMonoidElem IndetProd,
                         const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());///< Not yet implemented
  PolyList ComputeElim(const PolyList& G, ConstRefPPMonoidElem IndetProd,
                       const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
  RingElem ComputeElimFirst(const PolyList& G, ConstRefPPMonoidElem IndetProd,
                            const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());

  VectorList ComputeSyz(const FreeModule& FM, const VectorList& G,
                        const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
  VectorList ComputeSyz(const FreeModule& FM, const PolyList& G,
                        const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());

  VectorList ComputeIntersection(const VectorList& G1, const VectorList& G2,
                                 const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
  PolyList ComputeIntersection(const PolyList& G1, const PolyList& G2,
                               const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());

  PolyList ComputeColon(const VectorList& G1, const VectorList& G2,
                        const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
  PolyList ComputeColon(const PolyList& G1, const PolyList& G2,
                        const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit() );
  VectorList ComputeColon(const VectorList& G1, const PolyList& G2,
                          const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());

  PolyList ComputeSaturation(const PolyList& G1, const PolyList& G2,
                             const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
  VectorList ComputeSaturation(const VectorList& G1, const PolyList& G2,
                               const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());

  VectorList ComputeHomogenization(const VectorList&, const PolyList&,
                                   const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());///< Not yet implemented
  PolyList ComputeHomogenization(const PolyList&, const PolyList&,
                                 const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());///< Not yet implemented

  bool RadicalMembership(const PolyList&,ConstRefRingElem,
                         const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());


}// end namespace cocoa

#endif
