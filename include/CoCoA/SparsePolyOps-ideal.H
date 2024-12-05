#ifndef CoCoA_SparsePolyOps_ideal_H
#define CoCoA_SparsePolyOps_ideal_H

//   Copyright (c)  2017-2018  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/CpuTimeLimit.H"
#include "CoCoA/MatrixForOrdering.H"

#include <vector>

namespace CoCoA
{

  class ideal;     // fwd decl -- defined in ideal.H
  class RingElem;  // fwd decl -- defined in ring.H
  class PPMonoidElem;  // fwd decl -- defined in PPMonoid.H
  
  //-- Special functions on ideals owned by SparsePolyRing ----------
  //-- defined in SparsePolyOps-ideal.C

  bool IsZeroDim(const ideal& I);
  bool IsHomog(const ideal& I);
  const std::vector<RingElem>& GBasis(const ideal& I, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
  std::vector<RingElem> GBasisTrunc(const ideal& I, const long TruncDeg, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
  const std::vector<RingElem>& GBasisByHomog(const ideal& I, const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit());
  std::vector<RingElem> GBasisSelfSatCore(const ideal& I);
  std::vector<RingElem> GBasisRealSolve(const ideal& I);
  const std::vector<RingElem>& ReducedGBasis(const ideal& I);
  const std::vector<RingElem>& MinGens(const ideal& I);
  std::vector<ideal> PrimaryDecomposition(const ideal& I);
  ideal LT(const ideal& I);
  ideal LF(const ideal& I);
  ideal homog(const ideal& I, ConstRefRingElem x);
  ideal IdealOfGBasis(const ideal& I);
  ideal IdealOfMinGens(const ideal& I);
  ideal elim(const ideal& I, const std::vector<RingElem>& ElimIndets);
  std::vector<PPMonoidElem> QuotientBasis(const ideal& I);
  std::vector<PPMonoidElem> QuotientBasisSorted(const ideal& I);
  bool HasGBasis(const ideal& I);
  bool AreGensMonomial(const ideal& I);
  bool AreGensSqFreeMonomial(const ideal& I);
  //  bool3 AreGensMonomial3(const ideal& I);
  //  bool3 AreGensSqFreeMonomial3(const ideal& I);
  void SetGBasisAsGens(const ideal& I);
  RingElem DenSigma(const ideal& I);
  bool IsSigmaGoodPrime(const BigInt& p, const ideal& I);
  bool IsSigmaGoodPrime(const long p, const ideal& I);

  bool IsRadical_0dim(const ideal& I); // will be IsRadical
  //  ideal radical_0dim(const ideal& I);
  ideal radical_general(const ideal& I);
  ideal radical(const ideal& I);

  bool IsPrimary(const ideal& I);

} // end of namespace CoCoA

#endif
