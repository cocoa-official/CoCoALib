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



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/SparsePolyOps-ideal.H,v 1.18 2024/05/30 09:41:06 bigatti Exp $
// $Log: SparsePolyOps-ideal.H,v $
// Revision 1.18  2024/05/30 09:41:06  bigatti
// Summary: removed MinGens_almost (redmine #1740)
//
// Revision 1.17  2024/05/28 12:53:40  bigatti
// Summary: added GBasisTrunc, myGBasisTrunc, ComputeGBasisTrunc
//
// Revision 1.16  2024/05/24 14:38:41  bigatti
// Summary: added MinGens_almost
//
// Revision 1.15  2024/03/15 15:33:53  bigatti
// Summary: added function elim
//
// Revision 1.14  2024/02/06 10:59:02  bigatti
// Summary: removed wrong comment (radical_general)
//
// Revision 1.13  2024/02/05 13:09:37  bigatti
// Summary: now radical is in C++
//
// Revision 1.12  2023/07/04 09:16:42  abbott
// Summary: Changed TimeOut to Timeout
//
// Revision 1.11  2022/02/18 14:11:49  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.10  2019/10/03 13:33:56  bigatti
// -- implemented radical for monomial ideals (and used where useful)
//
// Revision 1.9  2018/08/06 09:38:29  bigatti
// -- renamed GBasisViaHomog --> GBasisByHomog
//
// Revision 1.8  2018/08/06 08:57:48  bigatti
// -- added timeout for GBasisViaHomog
// -- now using GBasisViaHomog in IsPrimary and radical
//
// Revision 1.7  2018/08/05 16:31:25  bigatti
// -- added GBasisViaHomog
//
// Revision 1.6  2018/06/27 12:15:18  abbott
// Summary: Renamed RealSolveCore to RealSolve
//
// Revision 1.5  2018/06/27 08:50:39  abbott
// Summary: Revised to work with new CpuTimeLimit
//
// Revision 1.4  2018/04/06 15:42:26  bigatti
// -- copyright year
//
// Revision 1.3  2018/03/20 15:37:28  bigatti
// -- changed Radical0dim --> radical_0dim (etc)
//
// Revision 1.2  2018/03/20 11:14:12  bigatti
// -- added IsSigmaGoodPrime, DenSigma, Radical0dim, IsRadical0dim, IsPrimary
//
// Revision 1.1  2018/03/15 14:56:54  bigatti
// -- from SparsePolyRing.H
//

#endif
