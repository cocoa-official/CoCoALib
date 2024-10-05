//   Copyright (c)  2014  John Abbott and Anna M. Bigatti

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

#include "CoCoA/QuasiPoly.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/VectorOps.H"
#include "CoCoA/error.H"
#include "CoCoA/utils.H"


namespace CoCoA
{

  QuasiPoly::QuasiPoly(const std::vector<RingElem>& v):
      myConstituents(v)
  {
    // Simply check that the input v is valid...
    if (v.empty()) CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
    const ring& R = owner(v[0]);
    const long n = len(v);
    for (long i=1; i < n; ++i)
      if (owner(v[i]) != R) CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (!IsPolyRing(R)) CoCoA_THROW_ERROR1(ERR::ReqPolyRing);
  }


  const std::vector<RingElem>& constituents(const QuasiPoly& p)
  { return p.myConstituents; }


  RingElem QuasiPoly::operator()(const MachineInt& n) const
  { return operator()(BigInt(n)); }

  
  RingElem QuasiPoly::operator()(const BigInt& N) const
  {
    const long i = ConvertTo<long>(LeastNNegRemainder(N, len(myConstituents)));
    const PolyRing P = owner(myConstituents[0]);
    const RingHom EvalAtN = EvalHom(P, N);
    return EvalAtN(myConstituents[i]);
  }

  
  std::ostream& operator<<(std::ostream& out, const QuasiPoly& p)
  {
    if (!out) return out;  // short-cut for bad ostreams
    out << "QuasiPoly(" << constituents(p) << ")";
    return out;
  }

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/QuasiPoly.C,v 1.11 2024/07/15 17:10:00 bigatti Exp $
// $Log: QuasiPoly.C,v $
// Revision 1.11  2024/07/15 17:10:00  bigatti
// Summary: Error codes:  added CoCoA_THROW_ERROR1(ErrID);
// Removed BadRowIndex, BadColIndex;  Renamed  ERR::Empty --> ERR::ReqNonEmpty
//
// Revision 1.10  2024/07/02 20:42:07  abbott
// Summary: Renaming of errors (redmine 308)
//
// Revision 1.9  2022/02/18 14:11:56  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.8  2021/01/07 15:16:52  abbott
// Summary: Corrected copyright
//
// Revision 1.7  2020/06/17 15:49:25  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.6  2018/05/18 12:15:04  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.5  2018/05/17 15:40:05  bigatti
// -- renamed VectorOperations --> VectorOps
//
// Revision 1.4  2016/11/11 14:15:33  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.3  2014/07/31 14:45:18  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.2  2014/07/14 13:16:46  abbott
// Summary: Minor cleaning
// Author: JAA
//
// Revision 1.1  2014/07/14 10:04:29  abbott
// Summary: Added new impl of quasipolynomials
// Author: JAA
//
