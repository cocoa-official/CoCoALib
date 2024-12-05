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
