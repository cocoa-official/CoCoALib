//   Copyright (c)  2020  John Abbott and Anna M. Bigatti
//   Authors:  2005-2020  John Abbott

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


// Source code for abstract class SparsePolyRing and friends

#include "CoCoA/SparsePolyIter.H"

#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/error.H"

namespace CoCoA
{

  SparsePolyIter BeginIter(ConstRefRingElem f)
  {
    if (!IsSparsePolyRing(owner(f)))
      CoCoA_THROW_ERROR2(ERR::ReqSparsePolyRing, "ring of f");
    return SparsePolyRingPtr(owner(f))->myBeginIter(raw(f));
  }


  SparsePolyIter EndIter(ConstRefRingElem f)
  {
    if (!IsSparsePolyRing(owner(f)))
      CoCoA_THROW_ERROR2(ERR::ReqSparsePolyRing, "ring of f");
    return SparsePolyRingPtr(owner(f))->myEndIter(raw(f));
  }

} // end of namespace CoCoA
