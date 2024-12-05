#ifndef CoCoA_TmpImplicit_H
#define CoCoA_TmpImplicit_H

//   Copyright (c)  2014  John Abbott, and Anna M. Bigatti

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

#include <vector>

namespace CoCoA
{

  RingElem ImplicitDirect(const std::vector<RingElem>& ParamDescr);
  RingElem ImplicitDirectLPP(const std::vector<RingElem>& ParamDescr);
  RingElem ImplicitDirectLPP2(const std::vector<RingElem>& ParamDescr);
  RingElem ImplicitDirectOrd2(const std::vector<RingElem>& ParamDescr);
  RingElem ImplicitDirectWithCond(const std::vector<RingElem>& ParamDescr, const std::vector<RingElem>& relations);
  RingElem ImplicitDirectWithCondLPP(const std::vector<RingElem>& ParamDescr, const std::vector<RingElem>& relations);

  RingElem ImplicitByPoints(const std::vector<RingElem>& ParamDescr);
  RingElem ImplicitByPoints2(const std::vector<RingElem>& ParamDescr);
  RingElem ImplicitByPoints3(const std::vector<RingElem>& ParamDescr);

  RingElem SliceCore(const std::vector<RingElem>& ParameDescr,
                     long RecDepth,
                     const std::string& FinalCalls);
  RingElem SliceCoreQQ(const std::vector<RingElem>& ParamDescr,
                       long RecDepth,
                       const std::string& FinalCall);
  

} // end of namespace CoCoA

#endif
