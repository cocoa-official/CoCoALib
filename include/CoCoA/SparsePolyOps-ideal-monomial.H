#ifndef CoCoA_ideal_monomial_H
#define CoCoA_ideal_monomial_H

//   Copyright (c)  2011-2019  John Abbott, Anna M. Bigatti

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


// Header file for operations on monomial ideals


#include <iosfwd>
//using std::ostream;
#include <vector>
//using std::vector;
#include <memory>
// using std::auto_ptr;


namespace CoCoA
{

  class ideal;        // fwd decl -- defined in ideal.H 

// the member functions still depend on SparsePolyRing::IdealImpl
// therefore they are declared in SparsePolyRing.H

//       const std::vector<RingElem>& myGBasis_MonId() const;
//       bool myTestIsRadical_MonId() const;
//       ideal myRadical_MonId() const;
//       std::vector<ideal> myPrimaryDecomposition_MonId() const;
//       void myMul_MonId(const ideal&);
//       void myIntersect_MonId(const ideal&);
//       void myColon_MonId(const ideal&);
//       void myElim_MonId(const std::vector<RingElem>& ELimIndets);

  ideal AlexanderDual(const ideal& I);

} // end of namespace CoCoA

#endif
