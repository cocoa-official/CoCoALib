#ifndef CoCoA_TmpMonomialFns_H
#define CoCoA_TmpMonomialFns_H

//   Copyright (c)  2008-2009  Anna Bigatti and Eduardo Sáenz-de-Cabezón

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


// Source code for functions on monomial ideals

#include <vector>
// using vector;

/**********************************************************/
namespace CoCoA
{
  class PPVector;       // forward declaration -- defined in PPVector.H
  class PPMonoidElem;   // forward declaration -- defined in PPMonoid.H

  void support(std::vector<long>& , const PPMonoidElem& );
  bool IsIrreducible(const PPVector& );
  bool IsPrime(const PPVector& );
  bool IsPrimary(const PPVector& );
  void ColonIdeal(PPVector& , const PPVector& , const PPVector& );
}

#endif /* TmpMonomialFunctions_H */

