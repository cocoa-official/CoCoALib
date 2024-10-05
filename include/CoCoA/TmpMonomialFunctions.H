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

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/TmpMonomialFunctions.H,v 1.6 2024/05/03 12:08:26 abbott Exp $
// $Log: TmpMonomialFunctions.H,v $
// Revision 1.6  2024/05/03 12:08:26  abbott
// Summary: Improved include guard (related to redmine 1811)
//
// Revision 1.5  2024/03/28 22:09:41  abbott
// Summary: Corrected include guards (redmine 1811)
//
// Revision 1.4  2022/02/18 14:11:50  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.3  2012/02/01 13:28:53  abbott
// Changed 1st arg of support into vector<long> (was vector<size_t>)
//
// Revision 1.2  2010/02/04 10:54:49  bigatti
// -- fix: moved forward calss declaration inside namespace CoCoA
//
// Revision 1.1  2010/02/03 18:41:29  bigatti
// -- first import
//

#endif /* TmpMonomialFunctions_H */

