#ifndef CoCoA_factor_H
#define CoCoA_factor_H

//   Copyright (c)  2006,2009,2012  John Abbott,  Anna M. Bigatti

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


#include "CoCoA/factorization.H"


namespace CoCoA
{
  class RingElem;         // fwd decl -- defined in ring.H
  class RingElemAlias;    // fwd decl -- defined in ring.H
  typedef const RingElemAlias& ConstRefRingElem;

  //TEMPORARILY AVAILABLE!!!!
  RingElem GCD_DUPFF(ConstRefRingElem f, ConstRefRingElem g);
  RingElem GCD_DMPZ(ConstRefRingElem f, ConstRefRingElem g);

  factorization<RingElem> factor(ConstRefRingElem f);
  bool IsIrredPoly(ConstRefRingElem f); // ***TEMPORARY*** only for PolyRingBase::myIsIrred

  RingElem RadicalOfPoly(ConstRefRingElem f);
  factorization<RingElem> SqFreeFactor(ConstRefRingElem f);
  factorization<RingElem> ContentFreeFactor(ConstRefRingElem f);

  RingElem PalindromicFactor(ConstRefRingElem f); // just gcd(f, reverse(f)); so x-1 is "palindromic"!!!

  //---- defines in factor-AlgExtn.C ---------------------------------
  factorization<RingElem> factor_AlgExtn(ConstRefRingElem f);

}

#endif
