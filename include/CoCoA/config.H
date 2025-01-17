#ifndef CoCoA_config_H
#define CoCoA_config_H

//   Copyright (c)  2005,2013  John Abbott,  Anna M. Bigatti

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


#include <climits>

namespace CoCoA
{

  // >>> EXPONENTS <<<

  // This is the type used internally for "small" exponents in some PPMonoids
  // (PPMonoidEv, PPMonoidEvOv).  Also implies the number of bits used for a
  // "compressed entry" in an order vector (see OrdvArith).
  typedef unsigned int SmallExponent_t;


  // -------------------------------------------------------
  // >>> ROUNDING HALVES <<<

  // The way CoCoALib rounds halves is determined by the fns RoundDiv in
  // IntOperations.C: one for machine ints, and one for BigInts.


  // -------------------------------------------------------
  // >>> INTEGER FACTORIZATION <<<

  // Trial division limit for factor(BigInt) defined in NumTheory-factor.C
  constexpr long FactorBigIntTrialLimit = 33554432; // arbitrary limit (happens to be 2^25)


  // -------------------------------------------------------
  // >>> INTEGER OVERFLOW <<<

  // Some functions check for BigInt overflow: power & factorial.
  // This constant sets the max bit-size allowed in result value.
  // Need to check whether platform is 32-bit (hence #if construct).
#if LONG_MAX/(1024L*1024L*1024L) >= 4
  constexpr long OVERFLOW_BITS = 4*1024L*1024L*1024L; // more than 32 bits
#else
  constexpr long OVERFLOW_BITS = 1024L*1024L*1024L; // 32-bitter
#endif

} // end of namespace CoCoA

#endif
