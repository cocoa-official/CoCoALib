#ifndef CoCoA_NumTheory_misc_H
#define CoCoA_NumTheory_misc_H

//   Copyright (c)  2009,2012,2020,2023  John Abbott,  Anna M. Bigatti
//   Some contributions from Nico Mexis

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


#include <vector>
// using std::vector;

namespace CoCoA
{
  
  class MachineInt; // fwd decl -- defined in MachineInt.H
  class BigInt;     // fwd decl -- defined in BigInt.H


  //--------------------------------------------
  // Euler's "totient" function, & inverse.

  long   EulerTotient(const MachineInt& n); // NOT noexcept with curr impl
  BigInt EulerTotient(const BigInt& N);


  // Compute all pre-images, or just sqfree preimages
  enum class InvTotientMode { AllPreimages, SqFreePreimages };

  std::vector<long> InvTotient(long n, InvTotientMode mode=InvTotientMode::AllPreimages);

  // Upper bounds for largest pre-image
  BigInt InvTotientBound(const BigInt& N);
  unsigned long InvTotientBound_ulong(const unsigned long n); // 0 means no pre-image; 1 means "overflow"

  // Upper bound for largest pre-image of any k <= n
  BigInt InvTotientBoundUpto(const BigInt& N);
  BigInt InvTotientBoundUpto(const unsigned long n);
  unsigned long InvTotientBoundUpto_ulong(const unsigned long n); // 0 means no pre-image; 1 means "overflow"

  
  bool IsEulerTotient(long n, long phi_n); // true iff phi_n == EulerTotient(n); error unless both args positive.


  //--------------------------------------------
  // Moebius function: 0 if n is not sqfr; o/w (-1)^num_prime_factors

  long MoebiusFn(const MachineInt& n);


  //-------------------------------------------------------
  // Repr of an int as sum of binomials:

  std::vector<BigInt> BinomialRepr(BigInt N, long r);
  BigInt BinomialReprShift(const BigInt& N, long r, long shift1, long shift2);

} // end of namespace CoCoA

#endif
