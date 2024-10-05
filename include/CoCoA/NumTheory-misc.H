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


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/NumTheory-misc.H,v 1.12 2024/01/25 08:38:36 abbott Exp $
// $Log: NumTheory-misc.H,v $
// Revision 1.12  2024/01/25 08:38:36  abbott
// Summary: Removed InvTotient iterator, and InvTotientUpto (no longer needed)
//
// Revision 1.11  2023/05/11 19:26:15  abbott
// Summary: New impl and interface to InvTotient
//
// Revision 1.10  2023/03/15 21:31:49  abbott
// Summary: Added new fn InvTotientUpto
//
// Revision 1.9  2023/03/09 22:32:58  abbott
// Summary: Changed BinomialReprShift to take arg by const ref
//
// Revision 1.8  2023/02/27 21:09:08  abbott
// Summary: Revised InvTotientBound (more complicated, but a bit better sometimes)
//
// Revision 1.7  2023/02/23 20:48:04  abbott
// Summary: Added InvTotientBoundUpto for BigInt; revised impls
//
// Revision 1.6  2023/02/16 20:41:57  abbott
// Summary: Added InvTotientIter and InvTotientBoundUpto
//
// Revision 1.5  2023/01/31 12:26:37  abbott
// Summary: Added InvTotient, InvTotientBound
//
// Revision 1.4  2022/11/30 15:15:13  abbott
// Summary: Added MoebiusFn (1st impl)
//
// Revision 1.3  2022/02/18 14:11:46  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.2  2021/02/10 19:32:00  abbott
// Summary: Added noexcept (redmine 1572)
//
// Revision 1.1  2020/01/26 14:14:31  abbott
// Summary: Finished splitting NumTheory into smaller pieces (redming 1161)
//
// Revision 1.47  2019/09/16 17:34:31  abbott
// Summary: Changed name EulerPhi into EulerTotient
//
// Revision 1.46  2019/03/18 11:24:57  abbott
// Summary: Split NumTheory into several smaller files
//
// Revision 1.45  2019/02/28 11:28:06  abbott
// Summary: Added KroneckerSymbol
//
// Revision 1.44  2018/04/20 18:51:25  abbott
// Summary: Changed ctors for BigInt/BigRat from string or from MPZ/MPQ
//
// Revision 1.43  2018/02/27 10:51:07  abbott
// Summary: Moved stuff about primes into NumTheory_prime
//
// Revision 1.42  2018/02/15 17:26:11  abbott
// Summary: Added EratosthenesRange, and PrimeSeq
//
// Revision 1.41  2018/01/17 10:25:18  abbott
// Summary: Changed InvModNoCheck into InvModNoArgCheck
//
// Revision 1.40  2018/01/16 11:41:18  abbott
// Summary: Changed NoThrow into RtnZeroOnError
//
// Revision 1.39  2017/11/08 14:03:30  abbott
// Summary: Added new fn SumOfFactors
//
// Revision 1.38  2017/10/17 15:44:26  abbott
// Summary: Added new fn IsCoprime
//
// Revision 1.37  2016/11/22 14:31:01  abbott
// Summary: Added SmallestNonDivisor
//
// Revision 1.36  2016/10/25 20:54:09  abbott
// Summary: Added new fn IsSqFree (for BigInt and ringelem of PolyRing over field)
//
// Revision 1.35  2016/10/19 13:41:49  abbott
// Summary: Added radical for integers
//
// Revision 1.34  2016/06/29 13:11:37  abbott
// Summary: Added "NoThrow" option to InvMod; much cleaning inside NumTheory.C
//
// Revision 1.33  2015/11/21 19:16:09  abbott
// Summary: Added SimplestBinaryRatBetween
//
// Revision 1.32  2015/11/05 14:16:22  abbott
// Summary: Added InvModNoCheck
//
// Revision 1.31  2015/06/29 13:23:26  abbott
// Summary: Changed name "valuation" --> "FactorMultiplicity"; added new fn "BadMFactor"
// Author: JAA
//
// Revision 1.30  2014/10/28 15:12:07  abbott
// Summary: Renamed modulus --> CombinedModulus, residue --> CombinedResidue (for CRTMill)
// Author: JAA
//
// Revision 1.29  2014/09/16 10:41:41  abbott
// Summary: Added new fn eratosthenes (with doc, example, test)
// Author: JAA
//
// Revision 1.28  2014/08/29 16:04:55  abbott
// Summary: Added optional 3rd arg to myAddInfo (so coprimality check is skipped)
// Author: JAA
//
// Revision 1.27  2014/05/02 13:54:05  abbott
// Summary: Simplified ctor interface for RatReconstruct* (need explicit arg 0 for default behaviour)
// Author: JAA
//
// Revision 1.26  2014/04/15 13:27:19  abbott
// Summary: Changed rtn type of PrimitiveRoot to long (for CoCoA-5/BuiltinOneLiners)
// Author: JAA
//
// Revision 1.25  2014/01/16 16:09:54  abbott
// Added NumPartitions.
//
// Revision 1.24  2013/10/15 16:20:03  abbott
// Added valuation.
//
// Revision 1.23  2013/05/21 14:31:45  abbott
// Added BinomialRepr and BinomialReprShift to CoCoALib and CoCoA-5.
//
// Revision 1.22  2013/05/20 15:47:09  abbott
// Added new fn BinomialRepr (placed in NumTheory).
//
// Revision 1.21  2013/02/26 11:29:17  abbott
// Added impl of RatReconstructWithBounds
//
// Revision 1.20  2013/02/22 22:43:06  abbott
// Changed syntax for getting result from CRTMill.
//
// Revision 1.19  2013/02/22 18:56:50  abbott
// Added feature that RatReconstructByContFrac & RatReconstructByLattice
// ctors accept arg 0 to mean "use default value".
//
// Revision 1.18  2013/02/19 18:48:15  abbott
// Added printing for CRTMill and RatReconstructByContFrac and RatReconstructByLattice.
//
// Revision 1.17  2013/02/15 17:46:00  abbott
// Added RatReconstructByContFrac and RatReconstructByLattice.
//
// Revision 1.16  2012/12/12 18:24:26  abbott
// Added new fn IsFinal for ContFracIter.
//
// Revision 1.15  2012/12/11 17:30:30  abbott
// Changed name from SimplestRationalInInterval to SimplestBigRatBetween.
// Also fixed a bug in the impl.
//
// Revision 1.14  2012/12/05 15:09:24  abbott
// Added new class ContFracApproximant.
// Added new fn SimplestRationalInInterval (NB name will change).
// Some minor cleaning.
//
// Revision 1.13  2012/12/04 20:14:11  abbott
// Added new class CRTMill.
// Improved impl of ContFracIter class.
// Fixed a bug in CFApproxIter class.
//
// Revision 1.12  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.11  2012/05/22 10:03:32  abbott
// Improved indentation/layout.
//
// Revision 1.10  2012/03/16 15:40:12  abbott
// Merged contents of NumTheoryQQ (continued fraction functions) into NumTheory.
// Merged the doc too.
//
// Revision 1.9  2011/11/09 13:47:56  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.8  2011/08/14 15:52:18  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.7  2011/03/23 21:00:46  abbott
// Removed FindPrimRoot from NumTheory.H because it was already
// available as PrimitiveRoot (a better name).
// Updated documentation for NumTheory.
//
// Revision 1.6  2011/03/22 20:17:18  abbott
// Added fn FindPrimRoot.
// Merged impls from obsolescent SmallPrime.C.
//
// Revision 1.5  2011/03/16 13:26:36  abbott
// Removed all "unsigned" from fn interfaces, and many unsigned from inside fn impls.
//
// Revision 1.4  2010/03/03 10:43:34  abbott
// Added PrimitiveRoot for big primes (might be very slow).
// Added MultiplicativeOrderMod (currently very SLUGGY implementation).
//
// Revision 1.3  2009/12/29 22:44:32  abbott
// Removed buggy proxy class ZZ::rtn.
// Consequent changes for function prototypes also in NumTheory.
// Corrected some minor buglets in NumTheory.
//
// Revision 1.2  2009/07/02 16:28:10  abbott
// Fairly extensive change to NumTheory (at least internally and philosophically).
// Improved and cleaned NumTheory.  Added documentation.
// Clarified the exact behaviour of most functions.
//
// Revision 1.1  2009/06/05 12:14:55  abbott
// Major change:
//   created new files NumTheory.H/C  which contain basic number theory operations
//   removed several basic number theory operations from ZZ.H/C
//   removed gcd from MachineInt.H/C
//   changed names of some basic fns:
//      IsPPrime -> IsProbPrime
//      invmod -> InvMod    (changed signature too)
//      powermod -> PowerMod  (changed signature too)
//   added new fns
//      NextProbPrime & PrevProbPrime
//   consequent changes to other code and tests and examples
//
//

#endif
