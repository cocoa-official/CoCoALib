#ifndef CoCoA_NumTheory_factor_H
#define CoCoA_NumTheory_factor_H

//   Copyright (c)  2009,2012  John Abbott,  Anna M. Bigatti

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

#include "CoCoA/CpuTimeLimit.H"

#include <vector>
//using std::vector;

namespace CoCoA
{
  
  class MachineInt; // fwd decl -- defined in MachineInt.H


  long radical(const MachineInt& n);  // NOT noexcept with curr impl
  BigInt radical(const BigInt& N);

  factorization<long>   factor(const MachineInt& N);
  factorization<BigInt> factor(const BigInt& N); // TimeOut???
  factorization<BigInt> factor(const BigInt& N, const CpuTimeLimit& CheckTime);

  factorization<long>   factor_TrialDiv(const MachineInt& N, const MachineInt& TrialLimit);
  factorization<BigInt> factor_TrialDiv(const BigInt& N,     const MachineInt& TrialLimit);
  factorization<BigInt> factor_TrialDiv(const BigInt& N,     const BigInt& TrialLimit);

  class PollardRhoSeq
  {
  public:
    explicit PollardRhoSeq(const BigInt& N, long StartVal=2, long incr=1);
    explicit PollardRhoSeq(const PollardRhoSeq&) = default;
    void myAdvance(long n); // advance by AT LEAST n iters
    PollardRhoSeq& myReduceN(const BigInt& fac); // PRS found a factor; fac is "coprime part" of orig N
    friend bool IsEnded(const PollardRhoSeq& PRS) noexcept;
    friend const BigInt& GetN(const PollardRhoSeq& PRS) noexcept;
    friend const BigInt& GetFactor(const PollardRhoSeq& PRS) noexcept;
    friend long GetNumIters(const PollardRhoSeq& PRS) noexcept;
    friend std::ostream& operator<<(std::ostream& out, const PollardRhoSeq& PRS);
  private: // data members
    long myNumIters; // counts how many steps
    long myAnchorIterNum;
    BigInt myAnchorVal;
    BigInt myN;
    BigInt myCurrVal;
    BigInt myGCD;
    long myStartVal;
    long myIncr;
    long myBlockSize; // how many 
  };

  bool IsEnded(const PollardRhoSeq& PRS) noexcept;
  const BigInt& GetN(const PollardRhoSeq& PRS) noexcept;
  const BigInt& GetFactor(const PollardRhoSeq& PRS) noexcept;
  long GetNumIters(const PollardRhoSeq& PRS) noexcept;
  std::ostream& operator<<(std::ostream& out, const PollardRhoSeq& PRS);
  factorization<BigInt> factor_PollardRho(BigInt N, long NumIters);


  std::vector<long> AllFactors(long n);
  BigInt SumOfFactors(const MachineInt& N, long k=1); // sum of k-th powers of all factors
//  BigInt SumOfFactors(const BigInt& N, long k=1);

  long SmallestNonDivisor(const MachineInt& N); // NOT noexcept -- arg check
  long SmallestNonDivisor(const BigInt& N);

  bool IsSqFree(const MachineInt& n); // NOT noexcept with curr impl
  bool3 IsSqFree(BigInt N);

  long FactorMultiplicity(const MachineInt& p, const MachineInt& n); // p need not be prime!
  long FactorMultiplicity(const MachineInt& p, BigInt N);            // p need not be prime!
  long FactorMultiplicity(const BigInt& p, BigInt N);                // p need not be prime!

  // CoprimeFactor(N,b) = N/gcd(N,b^infty)
  long CoprimeFactor(long N, long b);
  BigInt CoprimeFactor(BigInt N, long b);
  BigInt CoprimeFactor(BigInt N, const BigInt& b);

} // end of namespace CoCoA

#endif
