//   Copyright (c)  2015  John Abbott,  Anna M. Bigatti

//   This file is part of the source of CoCoALib, the CoCoA Library.

//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.

//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.

//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.


#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/error.H"
#include "CoCoA/NumTheory-modular.H"
#include "CoCoA/NumTheory-factor.H"
#include "CoCoA/NumTheory-prime.H"


#include <algorithm>
using std::min;
#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;


namespace CoCoA
{

  // For test_radical:
  // ad hoc impl -- no arg checks!
  long pwr(long base, int exp)
  {
    long ans = 1;
    for (int i=0; i < exp; ++i)
      ans *= base;
    return ans;
  }


  //----------------------------------------------------------------------
  // This test checks that radical gives the expected answer over a
  // range of sall numbers with prime factors 2,3,5,7.
  //----------------------------------------------------------------------
  void test_radical()
  {
    const long p1 = 2;
    const long p2 = 3;
    const long p3 = 5;
    const long p4 = 7;
    for (int e1 = 0; e1 < 4; ++e1)    
    for (int e2 = 0; e2 < 4; ++e2)    
    for (int e3 = 0; e3 < 4; ++e3)    
    for (int e4 = 0; e4 < 4; ++e4)    
    {
      const long N = pwr(p1,e1) * pwr(p2,e2) * pwr(p3,e3) * pwr(p4,e4);
      const long radN = pwr(p1,min(e1,1)) * pwr(p2,min(e2,1)) * pwr(p3,min(e3,1)) * pwr(p4,min(e4,1));
      CoCoA_ASSERT_ALWAYS(radical(N) == radN);
      CoCoA_ASSERT_ALWAYS(radical(-N) == radN);
      CoCoA_ASSERT_ALWAYS(radical(radical(N)) == radN);
    }
  }


  //----------------------------------------------------------------------
  // This test checks that IsSqFree gives the expected answer over the
  // range 1 to 100.
  //----------------------------------------------------------------------
  void test_IsSqFree()
  {
    for (int i=1; i < 100; ++i)
    {
      const factorization<long> FacInfo = factor(i);
      const vector<long>& mult = FacInfo.myMultiplicities();
      int MaxExp = 0;
      for (int j=0; j < len(mult); ++j)
        if (mult[j] > MaxExp) MaxExp = mult[j];
      const bool ans1 = IsSqFree(i);
      const bool3 ans2 = IsSqFree(BigInt(i));
      const bool ans3 = IsSqFree(-i);
      const bool3 ans4 = IsSqFree(BigInt(-i));
      CoCoA_ASSERT_ALWAYS(ans1 == (MaxExp <= 1));
      CoCoA_ASSERT_ALWAYS(ans1 == IsTrue3(ans2) && ans1 == !IsFalse3(ans2));
      CoCoA_ASSERT_ALWAYS(ans1 == ans3);
      CoCoA_ASSERT_ALWAYS(ans1 == IsTrue3(ans4) && ans1 == !IsFalse3(ans4));
    }

  }


  //-------------------------------------------------------
  // Test the various "prime sequence" implementations
  //-------------------------------------------------------
  void test_PrimeSeq()
  {
    // FastFinitePrimeSeq is really just table-lookup.
    // Check the sequence is coherent: starts with 2, and each
    // successive entry is NextPrime of previous entry.
    FastFinitePrimeSeq seq3;
    long p3 = *seq3;
    CoCoA_ASSERT_ALWAYS(p3 == 2);
    while (!IsEnded(++seq3))
    {
      CoCoA_ASSERT_ALWAYS(*seq3 == NextPrime(p3));
      p3 = *seq3;
    }


    // FastMostlyPrimeSeq is same as FastFinitePrimeSeq up to index 82024
    // then it generates some composite numbers (but with no factor < 29)
    FastMostlyPrimeSeq seq0;
    for (int i=0; i < 100000; ++i)
    {
      const long n = *seq0;
      if (n < 25) { CoCoA_ASSERT_ALWAYS(IsPrime(n)); ++seq0; continue; }
      CoCoA_ASSERT_ALWAYS(n%2 != 0);
      CoCoA_ASSERT_ALWAYS(n%3 != 0);
      CoCoA_ASSERT_ALWAYS(n%5 != 0);
      CoCoA_ASSERT_ALWAYS(n%7 != 0);
      CoCoA_ASSERT_ALWAYS(n%11 != 0);
      CoCoA_ASSERT_ALWAYS(n%13 != 0);
      CoCoA_ASSERT_ALWAYS(n%17 != 0);
      CoCoA_ASSERT_ALWAYS(n%19 != 0);
      CoCoA_ASSERT_ALWAYS(n%23 != 0);
      ++seq0;
    }


    // PrimeSeq generates the sequence of primes starting from 2.
    // Check the sequence is coherent: starts with 2, and each
    // successive entry is NextPrime of previous entry.
    // Internally it "changes gear" at around 82024; check for
    // correct behaviour around this point.
    PrimeSeq seq1;
    CoCoA_ASSERT_ALWAYS(*seq1 == 2);
    for (int i=0; i < 82000; ++i)
      ++seq1;
    long p1 = *seq1;
    for (int i=0; i < 100; ++i)
    {
      ++seq1;
      CoCoA_ASSERT_ALWAYS(*seq1 == NextPrime(p1));
      p1 = *seq1;
    }

    // PrimeSeqForCRT generates "large" small primes.
    // It uses a table for the first 50000ish entries then
    // generates on the fly.  We check for a smooth transition.
    PrimeSeqForCRT seq2;
    for (int i=0; i < 50000; ++i)
      ++seq2;
    long p2 = *seq2;
    for (int i=0; i < 100; ++i)
    {
      ++seq2;
      CoCoA_ASSERT_ALWAYS(*seq2 == NextPrime(p2));
      p2 = *seq2;
    }
  }


  //----------------------------------------------------------------------
  // Next test checks the function KroneckerSymbol
  // (for testing whether a value is a quadratic residue).
  // First 2 auxiliary fns.
  //----------------------------------------------------------------------

  // KroneckerSymbol: Check that sum of all values is 0 if p is odd
  void Kronecker_CheckPrime1(int p)
  {
    if (p == 2) return;
    long sum = 0;
    for (long r=0; r < p; ++r)
      sum += KroneckerSymbol(r, p);
    CoCoA_ASSERT_ALWAYS(sum == 0);

    for (long r=0; r < p; ++r)
      sum += KroneckerSymbol(r, BigInt(p));
    CoCoA_ASSERT_ALWAYS(sum == 0);

    for (long r=0; r < p; ++r)
      sum += KroneckerSymbol(BigInt(r), p);
    CoCoA_ASSERT_ALWAYS(sum == 0);

    for (long r=0; r < p; ++r)
      sum += KroneckerSymbol(BigInt(r), BigInt(p));
    CoCoA_ASSERT_ALWAYS(sum == 0);
  }


  // KroneckerSymbol: Check that 0 gives 0, and that all squares give 1
  void Kronecker_CheckPrime2(int p)
  {
    CoCoA_ASSERT_ALWAYS(KroneckerSymbol(0, p) == 0);
    for (long r=1; r < p; ++r)
      CoCoA_ASSERT_ALWAYS(KroneckerSymbol(r*r, p) == 1);

    CoCoA_ASSERT_ALWAYS(KroneckerSymbol(0, BigInt(p)) == 0);
    for (long r=1; r < p; ++r)
      CoCoA_ASSERT_ALWAYS(KroneckerSymbol(r*r, BigInt(p)) == 1);

    CoCoA_ASSERT_ALWAYS(KroneckerSymbol(BigInt(0), p) == 0);
    for (long r=1; r < p; ++r)
      CoCoA_ASSERT_ALWAYS(KroneckerSymbol(BigInt(r*r), p) == 1);

    CoCoA_ASSERT_ALWAYS(KroneckerSymbol(BigInt(0), BigInt(p)) == 0);
    for (long r=1; r < p; ++r)
      CoCoA_ASSERT_ALWAYS(KroneckerSymbol(BigInt(r*r), BigInt(p)) == 1);
  }


  void test_kronecker()
  {
    Kronecker_CheckPrime1(2);
    Kronecker_CheckPrime2(2);

    Kronecker_CheckPrime1(3);
    Kronecker_CheckPrime2(3);

    Kronecker_CheckPrime1(5);
    Kronecker_CheckPrime2(5);

    Kronecker_CheckPrime1(101);
    Kronecker_CheckPrime2(101);

    Kronecker_CheckPrime1(32003);
    Kronecker_CheckPrime2(32003);
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    test_radical();
    test_IsSqFree();
    test_PrimeSeq();
    test_kronecker();
  }

} // end of namespace CoCoA


// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
  try
  {
    CoCoA::program();
    return 0;
  }
  catch (const CoCoA::ErrorInfo& err)
  {
    cerr << "***ERROR***  UNCAUGHT CoCoA Error";
    ANNOUNCE(cerr, err);
  }
  catch (const std::exception& exc)
  {
    cerr << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
  }
  catch(...)
  {
    cerr << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
  }

  CoCoA::BuildInfo::PrintAll(cerr);
  return 1;
}
