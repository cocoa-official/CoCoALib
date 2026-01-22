//   Copyright (c)  2023,2025  John Abbott,  Anna M. Bigatti
//   Major contributions from Nico Mexis

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


#include "CoCoA/assert.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/combinatorics.H"
#include "CoCoA/error.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/SparsePolyOps-cyclotomic.H"

#include <algorithm>
using std::transform;
#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

namespace CoCoA
{

  void TestSpecialCases()
  {
    SparsePolyRing P = NewPolyRing(RingZZ(), symbols("x")); // ZZ[x];

    CoCoA_ASSERT_ALWAYS(CyclotomicTest(zero(P)) == 0);
    CoCoA_ASSERT_ALWAYS(CyclotomicIndex(zero(P)) == 0);

    CoCoA_ASSERT_ALWAYS(CyclotomicTest(one(P)) == 0);
    CoCoA_ASSERT_ALWAYS(CyclotomicIndex(one(P)) == 0);

    RingElem f = CyclotomicPoly(8, indet(P,0));
    CoCoA_ASSERT_ALWAYS(CyclotomicTest(f) == 8);
    CoCoA_ASSERT_ALWAYS(CyclotomicIndex(f) == 8);

    RingElem g = RingElem(P, "x^8-x^6+x^4-x^2+1");
    CoCoA_ASSERT_ALWAYS(CyclotomicTest(g) == 20);
    CoCoA_ASSERT_ALWAYS(CyclotomicIndex(g) == 20);

    RingElem h = RingElem(P, "x^19+x^10+1");
    CoCoA_ASSERT_ALWAYS(CyclotomicTest(h) == 0);
    CoCoA_ASSERT_ALWAYS(CyclotomicIndex(h) == 0);

    RingElem l = RingElem(P, "x^10+x^9-x^7-x^6-x^5-x^4-x^3+x+1");
    CoCoA_ASSERT_ALWAYS(CyclotomicTest(l) == 0);
    // Used to be a false positive!

    SparsePolyRing Pxy = NewPolyRing(RingZZ(), symbols("x,y")); // ZZ[x,y];
    RingElem m = RingElem(Pxy, "y^10+x^9-x^7-x^6-x^5-x^4-x^3+x+1");
    try
    {
      CyclotomicTest(m);
      CoCoA_ASSERT_ALWAYS(!"Non-univariate polynomials should throw in CyclotomicTest!");
    }
    catch (const CoCoA::ErrorInfo &err)
    {
      if (err != ERR::ReqUnivariate) throw;
    }
    // No corresponding multivariate check for CyclotomicIndex (since it may not "read" the whole poly!)
  }

  void TestUpto499()
  {
    const SparsePolyRing Px = NewPolyRing(RingZZ(), symbols("x"));
    const RingElem x = indet(Px, 0);

    for (unsigned long i = 1; i < 500; ++i)
    {
      const RingElem f = CyclotomicPoly(i, x);
      CoCoA_ASSERT_ALWAYS(CyclotomicTest(f) == i);
      CoCoA_ASSERT_ALWAYS(CyclotomicIndex(f) == i);
      CoCoA_ASSERT_ALWAYS(CyclotomicTest(2*f) == 0);
      CoCoA_ASSERT_ALWAYS(CyclotomicIndex(2*f) == 0);
      const long degf = deg(f);
      if (degf < 8) continue;
      RingElem g = f + (x*x-1)*(x*x-1)*power(x,degf/2-2);
      CoCoA_ASSERT_ALWAYS(i==20 || CyclotomicTest(g) == 0);
    }
  }

  void TestExbugs()
  {
    const SparsePolyRing Px = NewPolyRing(RingQQ(), symbols("x"));
    const RingElem x = indet(Px, 0);

    CoCoA_ASSERT_ALWAYS(CyclotomicFactorIndices(one(Px)).empty());
    CoCoA_ASSERT_ALWAYS(CyclotomicFactorIndices(2*one(Px)).empty());
    CoCoA_ASSERT_ALWAYS(CyclotomicFactorIndices(one(Px)/2).empty());

    RingElem NastyFactor1 = one(Px); // large fixed divisor
    RingElem NastyFactor2 = one(Px); // large fixed divisor & palindromic
    for (int i=0; i < 50; ++i)
    {
      NastyFactor1 *= (x*x-i*i);
      NastyFactor2 *= (x*x-i*i)*(i*i*x*x-1);
    }

    for (long i = 6; i < 100; ++i) // deliberately skip 1,2,3,4,5 (because 6 often appears as false positive)
    {
      const RingElem f = CyclotomicPoly(i, x);
      vector<long> CFI = CyclotomicFactorIndices(f);
      CoCoA_ASSERT_ALWAYS(len(CFI) == 1 && CFI[0] == i);
      CFI = CyclotomicFactorIndices(f*NastyFactor1);
      CoCoA_ASSERT_ALWAYS(len(CFI) >= 3 && CFI[0] == 1 && CFI[1] == 2 && CFI.back() == i);
      CFI = CyclotomicFactorIndices(f*NastyFactor2);
      CoCoA_ASSERT_ALWAYS(len(CFI) >= 3 && CFI[0] == 1 && CFI[1] == 2 && CFI.back() == i);
    }
  }

  void TestAllCombinations()
  {
    const SparsePolyRing Px = NewPolyRing(RingZZ(), symbols("x"));
    const RingElem x = indet(Px, 0);

    for (SubsetIter it(10); !IsEnded(it); ++it)
    {
      vector<long> expected(len(*it));
      std::transform((*it).begin(), (*it).end(), expected.begin(), [](int i){ return i+1;});
      RingElem f = one(Px);
      for (int j: expected)  f *= CyclotomicPoly(j,x);
      const vector<long> indices = CyclotomicFactorIndices(f);
      CoCoA_ASSERT_ALWAYS(indices == expected);
    }
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    TestSpecialCases();
    TestUpto499();
    TestExbugs();
    TestAllCombinations();
  }

} // end of namespace CoCoA

//----------------------------------------------------------------------
// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
  try
  {
    CoCoA::program();
    return 0;
  }
  catch (const CoCoA::ErrorInfo &err)
  {
    cerr << "***ERROR***  UNCAUGHT CoCoA Error";
    ANNOUNCE(cerr, err);
  }
  catch (const std::exception &exc)
  {
    cerr << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
  }
  catch (...)
  {
    cerr << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
  }

  CoCoA::BuildInfo::PrintAll(cerr);
  return 1;
}
