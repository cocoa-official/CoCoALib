//   Copyright (c)  2014  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRatOps.H"
#include "CoCoA/NumTheory-gcd.H"
#include "CoCoA/utils.H"


#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;


namespace CoCoA
{

  // tests for CmpAbs of BigRat, and ILogBase of BigRat

  // This test checks that CmpAbs(q1,q2) gives the same answer as cmp(abs(q1),abs(q2))
  void TestCmpAbs()
  {
    constexpr int MAX = 40; // magic number -- takes a reasonable time (~2s) on my computer.
    vector<BigRat> q;
    q.push_back(BigRat(0,1));
    for (int n=1; n <= MAX; ++n)
      for (int d=1; d <= MAX; ++d)
        if (gcd(n,d)==1) q.push_back(BigRat(n,d));

    const int N = len(q);
    for (int i=0; i < N; ++i)
      for (int j=0; j < N; ++j)
      {
        const int ans = cmp(q[i],q[j]); // we know q[i] & q[j] are non-neg
        CoCoA_ASSERT_ALWAYS(ans == CmpAbs(q[i],q[j]));
        CoCoA_ASSERT_ALWAYS(ans == CmpAbs(-q[i],q[j]));
        CoCoA_ASSERT_ALWAYS(ans == CmpAbs(q[i],-q[j]));
        CoCoA_ASSERT_ALWAYS(ans == CmpAbs(-q[i],-q[j]));
      }
  }


  // Check that FloorLogBase gives plausible answers
  void TestFloorLogBase()
  {
    constexpr int MaxBase = 100;
    constexpr int MaxExp = 100;
    const BigRat SlightlyBigger(101,100);
    const BigRat SlightlySmaller(99,100);
    for (int base=2; base <= MaxBase; ++base)
      for (int exp=-MaxExp; exp <= MaxExp; ++exp)
      {
        BigRat b(base,1);
        CoCoA_ASSERT_ALWAYS(FloorLogBase(power(b,exp),base) == exp);
        CoCoA_ASSERT_ALWAYS(FloorLogBase(SlightlyBigger*power(b,exp),base) == exp);
        CoCoA_ASSERT_ALWAYS(FloorLogBase(SlightlySmaller*power(b,exp),base) == exp-1);

        if (b == 2)
        {
          CoCoA_ASSERT_ALWAYS(FloorLog2(power(b,exp)) == exp);
          CoCoA_ASSERT_ALWAYS(FloorLog2(SlightlyBigger*power(b,exp)) == exp);
          CoCoA_ASSERT_ALWAYS(FloorLog2(SlightlySmaller*power(b,exp)) == exp-1);
        }

        if (b == 10)
        {
          CoCoA_ASSERT_ALWAYS(FloorLog10(power(b,exp)) == exp);
          CoCoA_ASSERT_ALWAYS(FloorLog10(SlightlyBigger*power(b,exp)) == exp);
          CoCoA_ASSERT_ALWAYS(FloorLog10(SlightlySmaller*power(b,exp)) == exp-1);
        }
      }
  }


  // Check that mantissa works plausibly
  void TestMantissa()
  {
    const BigRat q(37*power(13,100)-1,power(65536,2)*power(13,101));
    long e0;
    const double m0 = mantissa(e0, q);
    for (int j=0; j < 64; ++j)
    {
      long e;
      double m = mantissa(e, q*power(2,j));
      CoCoA_ASSERT_ALWAYS(m == m0);
      CoCoA_ASSERT_ALWAYS(e == e0+j);
    }
  }


  void program()
  {
    GlobalManager CoCoAFoundations;
    TestCmpAbs();
    TestFloorLogBase();
    TestMantissa();
  }

}

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
