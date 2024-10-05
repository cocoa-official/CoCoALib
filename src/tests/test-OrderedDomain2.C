//   Copyright (c)  2016  John Abbott,  Anna M. Bigatti

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


#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRatOps.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/RingTwinFloat.H"
#include "CoCoA/ring.H"
#include "CoCoA/utils.H"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <vector>
using std::vector;



namespace CoCoA
{

  void TestComparisons(const ring& R)
  {
    CoCoA_ASSERT_ALWAYS(IsOrderedDomain(R));
    const RingElem a = one(R);
    const RingElem b = -one(R);

    CoCoA_ASSERT_ALWAYS(cmp(a, b) > 0);
    CoCoA_ASSERT_ALWAYS(cmp(b, a) < 0);
    CoCoA_ASSERT_ALWAYS(cmp(a, a) == 0);
    CoCoA_ASSERT_ALWAYS(cmp(b, b) == 0);
    CoCoA_ASSERT_ALWAYS(cmp(a, b*b) == 0);

    CoCoA_ASSERT_ALWAYS(a > b);
    CoCoA_ASSERT_ALWAYS(a >= b);
    CoCoA_ASSERT_ALWAYS(b < a);
    CoCoA_ASSERT_ALWAYS(b <= a);

    CoCoA_ASSERT_ALWAYS(!(a > a));
    CoCoA_ASSERT_ALWAYS(a >= a);
    CoCoA_ASSERT_ALWAYS(!(a < a));
    CoCoA_ASSERT_ALWAYS(a <= a);

    CoCoA_ASSERT_ALWAYS(!(b > b));
    CoCoA_ASSERT_ALWAYS(b >= b);
    CoCoA_ASSERT_ALWAYS(!(b < b));
    CoCoA_ASSERT_ALWAYS(b <= b);

    CoCoA_ASSERT_ALWAYS(sign(a) == 1);
    CoCoA_ASSERT_ALWAYS(sign(b) == -1);
    CoCoA_ASSERT_ALWAYS(sign(a+b) == 0);
    CoCoA_ASSERT_ALWAYS(sign(a*b) == -1);
    CoCoA_ASSERT_ALWAYS(sign(b*b) == 1);
  }


// Test floor, ceil and NearestInteger
  void TestFloor(const ring& R)
  {
    CoCoA_ASSERT_ALWAYS(IsOrderedDomain(R));

    const BigRat half(1,2);

    vector<BigInt> BaseValues;
    BaseValues.push_back(BigInt(0));
    BaseValues.push_back(BigInt(1));
    BaseValues.push_back(BigInt(3));
    BaseValues.push_back(power(3,50));
    BaseValues.push_back(BigInt(-1));
    BaseValues.push_back(BigInt(-3));
    BaseValues.push_back(-power(3,50));

    for (int i=0; i < len(BaseValues); ++i)
    {
      const BigInt BaseValue = BaseValues[i];

      const RingElem exact(R, BaseValue);
      CoCoA_ASSERT_ALWAYS(floor(exact) == BaseValue);
      CoCoA_ASSERT_ALWAYS(ceil(exact) == BaseValue);
      CoCoA_ASSERT_ALWAYS(NearestInt(exact) == BaseValue);

      // Check for values of form integer+half:
      const RingElem UpHalf(R, BaseValue + half);
      const RingElem DownHalf(R, BaseValue - half);
      CoCoA_ASSERT_ALWAYS(floor(UpHalf) == BaseValue);
      CoCoA_ASSERT_ALWAYS(ceil(UpHalf) == BaseValue+1);
      CoCoA_ASSERT_ALWAYS(NearestInt(UpHalf) == BaseValue || NearestInt(UpHalf) == BaseValue+1);
      CoCoA_ASSERT_ALWAYS(floor(DownHalf) == BaseValue-1);
      CoCoA_ASSERT_ALWAYS(ceil(DownHalf) == BaseValue);
      CoCoA_ASSERT_ALWAYS(NearestInt(DownHalf) == BaseValue || NearestInt(DownHalf) == BaseValue-1);
      
      // Check for values of form integer +/- eps with 0 < eps < 1/2
      for (int j=2; j < 100; ++j)
      {
        const BigRat eps = power(half, j);
        const RingElem up(R, BaseValue + eps);
        const RingElem down(R, BaseValue - eps);
        CoCoA_ASSERT_ALWAYS(floor(up) == BaseValue);
        CoCoA_ASSERT_ALWAYS(ceil(up) == BaseValue+1);
        CoCoA_ASSERT_ALWAYS(NearestInt(up) == BaseValue);

        CoCoA_ASSERT_ALWAYS(floor(down) == BaseValue-1);
        CoCoA_ASSERT_ALWAYS(ceil(down) == BaseValue);
        CoCoA_ASSERT_ALWAYS(NearestInt(down) == BaseValue);
      }

    }
  }

  
  void CheckIsOrdered(const ring& R)
  {
    CoCoA_ASSERT_ALWAYS(IsOrderedDomain(R));
    TestComparisons(R);
    TestFloor(R);
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    // Low precision twin-floats will complain about InsufficientPrecision...
    bool threw = false;
    try
    {
      CheckIsOrdered(NewRingTwinFloat(64));
    }
    // We expect an InsufficientPrecision exception
    catch (const RingTwinFloat::InsufficientPrecision&) 
    {
      // Ignore expected exception.
      threw = true;
    }
    CoCoA_ASSERT_ALWAYS(threw);

    // Higher precision twin-floats do ot throw...
    CheckIsOrdered(NewRingTwinFloat(99));
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
