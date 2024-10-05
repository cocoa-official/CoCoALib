//   Copyright (c)  2012  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/convert.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/BigRatOps.H"
#include "CoCoA/NumTheory-SimplestRat.H"
#include "CoCoA/NumTheory-ContFrac.H"
#include "CoCoA/error.H"


#include <iostream>
using std::cerr;
using std::cout;
using std::endl;


// Tests CFApproximantsIter, CFApprox, ContFracIter &
//       SimplestBigRatBetween, SimplestBinaryRatBetween

namespace CoCoA
{

  void program()
  {
    // Test for functions related to continued fractions.
    GlobalManager CoCoAFoundations;

    // Generator/iterator for continued fraction "quotients"
    BigRat x(-101,100);
    while (x < 0)
    {
      ContFracIter cfi(x);
      cout << "The continued fraction quotients of " << x << " are:";
      while (!IsEnded(cfi))
        cout << "  " << *cfi++;
      cout << endl;
      x += BigRat(1,10);
    }

    // Make a rational from its cont frac quotients...
    ContFracApproximant approx;
    for (ContFracIter it(x); !IsEnded(it); ++it)
      approx.myAppendQuot(quot(it));
    CoCoA_ASSERT_ALWAYS(approx.myRational() == x);

    // Compute the approximants for a rational...
    x = BigRat(-123,100);
    for (CFApproximantsIter it(x); !IsEnded(it); ++it)
      cout << "Approximant: " << *it << endl;

    // CFApprox...
    CoCoA_ASSERT_ALWAYS(CFApprox(ConvertTo<BigRat>(1.4142135), BigRat(1)) == BigRat(1));
    CoCoA_ASSERT_ALWAYS(CFApprox(ConvertTo<BigRat>(-1.4142135), BigRat(1)) == BigRat(-2));
    CoCoA_ASSERT_ALWAYS(CFApprox(ConvertTo<BigRat>(1.4142135), ConvertTo<BigRat>(0.01)) == BigRat(17,12));
    CoCoA_ASSERT_ALWAYS(CFApprox(ConvertTo<BigRat>(-1.4142135), ConvertTo<BigRat>(0.01)) == BigRat(-17,12));

    CoCoA_ASSERT_ALWAYS(SimplestBigRatBetween(BigRat(1), BigRat(2)) == BigRat(1));
    CoCoA_ASSERT_ALWAYS(SimplestBigRatBetween(BigRat(-2), BigRat(-1)) == BigRat(-1));
    CoCoA_ASSERT_ALWAYS(SimplestBigRatBetween(BigRat(-1), BigRat(1)) == BigRat(0));
    CoCoA_ASSERT_ALWAYS(SimplestBigRatBetween(BigRat(0), BigRat(1)) == BigRat(0));
    CoCoA_ASSERT_ALWAYS(SimplestBigRatBetween(BigRat(-1), BigRat(0)) == BigRat(0));
    CoCoA_ASSERT_ALWAYS(SimplestBigRatBetween(BigRat(1), BigRat(4,3)) == BigRat(1));
    CoCoA_ASSERT_ALWAYS(SimplestBigRatBetween(BigRat(4,3), BigRat(2)) == BigRat(2));
    CoCoA_ASSERT_ALWAYS(SimplestBigRatBetween(BigRat(-4,3), BigRat(-1)) == BigRat(-1));
    CoCoA_ASSERT_ALWAYS(SimplestBigRatBetween(BigRat(-2), BigRat(-4,3)) == BigRat(-2));
    CoCoA_ASSERT_ALWAYS(SimplestBigRatBetween(BigRat(4,3), BigRat(10,3)) == BigRat(2));
    CoCoA_ASSERT_ALWAYS(SimplestBigRatBetween(BigRat(-10,3), BigRat(-4,3)) == BigRat(-2));
    CoCoA_ASSERT_ALWAYS(SimplestBigRatBetween(BigRat(123,100), BigRat(135,100)) == BigRat(4,3));
    CoCoA_ASSERT_ALWAYS(SimplestBigRatBetween(BigRat(-135,100), BigRat(-123,100)) == BigRat(-4,3));
    CoCoA_ASSERT_ALWAYS(SimplestBigRatBetween(BigRat(5,2), BigRat(299,100)) == BigRat(5,2));
    CoCoA_ASSERT_ALWAYS(SimplestBigRatBetween(BigRat(-299,100), BigRat(-5,2)) == BigRat(-5,2));
    CoCoA_ASSERT_ALWAYS(SimplestBigRatBetween(BigRat(4181, 2584), BigRat(6765,4181)) == BigRat(4181,2584));
    CoCoA_ASSERT_ALWAYS(SimplestBigRatBetween(BigRat(-6765,4181), BigRat(-4181, 2584)) == BigRat(-4181,2584));
    CoCoA_ASSERT_ALWAYS(SimplestBigRatBetween(BigRat(225,157), BigRat(268,187)) == BigRat(225,157));
    CoCoA_ASSERT_ALWAYS(SimplestBigRatBetween(BigRat(-268,187), BigRat(-225,157)) == BigRat(-225,157));
    CoCoA_ASSERT_ALWAYS(SimplestBigRatBetween(BigRat(1393,972), BigRat(268,187)) == BigRat(268,187));
    CoCoA_ASSERT_ALWAYS(SimplestBigRatBetween(BigRat(-268,187), BigRat(-1393,972)) == BigRat(-268,187));

    CoCoA_ASSERT_ALWAYS(SimplestBinaryRatBetween(BigRat(-1), BigRat(1)) == BigRat(0));
    CoCoA_ASSERT_ALWAYS(SimplestBinaryRatBetween(BigRat(-1), BigRat(2)) == BigRat(0));
    CoCoA_ASSERT_ALWAYS(SimplestBinaryRatBetween(BigRat(-2), BigRat(1)) == BigRat(0));
    CoCoA_ASSERT_ALWAYS(SimplestBinaryRatBetween(BigRat(-1), BigRat(3)) == BigRat(0));
    CoCoA_ASSERT_ALWAYS(SimplestBinaryRatBetween(BigRat(-2), BigRat(2)) == BigRat(0));
    CoCoA_ASSERT_ALWAYS(SimplestBinaryRatBetween(BigRat(-3), BigRat(1)) == BigRat(0));
    CoCoA_ASSERT_ALWAYS(SimplestBinaryRatBetween(BigRat(5), BigRat(10)) == BigRat(8));
    CoCoA_ASSERT_ALWAYS(SimplestBinaryRatBetween(BigRat(-10), BigRat(-5)) == BigRat(-8));
    CoCoA_ASSERT_ALWAYS(SimplestBinaryRatBetween(BigRat(9,10), BigRat(11,10)) == BigRat(1));
    CoCoA_ASSERT_ALWAYS(SimplestBinaryRatBetween(BigRat(-11,10), BigRat(-9,10)) == BigRat(-1));
    CoCoA_ASSERT_ALWAYS(SimplestBinaryRatBetween(BigRat(1,5), BigRat(1,3)) == BigRat(1,4));
    CoCoA_ASSERT_ALWAYS(SimplestBinaryRatBetween(BigRat(1,5), BigRat(1,4)) == BigRat(1,4));
    CoCoA_ASSERT_ALWAYS(SimplestBinaryRatBetween(BigRat(1,4), BigRat(1,3)) == BigRat(1,4));
    CoCoA_ASSERT_ALWAYS(SimplestBinaryRatBetween(BigRat(1,3), BigRat(1,2)) == BigRat(1,2));
    CoCoA_ASSERT_ALWAYS(SimplestBinaryRatBetween(BigRat(1,4), BigRat(1,2)) == BigRat(1,2));
    CoCoA_ASSERT_ALWAYS(SimplestBinaryRatBetween(BigRat(1,4), BigRat(1,2)) == BigRat(1,2));
    CoCoA_ASSERT_ALWAYS(SimplestBinaryRatBetween(BigRat(2,3), BigRat(5,6)) == BigRat(3,4));
    CoCoA_ASSERT_ALWAYS(SimplestBinaryRatBetween(BigRat(3,10), BigRat(4,10)) == BigRat(3,8));
    CoCoA_ASSERT_ALWAYS(SimplestBinaryRatBetween(BigRat(15), BigRat(20)) == BigRat(16));
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
