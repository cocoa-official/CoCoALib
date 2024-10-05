//   Copyright (c)  2010  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/error.H"

#include <cmath>
using std::abs;
#include <iostream>
using std::cerr;
using std::endl;


namespace CoCoA
{

  // test factorial, binomial, fibonacci, RoundDiv
  void program()
  {
    GlobalManager CoCoAFoundations(UseGMPAllocator); // for speed specify GMPAllocator (in case the default changes)

    constexpr int Nmax = 200;

    // factorial function
    CoCoA_ASSERT_ALWAYS(factorial(0) == 1);
    for (int i=1; i <= Nmax; ++i)
    {
      CoCoA_ASSERT_ALWAYS(factorial(i) == factorial(BigInt(i)));
      CoCoA_ASSERT_ALWAYS(i*factorial(i-1) == factorial(i));
    }

    // binomial function
    for (int i=-Nmax; i <= Nmax; ++i)
      for (int j=-Nmax-2; j <= Nmax+2; ++j)
      {
        if (j < 0 || j > std::abs(i))
        {
          CoCoA_ASSERT_ALWAYS(binomial(i,j) == 0);
          CoCoA_ASSERT_ALWAYS(binomial(BigInt(i),j) == 0);
          CoCoA_ASSERT_ALWAYS(binomial(i,BigInt(j)) == 0);
          CoCoA_ASSERT_ALWAYS(binomial(BigInt(i),BigInt(j)) == 0);
          continue;
        }
        if (j == 0)
        {
          CoCoA_ASSERT_ALWAYS(binomial(i,j) == 1);
          CoCoA_ASSERT_ALWAYS(binomial(BigInt(i),j) == 1);
          CoCoA_ASSERT_ALWAYS(binomial(i,BigInt(j)) == 1);
          CoCoA_ASSERT_ALWAYS(binomial(BigInt(i),BigInt(j)) == 1);
          continue;
        }

        // Now have j >= 1
        const BigInt b = binomial(i,j);
        CoCoA_ASSERT_ALWAYS(b == binomial(BigInt(i),j));
        CoCoA_ASSERT_ALWAYS(b == binomial(i,BigInt(j)));
        CoCoA_ASSERT_ALWAYS(b == binomial(BigInt(i),BigInt(j)));
        CoCoA_ASSERT_ALWAYS(b == binomial(i-1,j-1) + binomial(i-1,j));
      }

    // fibonacci
    CoCoA_ASSERT_ALWAYS(fibonacci(0) == 0);
    CoCoA_ASSERT_ALWAYS(fibonacci(1) == 1);
    for (int i=0; i <= Nmax; ++i)
    {
      const BigInt a = fibonacci(i);
      const BigInt b = fibonacci(i+1);
      const BigInt c = fibonacci(i+2);
      CoCoA_ASSERT_ALWAYS(a+b == c);
      CoCoA_ASSERT_ALWAYS(a == fibonacci(BigInt(i)));
    }

    // RoundDiv
    for (int i = -Nmax; i <= Nmax; ++i)
      for (int j = -Nmax; j <= Nmax; ++j)
      {
        if (j == 0) continue;
        const int q = RoundDiv(i,j);
        CoCoA_ASSERT_ALWAYS(q == RoundDiv(BigInt(i), j));
        CoCoA_ASSERT_ALWAYS(q == RoundDiv(i, BigInt(j)));
        CoCoA_ASSERT_ALWAYS(q == RoundDiv(BigInt(i), BigInt(j)));
      }
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
