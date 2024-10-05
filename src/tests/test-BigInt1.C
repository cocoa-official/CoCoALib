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

#include <iostream>
using std::cerr;
using std::endl;


namespace CoCoA
{

  // test BigInt ctor, equality & not-equality, cmp, basic arith
  void program()
  {
    GlobalManager CoCoAFoundations(UseGMPAllocator); // for speed specify GMPAllocator (in case the default changes)

    {
      // test ctors
      const BigInt N1;
      CoCoA_ASSERT_ALWAYS(IsZero(N1));

      const BigInt N2 = BigIntFromString("0");
      CoCoA_ASSERT_ALWAYS(N1 == N2);
      CoCoA_ASSERT_ALWAYS(N2 == N1);

      CoCoA_ASSERT_ALWAYS(IsOne(BigInt(1)));
      CoCoA_ASSERT_ALWAYS(IsMinusOne(BigInt(-1)));
    }

    constexpr int LIMIT = 100;

    // Check some equalities/inequalities
    for (int i = -LIMIT; i <= LIMIT; ++i)
      for (int j = -LIMIT; j <= LIMIT; ++j)
      {
        if (i == j)
        {
          CoCoA_ASSERT_ALWAYS(i == BigInt(j));
          CoCoA_ASSERT_ALWAYS(BigInt(i) == j);
          CoCoA_ASSERT_ALWAYS(BigInt(i) == BigInt(j));
        }
        if (i != j)
        {
          CoCoA_ASSERT_ALWAYS(i != BigInt(j));
          CoCoA_ASSERT_ALWAYS(BigInt(i) != j);
          CoCoA_ASSERT_ALWAYS(BigInt(i) != BigInt(j));
        }
        if (i < j)
        {
          CoCoA_ASSERT_ALWAYS(i < BigInt(j));
          CoCoA_ASSERT_ALWAYS(BigInt(i) < j);
          CoCoA_ASSERT_ALWAYS(BigInt(i) < BigInt(j));
        }
        if (i <= j)
        {
          CoCoA_ASSERT_ALWAYS(i <= BigInt(j));
          CoCoA_ASSERT_ALWAYS(BigInt(i) <= j);
          CoCoA_ASSERT_ALWAYS(BigInt(i) <= BigInt(j));
        }
        if (i > j)
        {
          CoCoA_ASSERT_ALWAYS(i > BigInt(j));
          CoCoA_ASSERT_ALWAYS(BigInt(i) > j);
          CoCoA_ASSERT_ALWAYS(BigInt(i) > BigInt(j));
        }
        if (i >= j)
        {
          CoCoA_ASSERT_ALWAYS(i >= BigInt(j));
          CoCoA_ASSERT_ALWAYS(BigInt(i) >= j);
          CoCoA_ASSERT_ALWAYS(BigInt(i) >= BigInt(j));
        }
      }

    // Check some cmp function calls.
    for (int i = -LIMIT; i <= LIMIT; ++i)
      for (int j = -LIMIT; j <= LIMIT; ++j)
      {
        if (cmp(i,j) > 0)
        {
          CoCoA_ASSERT_ALWAYS(cmp(i, BigInt(j)) > 0);
          CoCoA_ASSERT_ALWAYS(cmp(BigInt(i), j) > 0);
          CoCoA_ASSERT_ALWAYS(cmp(BigInt(i), BigInt(j)) > 0);
        }
        if (cmp(i,j) == 0)
        {
          CoCoA_ASSERT_ALWAYS(cmp(i, BigInt(j)) == 0);
          CoCoA_ASSERT_ALWAYS(cmp(BigInt(i), j) == 0);
          CoCoA_ASSERT_ALWAYS(cmp(BigInt(i), BigInt(j)) == 0);
        }
        if (cmp(i,j) < 0)
        {
          CoCoA_ASSERT_ALWAYS(cmp(i, BigInt(j)) < 0);
          CoCoA_ASSERT_ALWAYS(cmp(BigInt(i), j) < 0);
          CoCoA_ASSERT_ALWAYS(cmp(BigInt(i), BigInt(j)) < 0);
        }
      }


    // Check some basic arithmetic operations.
    for (int i = -LIMIT; i <= LIMIT; ++i)
      for (int j = -LIMIT; j <= LIMIT; ++j)
      {
        CoCoA_ASSERT_ALWAYS(i+j == i+BigInt(j));
        CoCoA_ASSERT_ALWAYS(i+j == BigInt(i)+j);
        CoCoA_ASSERT_ALWAYS(i+j == BigInt(i)+BigInt(j));

        CoCoA_ASSERT_ALWAYS(i-j == i-BigInt(j));
        CoCoA_ASSERT_ALWAYS(i-j == BigInt(i)-j);
        CoCoA_ASSERT_ALWAYS(i-j == BigInt(i)-BigInt(j));

        CoCoA_ASSERT_ALWAYS(i*j == i*BigInt(j));
        CoCoA_ASSERT_ALWAYS(i*j == BigInt(i)*j);
        CoCoA_ASSERT_ALWAYS(i*j == BigInt(i)*BigInt(j));

        // Restrict to non-negative for division as C++ is not well-defined with negative numbers.
        if (i >= 0 && j > 0)
        {
          CoCoA_ASSERT_ALWAYS(i/j == i/BigInt(j));
          CoCoA_ASSERT_ALWAYS(i/j == BigInt(i)/j);
          CoCoA_ASSERT_ALWAYS(i/j == BigInt(i)/BigInt(j));
        }
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
