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
#include "CoCoA/NumTheory-prime.H"
#include "CoCoA/NumTheory-modular.H"
#include "CoCoA/NumTheory-gcd.H"
#include "CoCoA/error.H"
#include "CoCoA/time.H"

#include <iostream>
using std::cerr;
using std::endl;


// Tests PrimitiveRoot, MultiplicativeOrderMod; assumes PowerMod works OK.

namespace CoCoA
{

  void TestPrime(long p)
  {
    CoCoA_ASSERT_ALWAYS(IsPrime(p));

    const long g = PrimitiveRoot(p);
    const long order = p-1;
    for (int i=1; i < p; ++i)
    {
      const long x = PowerMod(g,i,p);
      CoCoA_ASSERT_ALWAYS(MultiplicativeOrderMod(x,p) == order/long(gcd(order,i)));
    }
  }

  void TestPrime(const BigInt& P)
  {
    CoCoA_ASSERT_ALWAYS(IsProbPrime(P));

    const double StartTime = CpuTime();
    const long g = PrimitiveRoot(P);
    const BigInt order(P-1);
    for (int i=1; i <= 9; ++i)
    {
      if (CpuTime() > StartTime+5) break;
      const BigInt x = PowerMod(g,i,P);
      CoCoA_ASSERT_ALWAYS(MultiplicativeOrderMod(x,P) == order/gcd(order,i));
    }
  }

  void program()
  {
    GlobalManager CoCoAFoundations;

    int p = 2;
    while (p < 2000)
    {
      TestPrime(p);
      p = NextPrime(p);
    }


    BigInt P(2);
    while (P < 2000)
    {
      TestPrime(P);
      P = NextProbPrime(P);
    }

    TestPrime(NextProbPrime(power(2,32)));
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
