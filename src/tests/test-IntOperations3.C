//   Copyright (c)  2017  John Abbott,  Anna M. Bigatti

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

//----------------------------------------------------------------------
// This test checks that MultiplicityOf2 gives correct answer.
//----------------------------------------------------------------------

namespace CoCoA
{
  // Put your code inside namespace CoCoA to avoid possibile
  // ambiguities with STL fns sharing names with CoCoALib fns.

  void check(const BigInt& a, const BigInt& b)
  {
    const long multa = MultiplicityOf2(a);
    const long multb = MultiplicityOf2(b);
    CoCoA_ASSERT_ALWAYS(MultiplicityOf2(-a) == multa);
    CoCoA_ASSERT_ALWAYS(MultiplicityOf2(-b) == multb);
    const BigInt OddParta = a/power(2, multa);
    const BigInt OddPartb = b/power(2, multb);
    CoCoA_ASSERT_ALWAYS(IsOdd(OddParta));
    CoCoA_ASSERT_ALWAYS(a == OddParta*power(2, multa));
    CoCoA_ASSERT_ALWAYS(IsOdd(OddPartb));
    CoCoA_ASSERT_ALWAYS(b == OddPartb*power(2, multb));

    CoCoA_ASSERT_ALWAYS(MultiplicityOf2(a*b) == multa+multb);
  }

  void program()
  {
    // you may use CoCoALib functions AFTER creating a GlobalManager:
    GlobalManager CoCoAFoundations;

    for (int i=1; i < 1000; ++i)
      check(fibonacci(i), factorial(i));
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
