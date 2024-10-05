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
#include "CoCoA/BigRatOps.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/NumTheory-gcd.H"
#include "CoCoA/error.H"

#include <iostream>
using std::cerr;
using std::endl;


namespace CoCoA
{

  // This test checks that the BigRat::AlreadyReduced flag works correctly.
  void program()
  {
    GlobalManager CoCoAFoundations;

    constexpr int NMAX = 99;
    for (int n=-NMAX; n <= NMAX; ++n)
      for (int d=-NMAX; d <= NMAX; ++d)
      {
        if (d == 0 || gcd(n,d) > 1) continue;
        const BigRat q = BigRat(n,d);
        CoCoA_ASSERT_ALWAYS(q == BigRat(n,d,                BigRat::AlreadyReduced));
        CoCoA_ASSERT_ALWAYS(q == BigRat(BigInt(n),d,        BigRat::AlreadyReduced));
        CoCoA_ASSERT_ALWAYS(q == BigRat(n,BigInt(d),        BigRat::AlreadyReduced));
        CoCoA_ASSERT_ALWAYS(q == BigRat(BigInt(n),BigInt(d),BigRat::AlreadyReduced));

        CoCoA_ASSERT_ALWAYS(den(BigRat(n,d,                BigRat::AlreadyReduced)) > 0);
        CoCoA_ASSERT_ALWAYS(den(BigRat(BigInt(n),d,        BigRat::AlreadyReduced)) > 0);
        CoCoA_ASSERT_ALWAYS(den(BigRat(n,BigInt(d),        BigRat::AlreadyReduced)) > 0);
        CoCoA_ASSERT_ALWAYS(den(BigRat(BigInt(n),BigInt(d),BigRat::AlreadyReduced)) > 0);

        if (n != 0)
        {
          const BigRat recip = BigRat(d,n,BigRat::AlreadyReduced);
          CoCoA_ASSERT_ALWAYS(power(q,-1) == recip);
          for (int e=1; e < 9; ++e)
            CoCoA_ASSERT_ALWAYS(power(q,-e) == power(recip,e));
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
