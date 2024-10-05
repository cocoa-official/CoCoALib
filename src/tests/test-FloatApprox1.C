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


#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRatOps.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/FloatApprox.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/error.H"

#include <iostream>
using std::cerr;
using std::endl;


namespace CoCoA
{

  // This test checks consistency between FloatApprox for ints and rationals.
  void program()
  {
    GlobalManager CoCoAFoundations;

    const int MantBits = 8;
    const int LogDenom = 9;
    const long denom = SmallPower(2,LogDenom);
    for (int i=1; i < 65536; ++i)
    {
      const MantExp2 integer = MantissaAndExponent2(i,MantBits);

      const BigRat q(i, denom);
      const MantExp2 rational = MantissaAndExponent2(q,MantBits);
      CoCoA_ASSERT_ALWAYS(integer.myMantissa == rational.myMantissa);
      CoCoA_ASSERT_ALWAYS(integer.myExponent == LogDenom+rational.myExponent);
      CoCoA_ASSERT_ALWAYS(integer.mySign == rational.mySign);
    }


    int prec = 8;
    for (int denom=1; denom < 100; ++denom)
      for (int numer=-1024; numer <= 1024; ++numer)
      {
        BigRat q(numer,denom);
        const BigRat approx = FloatApprox(q, prec);
        CoCoA_ASSERT_ALWAYS(abs(approx-q) <= abs(approx)/power(2,prec));
      }

    for (int i=1; i < 100; ++i)
    {
      const BigRat q(fibonacci(i*i+1), fibonacci(i));
      for (int prec=20; prec <= 160; prec*=2)
      {
        const BigRat approx = FloatApprox(q, prec);
        CoCoA_ASSERT_ALWAYS(abs(approx-q) <= abs(approx)/power(2,prec));
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
