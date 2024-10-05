//   Copyright (c)  2007,2014  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingTwinFloat.H"
#include "CoCoA/RingZZ.H"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

namespace CoCoA
{

  // test IsInteger for various values in various rings
  // (incl. an undecidable case in a RingTwinFloat)

  void TestInteger(const RingElem& x, bool ItIsInteger)
  {
    BigInt N;
    try
    {
      if (IsInteger(N, x))
        CoCoA_ASSERT_ALWAYS(ItIsInteger && x == N);
     else
        CoCoA_ASSERT_ALWAYS(!ItIsInteger);
    }
    catch (const RingTwinFloat::InsufficientPrecision&)
    {
      cout << "Failed to decide whether the ring element " << x
           << " is or is not the image of an integer." << endl;
    }

  }


  void trial(const ring& R)
  {
    const RingElem x(R, 29);
    const bool PositiveChar = (characteristic(R) > 0);
    TestInteger(x, true);
    TestInteger(power(x, 10), true);
    if (IsDivisible_AllowFields(x, 2)) TestInteger(x/2, PositiveChar);
    if (IsDivisible_AllowFields(x, 3)) TestInteger(x/3, PositiveChar);
    if (IsDivisible_AllowFields(x, x+1)) TestInteger(x/(x+1), PositiveChar);
  }

  void program()
  {
    GlobalManager CoCoAFoundations(UseNonNegResidues);

    trial(RingZZ());
    trial(RingQQ());
    trial(NewZZmod(2));
    trial(NewZZmod(3));
    trial(NewZZmod(32003));   // large prime
    trial(NewZZmod(1000003)); // larger prime
    trial(NewZZmod(6*29)); // ring has zero divisors
    trial(NewZZmod(1048576)); // ring has zero divisors
    trial(NewRingTwinFloat(32));
    trial(NewRingTwinFloat(16)); // this will cause a "Failed..." message to be printed for 29^10
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
