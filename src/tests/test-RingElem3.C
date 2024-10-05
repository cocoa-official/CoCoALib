//   Copyright (c)  2013  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/error.H"

#include <iostream>
using std::cerr;
using std::endl;


namespace CoCoA
{

  void TestIsDivisible()
  {
    ring ZZmod6 = NewZZmod(6);
    const RingElem a(ZZmod6, 4);
    const RingElem b(ZZmod6, 2);
    CoCoA_ASSERT_ALWAYS(!IsDivisible(a,b));
    CoCoA_ASSERT_ALWAYS(!IsDivisible(a,a));
    CoCoA_ASSERT_ALWAYS(!IsDivisible(b,b));
    CoCoA_ASSERT_ALWAYS(!IsDivisible(0,a));
    CoCoA_ASSERT_ALWAYS(!IsDivisible(0,b));
    CoCoA_ASSERT_ALWAYS(!IsDivisible(0,zero(ZZmod6)));
  }

  void TestZeroPowerZero(ring R)
  {
    const RingElem z = zero(R);
    CoCoA_ASSERT_ALWAYS(IsOne(power(z,0)));
    CoCoA_ASSERT_ALWAYS(IsOne(power(z,BigInt(0))));
  }

  void program()
  {
    GlobalManager CoCoAFoundations;
    TestIsDivisible();
    TestZeroPowerZero(RingZZ());
    TestZeroPowerZero(RingQQ());
    TestZeroPowerZero(NewZZmod(3));
    TestZeroPowerZero(NewZZmod(4));
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
