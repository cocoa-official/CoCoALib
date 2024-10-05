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


#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/error.H"
#include "CoCoA/ring.H"

#include <iostream>
using std::cerr;
using std::endl;



namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    // Test default ring, and automatic mapping into ring when assigning.
    ring R;  // default is ZZ
    CoCoA_ASSERT_ALWAYS(R == RingZZ());
    RingElem x; // default is 0 in ZZ
    CoCoA_ASSERT_ALWAYS(owner(x) == RingZZ());
    x = 1;  // will not change ring
    CoCoA_ASSERT_ALWAYS(owner(x) == RingZZ());
    x = power(10,100);  // will not change ring
    CoCoA_ASSERT_ALWAYS(owner(x) == RingZZ());
    x = BigRat(2,1);    // will not change ring
    CoCoA_ASSERT_ALWAYS(owner(x) == RingZZ());

    // We expect error if rational is not an integer
    try { x = BigRat(2,3); CoCoA_ASSERT_ALWAYS(!"NEVER GET HERE!"); }
    catch (const ErrorInfo& err) { CoCoA_ASSERT_ALWAYS(err == ERR::EmbedBigRatFailed); }

    R = NewFractionField(R);  // now R is QQ
    CoCoA_ASSERT_ALWAYS(R == RingQQ());
    CoCoA_ASSERT_ALWAYS(owner(x) != R);

    RingElem y(R); // y belongs to QQ
    CoCoA_ASSERT_ALWAYS(owner(y) == RingQQ());
    y = BigRat(2,3);   // works because y is in QQ
    CoCoA_ASSERT_ALWAYS(owner(y) == RingQQ());

    CoCoA_ASSERT_ALWAYS(owner(x) != owner(y));
    swap(x, y);       // now x is in QQ, and y in ZZ
    CoCoA_ASSERT_ALWAYS(owner(x) == RingQQ());
    CoCoA_ASSERT_ALWAYS(owner(y) == RingZZ());

    x = y;             // x is now back in ZZ
    CoCoA_ASSERT_ALWAYS(owner(x) == owner(y));
    CoCoA_ASSERT_ALWAYS(owner(x) == RingZZ());
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
