//   Copyright (c)  2018  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/error.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/SparsePolyOps-SturmSeq.H"

#include <iostream>
using std::cerr;
using std::endl;

//----------------------------------------------------------------------
// This test checks for a bug that was in the Sturm seq code
// (the signs in the remainder sequence were incorrect).
//----------------------------------------------------------------------

namespace CoCoA
{
  // Put your code inside namespace CoCoA to avoid possibile
  // ambiguities with STL fns sharing names with CoCoALib fns.

  void program()
  {
    GlobalManager CoCoAFoundations;

    ring P = NewPolyRing(RingQQ(), symbols("x"));
    RingElem f1(P, "(x^2+1)^2+1");
    const long NumRoots1 = NumRealRoots(f1);
    CoCoA_ASSERT_ALWAYS(NumRoots1 == 0);

    RingElem f2(P, "(4*x^3-6*x+3)^2+1");
    const long NumRoots2 = NumRealRoots(f2);
    CoCoA_ASSERT_ALWAYS(NumRoots2 == 0);

    RingElem f3(P, "(3*x^4+2*x)^2+1");
    const long NumRoots3 = NumRealRoots(f3);
    CoCoA_ASSERT_ALWAYS(NumRoots3 == 0);

    RingElem f4(P, "(3*x^3-3*x^2+3*x)^2+1");
    const long NumRoots4 = NumRealRoots(f4);
    CoCoA_ASSERT_ALWAYS(NumRoots4 == 0);

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
