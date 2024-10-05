//   Copyright (c)  2007  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/FreeModule.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/RingTwinFloat.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/submodule.H"
#include "CoCoA/utils.H"

#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

namespace CoCoA
{

// This trial is still very rudimentary..
  void trial(ring R)
  {
    const int dim = 4;
    const FreeModule F = NewFreeModule(R, dim);
    const vector<ModuleElem>& Fgens = gens(F);
    CoCoA_ASSERT_ALWAYS(len(Fgens) == dim);
// This does not yet compile
//   for (int i=0; i < dim; ++i)
//     for (int j=0; j < dim; ++j)
//       CoCoA_ASSERT_ALWAYS(Fgens[i][j] == (i==j));
    CoCoA_ASSERT_ALWAYS(IsFGModule(F));

    const ModuleElem& e0 = Fgens[0];
    const ModuleElem& e1 = Fgens[1];
    const ModuleElem& e2 = Fgens[2];
    const ModuleElem& e3 = Fgens[3];

    const FGModule M = submodule(1*e0 +2*e1 +3*e2 +4*e3, 4*e0 +3*e1 +2*e2 +1*e3);
    vector<ModuleElem> Mgens = gens(M);

    CoCoA_ASSERT_ALWAYS(gens(M).size() == 2);
    CoCoA_ASSERT_ALWAYS(e2 != e3);
    CoCoA_ASSERT_ALWAYS(Mgens[0]+Mgens[1] == 5*(e0 + e1 + e2 + e3));
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    trial(RingZZ());
    trial(NewRingTwinFloat(32));
    trial(NewZZmod(2));
    trial(NewZZmod(1048576)); // ring has zero divisors
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
