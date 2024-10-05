//   Copyright (c)  2022  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/BigIntOps.H"
#include "CoCoA/combinatorics.H"
#include "CoCoA/utils.H"

#include <iostream>
using std::cerr;
using std::endl;

//#include <vector>
using std::vector;

//----------------------------------------------------------------------
// This test does virtually nothing, but is a handy template if you want
// to create your own test code: just copy this file and add your code
// after the line with "PUT YOUR CODE HERE"
//----------------------------------------------------------------------

namespace CoCoA
{
  // Put your code inside namespace CoCoA to avoid possibile
  // ambiguities with STL fns sharing names with CoCoALib fns.

  void program()
  {
    // you may use CoCoALib functions AFTER creating a GlobalManager:
    GlobalManager CoCoAFoundations;

    CoCoA_ASSERT_ALWAYS(NumPartitions(9) == 30);
    CoCoA_ASSERT_ALWAYS(CatalanNumber(9) == 4862);

    const vector<long> pi = RandomPermutation(9);
    const int sigma = signature(pi);
    CoCoA_ASSERT_ALWAYS(sigma == 1 || sigma == -1);

    vector<int> pi2(9);
    for (int i=0; i < 9; ++i) pi2[i] = pi[i];
    const int sigma2 = signature(pi2);
    CoCoA_ASSERT_ALWAYS(sigma == sigma2);

    vector<long> S = RandomSubsetIndices(9);
    CoCoA_ASSERT_ALWAYS(len(S) >= 0 && len(S) <= 9);

    S = RandomSubsetIndices(9,3);
    CoCoA_ASSERT_ALWAYS(len(S) == 3);

    SubsetIter it(3);
    CoCoA_ASSERT_ALWAYS(!IsEnded(it));
    CoCoA_ASSERT_ALWAYS(len(*it) == 0);
    ++it;
    CoCoA_ASSERT_ALWAYS(!IsEnded(it));
    CoCoA_ASSERT_ALWAYS(len(*it) == 1);
    for (int i=1; i <= 7; ++i) ++it;
    CoCoA_ASSERT_ALWAYS(IsEnded(it));
    CoCoA_ASSERT_ALWAYS(len(*it) == 0);

    SubsetIter it0(0);
    CoCoA_ASSERT_ALWAYS(!IsEnded(it0));
    CoCoA_ASSERT_ALWAYS(len(*it0) == 0);
    ++it0;
    CoCoA_ASSERT_ALWAYS(IsEnded(it0));
    CoCoA_ASSERT_ALWAYS(len(*it0) == 0);
    
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
