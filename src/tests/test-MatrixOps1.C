//   Copyright (c)  2021  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/matrix.H"
#include "CoCoA/MatrixOps.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/RingZZ.H"

#include <iostream>
using std::cerr;
using std::endl;

//----------------------------------------------------------------------
// This test does virtually nothing, but is a handy template if you want
// to create your own test code: just copy this file and add your code
// after the line with "PUT YOUR CODE HERE"
//----------------------------------------------------------------------

namespace CoCoA
{
  // Put your code inside namespace CoCoA to avoid possibile
  // ambiguities with STL fns sharing names with CoCoALib fns.

  void redmine1590()
  {
    const ring ZZ = RingZZ();
    matrix M = NewDenseMat(ZZ,3,2);
    SetEntry(M,0,0, 1);
    SetEntry(M,0,1, 1);
    SetEntry(M,1,0, 2);
    SetEntry(M,1,1, 0);
    SetEntry(M,2,0, 0);
    SetEntry(M,2,1, 3);
    const matrix K = LinKerZZ(M); // used to crash
    CoCoA_ASSERT_ALWAYS(NumRows(K)==0  && NumCols(K)==0);
  }

  void program()
  {
    // you may use CoCoALib functions AFTER creating a GlobalManager:
    GlobalManager CoCoAFoundations;

    // *** PUT YOUR CODE HERE ***
    CoCoA_ASSERT_ALWAYS(CoCoA::IsInitialized());
    redmine1590();
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
