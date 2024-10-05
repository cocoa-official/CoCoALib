//   Copyright (c)  2016  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/MatrixOps.H"
#include "CoCoA/MatrixSpecial.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/error.H"
#include "CoCoA/matrix.H"


#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

//----------------------------------------------------------------------
// Test for special matrices:
// functions: RandomUnimodularMat
//----------------------------------------------------------------------
namespace CoCoA
{


  void TestRandomUnimodularMat(const ring& R)
  {
    for (int i=0; i < 10; ++i)
      for (int niters = 0; niters < 10; ++niters)
      {
        const matrix M = RandomUnimodularMat(R,i,niters);
        CoCoA_ASSERT_ALWAYS(NumRows(M) == i);
        CoCoA_ASSERT_ALWAYS(NumCols(M) == i);
        CoCoA_ASSERT_ALWAYS(det(M) == 1 || det(M) == -1); // det is +1 or -1
      }
  }

void program()
{
  GlobalManager CoCoAFoundations;

  TestRandomUnimodularMat(RingZZ());
  TestRandomUnimodularMat(RingQQ());
  TestRandomUnimodularMat(NewZZmod(2));
  TestRandomUnimodularMat(NewZZmod(3));
  TestRandomUnimodularMat(NewZZmod(32003));
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
