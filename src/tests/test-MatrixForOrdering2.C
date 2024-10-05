//   Copyright (c)  2012  John Abbott, Anna M Bigatti

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


#include "CoCoA/BigRat.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/MatrixForOrdering.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/error.H"
#include "CoCoA/VectorOps.H"
#include "CoCoA/matrix.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/utils.H"


#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

//----------------------------------------------------------------------
// Test for matrix fro orderings
// functions: ElimMat, ElimHomogMat
// environments: ZZ, QQ
//----------------------------------------------------------------------


namespace CoCoA
{

// convention: a function containing a "new" should be named "New.."
  matrix NewMatrixFromC(ring K, int* cmat, long NumRows, long NumCols)
  {
    matrix M(NewDenseMat(K,NumRows,NumCols));

    for (long i=0; i < NumRows; ++i)
      for (long j=0; j < NumCols; ++j)
        SetEntry(M, i, j, cmat[i*NumCols+j]);
    return M;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    const ring ZZ = RingZZ();
    const ring QQ = RingQQ();

    int M1[4*4] = {1, 0, 0, 4,
                   0, 3, 2, 0,
                   0,-3, 2, 1,
                   1, 0, 0, 0};

    const long N = 4;
    const matrix M = NewMatrixFromC(ZZ, M1, N,N);

    const std::vector<long> elim23 = LongRange(2,3);
  
    std::cout << ElimMat(elim23, N) << endl << endl;
    std::cout << ElimMat(LongRange(0,2), N) << endl << endl;
    std::cout << ElimMat(elim23, FirstRows(M, 0)) << endl << endl;
    std::cout << ElimMat(elim23, FirstRows(M, 2)) << endl << endl;
    std::cout << ElimHomogMat(elim23, FirstRows(M, 2)) << endl << endl;

    std::cout << MakeTermOrdMat(RowMat(M,0)) << endl << endl;
    std::cout << ElimHomogMat(LongRange(2,3), RowMat(M,0)) << endl << endl;
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
