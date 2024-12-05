//   Copyright (c)  2011 Anna M. Bigatti

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
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/RingDistrMPolyClean.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/RingDenseUPolyClean.H"
#include "CoCoA/DenseUPolyRing.H"
#include "CoCoA/TmpToric.H"
#include "CoCoA/matrix.H"
#include "CoCoA/symbol.H"
#include "CoCoA/time.H"
#include "CoCoA/MatrixView.H"


#include <algorithm>
using std::min;
#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

//----------------------------------------------------------------------
// First test for Hilbert.
//----------------------------------------------------------------------

namespace CoCoA
{

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

    //  double t0 = CpuTime();

    SparsePolyRing P = NewPolyRing(NewZZmod(2), symbols("w,x,y,z"));
    RingElem w = RingElem(P, symbol("w"));
    RingElem x = RingElem(P, symbol("x"));
    RingElem y = RingElem(P, symbol("y"));
    RingElem z = RingElem(P, symbol("z"));

    std::vector<long> indices;
    indices.push_back(2);
  
    ideal I = ideal(x*z-y*y, x*w-y*z);
    CoCoA_ASSERT_ALWAYS(SequentialToric_C(I,indices) == ideal(y*y+x*z, w*x+y*z, w*y+z*z));
   
    int M0[1*3] = { 1,3,2 };
    matrix M0CC = NewMatrixFromC(RingQQ(),M0, 1,3);
    CoCoA_ASSERT_ALWAYS(SequentialToric_C(P, M0CC) == ideal(w*w +y, w*w*w +x));

    int M1[2*3] = { 1,3,2, 3,4,8 };
    matrix M1CC = NewMatrixFromC(RingQQ(),M1, 2,3);
    CoCoA_ASSERT_ALWAYS(SequentialToric_C(P, M1CC) == ideal(power(w,16) +x*x*power(y,5)));

    ideal T = SequentialToric_C(P, M1CC);
  
    int M2[5*4] = {0,1,0,0,  2,3,0,0,  3,1,0,0,  0,1,0,1,  1,0,1,0};
    matrix M2CC = NewMatrixFromC(RingQQ(),M2, 5,4);  
    CoCoA_ASSERT_ALWAYS(SequentialToric_C(P, M2CC) == ideal(w*0));

    int M3[5*4] = {0,1,0,0,  2,3,0,0,  3,1,0,0,  0,1,0,1,  1,0,1,0};
    matrix M3CC = NewMatrixFromC(RingQQ(),M3, 5,4);  
    CoCoA_ASSERT_ALWAYS(SequentialToric_C(P, M2CC) == ideal(w*0));

    EndToric_C();  // calls EndPoincare_C
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
    cerr << "***ERROR***  UNCAUGHT CoCoA error";
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
