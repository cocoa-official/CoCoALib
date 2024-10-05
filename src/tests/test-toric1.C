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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/tests/test-toric1.C,v 1.10 2016/09/16 16:36:42 abbott Exp $
// $Log: test-toric1.C,v $
// Revision 1.10  2016/09/16 16:36:42  abbott
// Summary: Changed TEST_ASSERT into CoCoA_ASSERT_ALWAYS; removed all include assert.H lines
//
// Revision 1.9  2015/07/27 12:45:47  bigatti
// -- now using "symbols(string)" comma separated symbols
//
// Revision 1.8  2015/05/11 13:00:53  abbott
// Summary: Put all code into namespace CoCoA (redmine #642)
// Author: JAA
//
// Revision 1.7  2014/04/17 13:47:21  bigatti
// -- MatrixViews --> MatrixView
//
// Revision 1.6  2013/02/01 17:31:27  bigatti
// -- uncommented tests
//
// Revision 1.5  2013/01/31 12:55:48  bigatti
// -- cleaner structure
//
// Revision 1.4  2013/01/25 14:53:15  bigatti
// -- added ex-bug
//
// Revision 1.3  2012/02/10 11:57:12  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.2  2012/02/08 17:38:15  bigatti
// -- changed Z,Q --> ZZ,QQ
//
// Revision 1.1  2011/05/09 14:50:06  bigatti
// -- added test-toric1.C
//
// Revision 1.15  2011/04/27 09:48:37  bigatti
// -- updated syntax
//
// Revision 1.14  2011/04/26 10:32:21  bigatti
// -- added tests for multigraded case
//
// Revision 1.13  2011/04/08 14:07:37  bigatti
// -- renamed HilbertNumeratorMod into HilbertNumQuot
//
// Revision 1.12  2011/03/10 17:58:33  bigatti
// -- using long instead of size_t
//
// Revision 1.11  2010/12/17 16:06:25  abbott
// Ensured that all i/o is on standard C++ streams (instead of GlobalInput, etc)
//
// Revision 1.10  2010/10/29 09:43:09  bigatti
// -- manually freeing global memory for C implementation
//    to removed complaints when compiled with -DCoCoA_MEMPOOL_DEBUG
//
// Revision 1.9  2010/10/08 08:19:51  bigatti
// -- RingDistrMPoly.H --> RingDistrMPolyClean.H
//
// Revision 1.8  2009/07/30 15:41:25  bigatti
// -- now using the new nice constructors for ideals
//
// Revision 1.7  2007/10/30 17:14:05  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.6  2007/10/19 10:04:23  bigatti
// -- RingDenseUPolyClean now allow to specify the MinCapacity for all
//    coeff vectors (to avoid too many reallocations)
//
// Revision 1.5  2007/10/18 11:37:33  bigatti
// -- added Q9 for timings tests
//
// Revision 1.4  2007/10/15 12:59:58  bigatti
// -- added include
//
// Revision 1.3  2007/10/15 12:45:59  bigatti
// -- HP computed in Z[lambda] instead of Q[lambda]
//
// Revision 1.2  2007/10/10 14:40:48  bigatti
// -- added: test for comparing old C code with development version using CoCoAlib
//
// Revision 1.1.1.1  2007/03/09 15:16:12  abbott
// Imported files
//
// Revision 1.9  2007/03/08 17:41:36  bigatti
// -- improved NewPolyRing
//
// Revision 1.8  2007/03/08 14:38:07  cocoa
// Added new range function in symbol.H, and tidied many calls to PolyRing
// pseudo ctors (as a consequence).
//
// Revision 1.7  2007/03/03 14:13:21  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.6  2007/03/02 17:46:40  bigatti
// -- unique RingZ and RingQ
// -- requires foundations.H ;  foundations blah;  (thik of a better name)
//
// Revision 1.5  2007/02/26 17:11:58  bigatti
// -- getting ready for unique ring Z: using NewZmod(N), NewRingQ()
//
// Revision 1.4  2007/02/10 18:44:02  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.3  2007/01/17 17:38:11  bigatti
// -- moved all cocoa-4 code for hilbert into src/TmpHilbertDir
//
// Revision 1.2  2006/12/07 17:25:32  cocoa
// -- minimal set of #include's instead of library.H
//
// Revision 1.1  2006/11/16 18:16:10  cocoa
// -- added test for HilbertNumQuot
//
