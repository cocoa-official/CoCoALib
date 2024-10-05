//   Copyright (c)  2013-2015  John Abbott,  Anna M. Bigatti
//   Original author: 2013-2015  Mario Albert

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
#include "CoCoA/GlobalManager.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/TmpResolutionMinimization.H"
#include "CoCoA/matrix.H"
#include "CoCoA/symbol.H"
#include "CoCoA/error.H"

using namespace CoCoA;

#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

// matrix NewMatrixFromC(ring R, RingElem cmat[4][4])
// {
//   matrix M(NewDenseMat(R,4,4));

//   for (int i=0; i < 4; ++i)
//     for (int j=0; j < 4; ++j)
//       SetEntry(M, i, j, cmat[i][j]);
//   return M;
// }

// void test_myFindPivot()
// {
//   ring Q = RingQQ();
//   SparsePolyRing PRing = NewPolyRing(Q, SymbolRange("x",0,4), StdDegLex);

//   const vector<RingElem> x = indets(PRing);

//   RingElem C_matrix[4][4] = {{RingElem(PRing, x[0] + x[1]), RingElem(PRing, x[0]), RingElem(PRing, x[0] + 0), RingElem(PRing, x[0])},
//                              {RingElem(PRing, x[0] + x[2]), RingElem(PRing, x[0]), RingElem(PRing, x[0] + 0), RingElem(PRing, x[0])},
//                              {RingElem(PRing, x[0] + x[3]), RingElem(PRing, x[0]), RingElem(PRing, x[0] + 0), RingElem(PRing, x[0])},
//                              {RingElem(PRing, x[0]), RingElem(PRing, x[0]), RingElem(PRing, x[0] + 0), RingElem(PRing, x[0])}};

// ???? MAKE myFindPivot PRIVATE ????
//   ResolutionMinimization min(PRing, vector<matrix>());
//   matrix test(NewMatrixFromC(PRing, C_matrix));

//   std::pair<long, long> pivot(min.myFindPivot(test));
//   CoCoA_ASSERT_ALWAYS(pivot.first == -1 && pivot.second == -1);

//   SetEntry(test, 2, 2, zero(PRing));
//   SetEntry(test, 1, 1, zero(PRing));
//   SetEntry(test, 3, 3, zero(PRing));
//   pivot = min.myFindPivot(test);
//   CoCoA_ASSERT_ALWAYS(pivot.first == -1 && pivot.second == -1);

//   SetEntry(test, 0, 2, x[0] + x[1] * x[2]);
//   SetEntry(test, 0, 3, 1);
//   SetEntry(test, 1, 0, 1);
//   pivot = min.myFindPivot(test);
//   CoCoA_ASSERT_ALWAYS(pivot.first == 0 && pivot.second == 3);

//   SetEntry(test, 0, 3, x[0] + x[1]);
//   SetEntry(test, 1, 0, one(PRing));
//   pivot = min.myFindPivot(test);
//   CoCoA_ASSERT_ALWAYS(pivot.first == 1 && pivot.second == 0);
// }

//  void test_myManipulateMatrix()
// {
//   //assume that myCancelColumn is correct
//   ring Q = RingQQ();
//   SparsePolyRing PRing = NewPolyRing(Q, SymbolRange("x",0,4), StdDegLex);

//   const vector<RingElem> x = indets(PRing);

//   RingElem C_matrix[4][4] = {{RingElem(PRing, x[0]), RingElem(PRing, x[0]), RingElem(PRing, x[0]), RingElem(PRing, x[0])},
//                              {RingElem(PRing, x[0]), RingElem(PRing, x[0]), RingElem(PRing, x[0]), RingElem(PRing, x[0])},
//                              {RingElem(PRing, x[0]), RingElem(PRing, x[0]), RingElem(PRing, x[0]), RingElem(PRing, x[0])},
//                              {RingElem(PRing, x[0]), RingElem(PRing, x[0]), RingElem(PRing, x[0]), RingElem(PRing, x[0])}};
//   RingElem Test_C_matrix[4][4] = {{RingElem(PRing, x[0]), RingElem(PRing, x[0]), RingElem(PRing, x[0]), RingElem(PRing, x[0])},
//                              {RingElem(PRing, 0), RingElem(PRing, 0), RingElem(PRing, 0), RingElem(PRing, 0)},
//                              {RingElem(PRing, 0), RingElem(PRing, 0), RingElem(PRing, 0), RingElem(PRing, 0)},
//                              {RingElem(PRing, 0), RingElem(PRing, 0), RingElem(PRing, 0), RingElem(PRing, 0)}};

// 2013-06-30 JAA DISABLES THIS TEST because PivotColumn is IMPOSSIBLE!!
// PivotColumn must have a "constant" entry!!!
//   ResolutionMinimization min(PRing, vector<matrix>());
//   matrix manipulating(NewMatrixFromC(PRing, C_matrix));
//   matrix test(NewMatrixFromC(PRing, Test_C_matrix));
//   std::vector<RingElem> PivotColumn;
//   PivotColumn.push_back(x[0]);
//   PivotColumn.push_back(x[0]);
//   PivotColumn.push_back(x[0]);
//   PivotColumn.push_back(x[0]);

//   min.myManipulateMatrix(manipulating, 0, PivotColumn);

  // matrix ResMan(manipulating);
  // CoCoA_ASSERT_ALWAYS(test == ResMan);

  // RingElem CC_matrix[4][4] = {{RingElem(PRing, x[1]), RingElem(PRing, x[0]), RingElem(PRing, x[0]), RingElem(PRing, x[0])},
  //                             {RingElem(PRing, x[0]), RingElem(PRing, x[2]), RingElem(PRing, x[3]), RingElem(PRing, x[2])},
  //                             {RingElem(PRing, x[0]), RingElem(PRing, x[2]), RingElem(PRing, x[3]), RingElem(PRing, x[2])},
  //                             {RingElem(PRing, x[0]), RingElem(PRing, x[2]), RingElem(PRing, x[3]), RingElem(PRing, x[2])}};
  // RingElem Test_CC_matrix[4][4] = {{RingElem(PRing, x[1])       , RingElem(PRing, x[0]), RingElem(PRing, x[0]), RingElem(PRing, x[0])},
  //                             {RingElem(PRing, x[0] - x[1]), RingElem(PRing, x[2] - x[0]), RingElem(PRing, x[3] - x[0]), RingElem(PRing, x[2] - x[0])},
  //                             {RingElem(PRing, x[0] - x[1]), RingElem(PRing, x[2] - x[0]), RingElem(PRing, x[3] - x[0]), RingElem(PRing, x[2] - x[0])},
  //                             {RingElem(PRing, x[0] - x[1]), RingElem(PRing, x[2] - x[0]), RingElem(PRing, x[3] - x[0]), RingElem(PRing, x[2] - x[0])}};


  // matrix manipulating2(NewMatrixFromC(PRing, CC_matrix));
  // matrix test2(NewMatrixFromC(PRing, Test_CC_matrix));
  // min.myManipulateMatrix(manipulating2, 0, PivotColumn);
  // ResMan = manipulating2;
  // CoCoA_ASSERT_ALWAYS(test2 == ResMan);
// }


void test_myMinimization()
{
  //example script Computer Algebra (Janko Boehm) p. 153 ff.
  ring Q = RingQQ();
  SparsePolyRing PRing = NewPolyRing(Q, SymbolRange("x",0,6), StdDegLex);

  const vector<RingElem> x = indets(PRing);

  matrix phi1(NewDenseMat(PRing, 1, 5));
  SetEntry(phi1, 0, 0, x[0] * x[1] * x[4]);
  SetEntry(phi1, 0, 1, x[0] * x[1] * x[5]);
  SetEntry(phi1, 0, 2, x[2] * x[3] * x[5]);
  SetEntry(phi1, 0, 3, x[2] * x[3] * x[6]);
  SetEntry(phi1, 0, 4, x[4] * x[6]);
  matrix phi2(NewDenseMat(PRing, 5, 6));
  SetEntry(phi2, 0, 0, -x[5] );
  SetEntry(phi2, 0, 3, -x[6]);
  SetEntry(phi2, 1, 0, x[4]);
  SetEntry(phi2, 1, 4, -x[2] * x[3]);
  SetEntry(phi2, 1, 5, -x[4] * x[6]);
  SetEntry(phi2, 2, 1, -x[6]);
  SetEntry(phi2, 2, 4, x[0] * x[1]);
  SetEntry(phi2, 3, 1, x[5]);
  SetEntry(phi2, 3, 2, -x[4]);
  SetEntry(phi2, 4, 2, x[2] * x[3]);
  SetEntry(phi2, 4, 3, x[0] * x[1]);
  SetEntry(phi2, 4, 5, x[0] * x[1] * x[5]);
  matrix phi3(NewDenseMat(PRing, 6, 2));
  SetEntry(phi3, 0 , 0, x[6]);
  SetEntry(phi3, 1 , 1, -x[0] * x[1] * x[4]);
  SetEntry(phi3, 2 , 1, -x[0] * x[1] * x[5]);
  SetEntry(phi3, 3 , 0, -x[5]);
  SetEntry(phi3, 4 , 1, -x[4] * x[6]);
  SetEntry(phi3, 5 , 0, one(PRing));
  SetEntry(phi3, 5 , 1, x[2] * x[3]);

  vector<matrix> resolution;
  resolution.push_back(phi1);
  resolution.push_back(phi2);
  resolution.push_back(phi3);

  ResolutionMinimization min(PRing, resolution);
  min.myMinimization();

  matrix test_phi1(min.myGetResolution(0));
  matrix test_phi2(min.myGetResolution(1));
  matrix test_phi3(min.myGetResolution(2));

  CoCoA_ASSERT_ALWAYS(test_phi1 == phi1);

  CoCoA_ASSERT_ALWAYS(test_phi2(0, 0) == -x[5] && test_phi2(0, 1) == 0      && test_phi2(0, 2) == 0           && test_phi2(0, 3) == -x[6]       && test_phi2(0, 4) == 0           );
  CoCoA_ASSERT_ALWAYS(test_phi2(1, 0) ==  x[4] && test_phi2(1, 1) == 0      && test_phi2(1, 2) == 0           && test_phi2(1, 3) == 0           && test_phi2(1, 4) == -x[2] * x[3]);
  CoCoA_ASSERT_ALWAYS(test_phi2(2, 0) == 0     && test_phi2(2, 1) == -x[6]  && test_phi2(2, 2) == 0           && test_phi2(2, 3) == 0           && test_phi2(2, 4) == x[0] * x[1] );
  CoCoA_ASSERT_ALWAYS(test_phi2(3, 0) == 0     && test_phi2(3, 1) ==  x[5]  && test_phi2(3, 2) == -x[4]       && test_phi2(3, 3) == 0           && test_phi2(3, 4) == 0           );
  CoCoA_ASSERT_ALWAYS(test_phi2(4, 0) == 0     && test_phi2(4, 1) == 0      && test_phi2(4, 2) == x[2] * x[3] && test_phi2(4, 3) == x[0] * x[1] && test_phi2(4, 4) == 0           );

  CoCoA_ASSERT_ALWAYS(test_phi3(0, 0) == -x[2] * x[3] * x[6]);
  CoCoA_ASSERT_ALWAYS(test_phi3(1, 0) == -x[0] * x[1] * x[4]);
  CoCoA_ASSERT_ALWAYS(test_phi3(2, 0) == -x[0] * x[1] * x[5]);
  CoCoA_ASSERT_ALWAYS(test_phi3(3, 0) ==  x[2] * x[3] * x[5]);
  CoCoA_ASSERT_ALWAYS(test_phi3(4, 0) == -x[4] * x[6]       );
}

void program()
{
  GlobalManager CoCoAFoundations;
//   test_myFindPivot();
//   test_myManipulateMatrix();
  test_myMinimization();
}


// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
  try
  {
    program();
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

  BuildInfo::PrintAll(cerr);
  return 1;
}
