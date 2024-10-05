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


#include "CoCoA/BuildInfo.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/MatrixOps.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/error.H"
#include "CoCoA/matrix.H"


#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

//----------------------------------------------------------------------
// Test for matrix solvers
// environments: DenseMatrix over QQ
//----------------------------------------------------------------------
namespace CoCoA
{

// convention: a function containing a "new" should be named "New.."
  matrix MatByRows(ring K, int nrows, int ncols, int* entries)
  {
    matrix M(NewDenseMat(K,nrows,ncols));

    for (int i=0; i < nrows; ++i)
      for (int j=0; j < ncols; ++j)
        SetEntry(M, i, j, *entries++);
    return M;
  }


  void TestLinSolveByGauss(ConstMatrixView M, ConstMatrixView rhs)
  {
    const matrix x = LinSolveByGauss(M, rhs);
    CoCoA_ASSERT_ALWAYS(M*x == rhs);
  }

  void TestLinSolveByGauss_NoSoln(ConstMatrixView M, ConstMatrixView rhs)
  {
    const matrix x = LinSolveByGauss(M, rhs);
    CoCoA_ASSERT_ALWAYS(IsMat0x0(x));
  }


  void invertible()
  {
    const ring QQ = RingQQ();

    // M is vandermonde so surely invertible.
    int M[] = {1, 1, 1, -1,
               0, 1, 2, 1,
               0, 1, 4, -1,
               0, 1, 8, 1};

    int rhs[] = {13, 37,
                 13, -37,
                 -13, 37,
                 -13, -37};

    matrix MM(MatByRows(QQ,4,4, M));
    CoCoA_ASSERT_ALWAYS(NumCols(LinKer(MM)) == 0);
    TestLinSolveByGauss(MM, MatByRows(QQ,4,2, rhs));
  }


  void NotInvertible_SolnExists()
  {
    const ring QQ = RingQQ();

    // M has rank 3
    int M[] = {1, 2, 3, 4,
               5, 6, 7, 8,
               1, 3, 5, 7,
               3, 5, 7, 9};

    int rhs[] = {8,
                 20,
                 13,
                 19};

    matrix MM(MatByRows(QQ,4,4, M));
    CoCoA_ASSERT_ALWAYS(MM*LinKer(MM) == ZeroMat(QQ,4,2));
    TestLinSolveByGauss(MM, MatByRows(QQ,4,1, rhs));
  }


  void NotInvertible_NoSolnExists()
  {
    const ring QQ = RingQQ();

    // M has rank 3
    int M[] = {1, 2, 3, 4,
               5, 6, 7, 8,
               1, 3, 5, 7,
               3, 5, 7, 9};

    int rhs[] = {8,
                 13,
                 19,
                 20};

    TestLinSolveByGauss_NoSoln(MatByRows(QQ,4,4, M),
                               MatByRows(QQ,4,1, rhs));
  }


  void MoreRowsThanCols_SolnExists()
  {
    const ring QQ = RingQQ();

    // M has rank 3
    int M[] = {3, 1, 4,
               1, 5, 9,
               2, 6, 5,
               3, 5, 8,
               9, 7, 9};


    int rhs[] = {17,
                 46,
                 51,
                 49,
                 76};

    TestLinSolveByGauss(MatByRows(QQ,5,3, M),
                        MatByRows(QQ,5,1, rhs));
  }


  void MoreRowsThanCols_NoSolnExists()
  {
    const ring QQ = RingQQ();

    // M has rank 3
    int M[] = {3, 1, 4,
               1, 5, 9,
               2, 6, 5,
               3, 5, 8,
               9, 7, 9};

    int rhs[] = {17,
                 46,
                 51,
                 49,
                 77};

    matrix MM(MatByRows(QQ,5,3, M));
    CoCoA_ASSERT_ALWAYS(NumCols(LinKer(MM)) == 0);
    TestLinSolveByGauss_NoSoln(MM, MatByRows(QQ,5,1, rhs));
  }


  void MoreColsThanRows_SolnExists()
  {
    const ring QQ = RingQQ();

    // M has rank 3
    int M[] = {3, 1, 4, 1, 5,
               9, 2, 6, 5, 3,
               5, 8, 9, 7, 9};

    int rhs[] = {34,
                 18,
                 33};

    matrix MM(MatByRows(QQ,3,5, M));
    CoCoA_ASSERT_ALWAYS(MM*LinKer(MM) == ZeroMat(QQ,3,2));
    TestLinSolveByGauss(MM, MatByRows(QQ,3,1, rhs));
  }


  void MoreColsThanRows_NoSolnExists()
  {
    const ring QQ = RingQQ();

    // M has rank 2
    int M[] = {3, 1, 4, 1, 5,
               9, 2, 6, 5, 3,
               6, 1, 2, 4, -2};

    int rhs[] = {34,
                 18,
                 33};

    TestLinSolveByGauss_NoSoln(MatByRows(QQ,3,5, M),
                               MatByRows(QQ,3,1, rhs));
  }


  void MoreColsThanRows_SolnExists2()
  {
    const ring QQ = RingQQ();

    // M has rank 3
    int M[] = {1, 1, 1, 1, 1,
               1, 1, 1, 1, 2,
               1, 1, 2, 2, 1};

    int rhs[] = {-1,
                 0,
                 1};

    matrix MM(MatByRows(QQ,3,5, M));
    CoCoA_ASSERT_ALWAYS(MM*LinKer(MM) == ZeroMat(QQ,3,2));
    TestLinSolveByGauss(MM, MatByRows(QQ,3,1, rhs));
  }

  void MoreColsThanRows_NoSolnExists2()
  {
    const ring QQ = RingQQ();

    // M has rank 3
    int M[] = {1, 1, 1, 1, 1,
               1, 1, 1, 1, 2,
               2, 2, 2, 2, 3};

    int rhs[] = {-1,
                 0,
                 1};

    matrix MM(MatByRows(QQ,3,5, M));
    CoCoA_ASSERT_ALWAYS(MM*LinKer(MM) == ZeroMat(QQ,3,3));
    TestLinSolveByGauss_NoSoln(MM, MatByRows(QQ,3,1, rhs));
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    invertible();
    NotInvertible_SolnExists();
    NotInvertible_NoSolnExists();
    MoreRowsThanCols_SolnExists();
    MoreRowsThanCols_NoSolnExists();
    MoreColsThanRows_SolnExists();
    MoreColsThanRows_NoSolnExists();
    MoreColsThanRows_SolnExists2();
    MoreColsThanRows_NoSolnExists2();
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
