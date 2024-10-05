//   Copyright (c)  2007,2017  John Abbott and Anna M. Bigatti

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
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/PPMonoidEv.H"
#include "CoCoA/PPMonoidEvOv.H"
#include "CoCoA/PPMonoidOv.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/degree.H"
#include "CoCoA/error.H"
#include "CoCoA/matrix.H"
#include "CoCoA/symbol.H"

#include<iostream>
using std::cout;
using std::cerr;
using std::endl;
#include<vector>
using std::vector;

//----------------------------------------------------------------------
// Test for PPMonoidOv functions on PPMonoidElem
// functions: *, wdeg, CmpWDeg, >, IsOne, ==, IsDivisible, /, IsCoprime, colon, 
//            gcd, lcm, power, cmp, log, IndetPower, wdeg
// environments: PPMonoidEv, PPMonoidOv, PPMonoidEvOv, PPMonoidBigEv;
//               lex, DegLex, DegRevLex, MatOrd
//----------------------------------------------------------------------

namespace CoCoA
{

// convention: a function containing a "new" should be named "New.."
  matrix NewMatrixFromC(ring Z, int cmat[6][6])
  {
    matrix M(NewDenseMat(Z,6,6));
  
    for (int i=0; i < 6; ++i)
      for (int j=0; j < 6; ++j)
        SetEntry(M, i, j, cmat[i][j]);
    return M;
  }


//-- TestPPMonoid ------------------------------------------------------
// behaviour of different orderings and gradings on PPMonoid and PolyRing

  void TestPPMonoid(PPMonoid PPM)
  {
    cout << "TESTING: " << PPM << endl << endl;
    cout << std::boolalpha; // prints true/false for bool

    vector<long> expv(6);
    for (expv[0]=0; expv[0] <= 2; ++expv[0])
      for (expv[1]=0; expv[1] <= 2; ++expv[1])
        for (expv[2]=0; expv[2] <= 2; ++expv[2])
          for (expv[3]=0; expv[3] <= 2; ++expv[3])
            for (expv[4]=0; expv[4] <= 2; ++expv[4])
              for (expv[5]=0; expv[5] <= 2; ++expv[5])
              {
                PPMonoidElem t = PPMonoidElem(PPM, expv);
                PPMonoidElem radt = radical(t);
                CoCoA_ASSERT_ALWAYS(IsSqFree(radt));
                CoCoA_ASSERT_ALWAYS(radical(radt*power(t,3)) == radt);
                if (t == radt)
                {
                  CoCoA_ASSERT_ALWAYS(IsSqFree(t));
                  continue;
                }
                // t is not radical, but since max exp is 2, t/radt must be radical
                CoCoA_ASSERT_ALWAYS(!IsSqFree(t));
                CoCoA_ASSERT_ALWAYS(IsSqFree(t/radt));
              }


    cout << "------------------------------------------------" << endl << endl;
  }


//-- program --------------------------------------------------------------
// we run TestPolyRing on predefined and user-defined orderings

  void program()
  {
    GlobalManager CoCoAFoundations;

    // each ordering is degree-compatible with grading over Z^GradingDim
    // i.e. the grading is given by the first GradingDim rows
    // of the ordering matrix
 
    const int NumIndets = 6;
    const vector<symbol> X = SymbolRange("x", 0, NumIndets-1);

    // user-defined ordering and grading
    const int GradingDim = 2;
    // the first 2 rows represent the degree matrix
    int M[6][6] = {{ 1,  0,  0,  4, 0,  2},
                   { 0,  3,  5,  0, 2,  1},
                   { 1,  0,  0,  0, 0,  0},
                   { 0,  1,  0,  0, 0, -2},
                   {-1, -1, -1, -1, 0, -1},
                   {-1,  0, -2, -4, 1,  0}};
    PPOrdering MatOrd = NewMatrixOrdering(NewMatrixFromC(RingZZ(),M),
                                          GradingDim);
                                          

    TestPPMonoid(NewPPMonoidEvOv(X, lex));
    TestPPMonoid(NewPPMonoidEvOv(X, StdDegLex));
    TestPPMonoid(NewPPMonoidEvOv(X, StdDegRevLex));
    TestPPMonoid(NewPPMonoidEvOv(X, MatOrd));

    TestPPMonoid(NewPPMonoidEv(X, lex));
    TestPPMonoid(NewPPMonoidEv(X, StdDegLex));
    TestPPMonoid(NewPPMonoidEv(X, StdDegRevLex));
    TestPPMonoid(NewPPMonoidEv(X, MatOrd));

    TestPPMonoid(NewPPMonoidOv(X, lex));
    TestPPMonoid(NewPPMonoidOv(X, StdDegLex));
    TestPPMonoid(NewPPMonoidOv(X, StdDegRevLex));
    TestPPMonoid(NewPPMonoidOv(X, MatOrd));

    TestPPMonoid(NewPPMonoidEv(X, lex, PPExpSize::big));
    TestPPMonoid(NewPPMonoidEv(X, StdDegLex, PPExpSize::big));
    TestPPMonoid(NewPPMonoidEv(X, StdDegRevLex, PPExpSize::big));
    TestPPMonoid(NewPPMonoidEv(X, MatOrd, PPExpSize::big));
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
