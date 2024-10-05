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
#include "CoCoA/PPMonoidSparse.H"
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
  matrix NewMatrixFromC(ring Z, int cmat[4][4])
  {
    matrix M(NewDenseMat(Z,4,4));
  
    for (int i=0; i < 4; ++i)
      for (int j=0; j < 4; ++j)
        SetEntry(M, i, j, cmat[i][j]);
    return M;
  }


//-- TestPPMonoid ------------------------------------------------------
// behaviour of different orderings and gradings on PPMonoid and PolyRing

  void TestPPMonoid(const PPMonoid& PPM)
  {
    cout << "TESTING: " << PPM << endl << endl;
    cout << std::boolalpha; // prints true/false for bool

    const vector<PPMonoidElem>& x = indets(PPM);
  
    const PPMonoidElem t1 = x[0]*x[314];
    const PPMonoidElem t2 = x[200];
    cout << "Given t1 = " << t1 << endl;
    cout << "  and t2 = " << t2 << endl << endl;

    cout << "t1*t2*t2*t1   gives  " << t1*t2*t2*t1 << endl;
    cout << "t1*t2*t2*t1*x[385]   gives  " << t1*t2*t2*t1*x[385] << endl;
//   cout << "gcd(t1,t2)   gives  " << gcd(t1,t2) << endl;
//   cout << "lcm(t1,t2)   gives  " << lcm(t1,t2) << endl;
//   cout << "wdeg(t1)   gives  " << wdeg(t1) << endl;
//   cout << "wdeg(t2)   gives  " << wdeg(t2) << endl;
//   cout << "CmpWDeg(t1,t2)   gives  " << CmpWDeg(t1,t2) << endl;
//   if ( GradingDim(PPM) > 0 )
//     cout << "CmpWDegPartial(t1,t2,1)   gives  " << CmpWDegPartial(t1,t2,1) << endl;
//   cout << "t1 > t2   gives  " << (t1 > t2) << endl;

    // Check some computations (without printing any results).
    PPMonoidElem one(PPM);
    CoCoA_ASSERT_ALWAYS(IsOne(one));
    CoCoA_ASSERT_ALWAYS(!IsOne(t1));
    CoCoA_ASSERT_ALWAYS(!IsOne(t2));

    CoCoA_ASSERT_ALWAYS(IsIndet(x[2]));
    CoCoA_ASSERT_ALWAYS(IsIndetPosPower(x[2]));
    CoCoA_ASSERT_ALWAYS(!IsIndet(t1));
    CoCoA_ASSERT_ALWAYS(!IsIndetPosPower(t1));
    PPMonoidElem t1t2 = t1*t2;
    PPMonoidElem t2t1 = t2*t1;
    CoCoA_ASSERT_ALWAYS(t1t2 == t2t1);
    CoCoA_ASSERT_ALWAYS(t1 != t2);
    CoCoA_ASSERT_ALWAYS(t1t2 != t1);
    CoCoA_ASSERT_ALWAYS(t1t2 != t2);
    CoCoA_ASSERT_ALWAYS(IsDivisible(t1t2, t1));
    CoCoA_ASSERT_ALWAYS(IsDivisible(t1t2, t2));
    CoCoA_ASSERT_ALWAYS(!IsDivisible(t1, t2));
    CoCoA_ASSERT_ALWAYS(!IsDivisible(t2, t1));
//   CoCoA_ASSERT_ALWAYS(t1t2/t1 == t2);
//   CoCoA_ASSERT_ALWAYS(t1t2/t2 == t1);
    CoCoA_ASSERT_ALWAYS(IsCoprime(t1,t2));
    CoCoA_ASSERT_ALWAYS(!IsCoprime(t1t2,t1));
    CoCoA_ASSERT_ALWAYS(!IsCoprime(t1t2, t2));
//   CoCoA_ASSERT_ALWAYS(colon(t1, t2) == t1);
//   CoCoA_ASSERT_ALWAYS(colon(t2, t1) == t2);
//   CoCoA_ASSERT_ALWAYS(colon(t1t2, t1) == t2);
//   CoCoA_ASSERT_ALWAYS(colon(t1t2, t2) == t1);
//   CoCoA_ASSERT_ALWAYS(gcd(t1*t1t2, t2*t1t2) == t1t2);
//   CoCoA_ASSERT_ALWAYS(lcm(t1*t1t2, t2*t1t2) == power(t1t2,2));

//   CoCoA_ASSERT_ALWAYS(cmp(t1, one) > 0);
//   CoCoA_ASSERT_ALWAYS(cmp(t2, one) > 0);
//   CoCoA_ASSERT_ALWAYS(cmp(one, t1) < 0);
//   CoCoA_ASSERT_ALWAYS(cmp(one, t2) < 0);
//   CoCoA_ASSERT_ALWAYS(cmp(t1t2, t1) > 0);
//   CoCoA_ASSERT_ALWAYS(cmp(t1t2, t2) > 0);
//   CoCoA_ASSERT_ALWAYS(cmp(t1, t1t2) < 0);
//   CoCoA_ASSERT_ALWAYS(cmp(t2, t1t2) < 0);
//   CoCoA_ASSERT_ALWAYS((t1 > t2)^(t1 <= t2));
//   CoCoA_ASSERT_ALWAYS((t2 > t1)^(t2 <= t1));
//   CoCoA_ASSERT_ALWAYS((t1 >= t2)^(t1 < t2));
//   CoCoA_ASSERT_ALWAYS((t2 >= t1)^(t2 < t1));

//   CoCoA_ASSERT_ALWAYS(power(t1t2, 3) > t1t2);
//   CoCoA_ASSERT_ALWAYS(exponent(power(t1,3),0) == 3*exponent(t1,0));
//   CoCoA_ASSERT_ALWAYS(exponent(power(t1,3),1) == 3*exponent(t1,1));
//   CoCoA_ASSERT_ALWAYS(exponent(power(t1,3),2) == 3*exponent(t1,2));
//   CoCoA_ASSERT_ALWAYS(exponent(power(t1,3),3) == 3*exponent(t1,3));
//   CoCoA_ASSERT_ALWAYS(t1 == IndetPower(PPM, 0, exponent(t1,0))*
//               IndetPower(PPM, 314, exponent(t1,314)));
  
    AssignOne(t2t1);
    CoCoA_ASSERT_ALWAYS(IsOne(t2t1));

//   CoCoA_ASSERT_ALWAYS(wdeg(t1) + wdeg(one) == wdeg(t1));
//   CoCoA_ASSERT_ALWAYS(wdeg(t2) + wdeg(one) == wdeg(t2));
//   CoCoA_ASSERT_ALWAYS(wdeg(t1) + wdeg(t2) == wdeg(t1t2));

    cout << "------------------------------------------------" << endl << endl;
  }


//-- program --------------------------------------------------------------
// we run TestPolyRing on predefined and user-defined orderings

  void program()
  {
    GlobalManager CoCoAFoundations;

    // Each ordering is degree-compatible with grading over Z^GradingDim
    // i.e. the grading is given by the first GradingDim rows of the
    // ordering matrix.
 
    constexpr int n = 400;  // must be >= 386 (see test above)
    const vector<symbol> X = SymbolRange("x", 0, n-1);

    TestPPMonoid(NewPPMonoidEv(X, lex));
    TestPPMonoid(NewPPMonoidSparse(X, lex));
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
