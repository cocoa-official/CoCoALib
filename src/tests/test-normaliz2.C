//   Copyright (c)  2015  John Abbott, Anna Bigatti
//   Orig author:  2012-2015  Christof Soeger

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


//#include "CoCoA/BigInt.H" // included in ??.H
#include "CoCoA/BuildInfo.H"
#include "CoCoA/ExternalLibs-Normaliz.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/TmpPPVector.H"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

//----------------------------------------------------------------------
// Second test for Normaliz library.
// This tests is for the functions working with PPVectors.
//----------------------------------------------------------------------

using namespace CoCoA;

vector<BigInt> BigIntVec(const int* CVector, int len)
{
  vector<BigInt> v(len);
  for (int i=0; i<len; ++i)  v[i] = (BigInt(CVector[i]));
  return v;
}

void program()
{
  GlobalManager CoCoAFoundations;

#ifdef CoCoA_WITH_NORMALIZ
  using namespace CoCoA::Normaliz;
/*
// this is the polytope example from the Normaliz examples
  const int M[4][4] = {{0, 0, 0, 1},
                       {2, 0, 0, 1},
                       {0, 3, 0, 1},
                       {0, 0, 5, 1}};
  const int g[4] = {0, 0, 0, 1};
  vector<vector<BigInt> > l;
  for (int i=0; i<4; ++i)
    l.push_back(BigIntVec(M[i], 4));

  CoCoA_ASSERT_ALWAYS(den.myFactors() == vector<RingElem>(1,RingElem(QQt, "1-t")));
  CoCoA_ASSERT_ALWAYS(den.myMultiplicities() == std::vector<long>(1,4));
  CoCoA_ASSERT_ALWAYS(den.myRemainingFactor() == one(QQt));

  RingElem hp = HilbertPoly(C);
  // compare with Hilbert polynomial:   1 +4*t +8*t^2 +5*t^3
  RingElem ref_hp = RingElem(QQt, "1 + 4*t + 8*t^2 +5*t^3");
  CoCoA_ASSERT_ALWAYS(hp == ref_hp);
*/

  SparsePolyRing R = NewPolyRing(RingQQ(),SymbolRange("x",1,4));
  PPMonoid R_ppm = PPM(R);
  PPVector ppv(R_ppm, NewDivMaskNull());
//JAA  const vector<PPMonoidElem>& x = indets(R_ppm);
  const int E[2][4] = {{-1,-1,2,0}, {1,1,-2,-1}};
  const int C[2][5] = {{1,1,1,1,5}, {1,0,2,0,7}};

  vector<vector<BigInt> > vv_E, vv_C;
  vv_E.push_back(BigIntVec(E[0], 4));
  vv_E.push_back(BigIntVec(E[1], 4));
  vv_C.push_back(BigIntVec(C[0], 5));
  vv_C.push_back(BigIntVec(C[1], 5));
  DiagInvariants(vv_E,vv_C, R_ppm);


#endif // #ifdef CoCoA_WITH_NORMALIZ
}


//----------------------------------------------------------------------
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
  return 1;
}
