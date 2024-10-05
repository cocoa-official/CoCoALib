//   Copyright (c)  2013  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/CanonicalHom.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/QuotientRing.H" // for NewZZmod
#include "CoCoA/SparsePolyOps-ideal-points.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/error.H"
#include "CoCoA/ideal.H"
#include "CoCoA/matrix.H"
#include "CoCoA/ring.H"
#include "CoCoA/utils.H"


#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;



namespace CoCoA
{

  void CheckEvalToZero(const vector<RingElem>& g, const matrix& pts)
  {
    const PolyRing& P = owner(g[0]);
    const ring& K = CoeffRing(P);
    const int NumVars = NumIndets(P);
    const int NumPts = NumRows(pts);

    vector<RingHom> EvalMap;
    for (int i=0; i < NumPts; ++i)
    {
      vector<RingElem> images(NumVars, zero(K));
      for (int j=0; j < NumCols(pts); ++j)
        images[j] = pts(i,j);
      EvalMap.push_back(PolyRingHom(P, K, IdentityHom(K), images));
    }

    for (int i=0; i < len(g); ++i)
      for(int j=0; j < len(EvalMap); ++j)
        CoCoA_ASSERT_ALWAYS(IsZero(EvalMap[j](g[i])));
  }

  // Just 2 points: (1,2) and (3,4)
  void test1()
  {
    const ring K = RingQQ();
    const long SpaceDim = 2;
    const SparsePolyRing P = NewPolyRing(K, SymbolRange("x", 1,SpaceDim));
    matrix pts = NewDenseMat(K,2,SpaceDim);
    SetEntry(pts,0,0,1);
    SetEntry(pts,0,1,2);
    SetEntry(pts,1,0,3);
    SetEntry(pts,1,1,4);

    const ideal I = IdealOfPoints(P, pts);
    CheckEvalToZero(gens(I), pts);
  }

  // 10 random points in 3-space
  void test2()
  {
    const ring K = RingQQ();
    constexpr long NumPts = 10;
    constexpr long SpaceDim = 3;
    const SparsePolyRing P = NewPolyRing(K, SymbolRange("x", 1,SpaceDim));
    matrix pts = NewDenseMat(K, NumPts, SpaceDim);
    for (long i=0; i < NumPts; ++i)
      for (long j=0; j < SpaceDim; ++j)
        SetEntry(pts,i,j, RandomLong(-99,99));

    const ideal I = IdealOfPoints(P, pts);
    CheckEvalToZero(gens(I), pts);
  }

  // 10 random points in 3-space -- mapping points from QQ
  void test3()
  {
    const ring K = NewZZmod(32003);
    constexpr long NumPts = 10;
    constexpr long SpaceDim = 3;
    const SparsePolyRing P = NewPolyRing(K, SymbolRange("x", 1,SpaceDim));
    matrix pts = NewDenseMat(RingQQ(), NumPts, SpaceDim);
    for (long i=0; i < NumPts; ++i)
      for (long j=0; j < SpaceDim; ++j)
        SetEntry(pts,i,j, RandomLong(-99,99));

    const ideal I = IdealOfPoints(P, pts);
    CheckEvalToZero(gens(I), CanonicalHom(RingQQ(),K) (pts));
  }

  void program()
  {
    GlobalManager CoCoAFoundations;

    test1();
    test2();
    test3();
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
