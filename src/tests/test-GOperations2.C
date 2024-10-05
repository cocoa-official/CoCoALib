//   Copyright (c)  2012 Anna M. Bigatti

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
#include "CoCoA/CanonicalHom.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/FreeModule.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/ModuleOrdering.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingDistrMPolyClean.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyOps-ideal.H"
#include "CoCoA/TmpGOperations.H"
#include "CoCoA/VectorOps.H"
#include "CoCoA/submodule.H"
#include "CoCoA/symbol.H"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

//----------------------------------------------------------------------
// Test for ideal/module operations (ex-bugs)
// functions: gcd, colon
// environments: DMPI (QQ)
//----------------------------------------------------------------------
namespace CoCoA
{

  // ex-bug: solved with anonymous symbols
  void test1()
  {
    SparsePolyRing P(NewPolyRing(RingQQ(),symbols("x,y,z")));
    ring K = NewFractionField(P);
    SparsePolyRing Kab(NewPolyRing(K,symbols("a,b")));  
    const RingElem a = indet(Kab,0);
    const RingElem b = indet(Kab,1);
    CoCoA_ASSERT_ALWAYS(gcd(a,b)==1);
  }

  // ex-bug: ideal with 0 generator
  void test2()
  {
    SparsePolyRing P(NewPolyRing(RingQQ(), symbols("x,y")));
    const RingElem x = indet(P,0);
    const RingElem y = indet(P,1);
    CoCoA_ASSERT_ALWAYS(colon(ideal(zero(P)), ideal(x))==ideal(zero(P)));
    CoCoA_ASSERT_ALWAYS(colon(ideal(x*y), ideal(x))==ideal(y));
    CoCoA_ASSERT_ALWAYS(colon(ideal(x*y, x*x), ideal(x))==ideal(y,x));
    CoCoA_ASSERT_ALWAYS(colon(ideal((x+1)*y, y*y), ideal(x+1))==ideal(y));
    CoCoA_ASSERT_ALWAYS(colon(ideal((x+1)*y, (x+1)*(x+1)), ideal(x+1))==ideal(y,x+1));
    CoCoA_ASSERT_ALWAYS(colon(ideal(x+1, zero(P)), ideal(one(P)))==ideal(x+1));
    CoCoA_ASSERT_ALWAYS(colon(ideal((x+1)*y, zero(P)), ideal(zero(P)))==ideal(one(P)));
    CoCoA_ASSERT_ALWAYS(colon(ideal(zero(P)), ideal(zero(P)))==ideal(one(P)));
  }

  void TestMinGens()
  {
    SparsePolyRing P(NewPolyRing(RingQQ(), symbols("x,y,z")));
    const RingElem x = indet(P,0);
    const RingElem y = indet(P,1);
    const RingElem z = indet(P,2);
    CoCoA_ASSERT_ALWAYS(len(MinGens(ideal(x-y, x*(x-y), x*x-y*y)))== 1);
    CoCoA_ASSERT_ALWAYS(len(MinGens(ideal(x-y, x*(x-y), x-z, y-z)))== 2);
  }

  void TestIntersection()
  {
    SparsePolyRing P(NewPolyRing(RingQQ(), symbols("x,y,z")));
    const RingElem x = indet(P,0);
    const RingElem y = indet(P,1);
    const RingElem z = indet(P,2);
    const ideal I = ideal(x-z, x*x-y*y);
    const ideal J = ideal(x-z, x-y);
    const ideal Z = ideal(zero(P));
    CoCoA_ASSERT_ALWAYS(intersect(I,J) == I);
    CoCoA_ASSERT_ALWAYS(intersect(Z,J) == Z);
    CoCoA_ASSERT_ALWAYS(intersect(J,Z) == Z);
  }

  void TestSubmodule1()
  {
    SparsePolyRing P(NewPolyRing(RingQQ(),symbols("x,y,z")));
    const RingElem x = indet(P,0);
    const RingElem y = indet(P,1);
    const RingElem z = indet(P,2);
    const FreeModule FM = NewFreeModule(P, 2, WDegPosnOrd);
    const vector<ModuleElem>& e = gens(FM);
    CoCoA_ASSERT_ALWAYS(IsHomog(submodule(x*e[0] +z*e[1], z*e[0] +y*e[1])));
    CoCoA_ASSERT_ALWAYS(!IsHomog(submodule(x*e[0] +2*e[1], z*e[0] +y*e[1])));
    CoCoA_ASSERT_ALWAYS(IsHomog(submodule(x*e[0] +2*e[1], 3*e[0])));
    CoCoA_ASSERT_ALWAYS(IsHomog(submodule(x*e[0] +2*y*e[1], x*e[0])));
    //  std::cout << TidyGens(submodule(x*(x*e[0] +2*y*e[1]), x*e[0]+2*y*e[1]-e[1]));
  }



  void program()
  {
    GlobalManager CoCoAFoundations;
  
    test1();
    test2();
    TestMinGens();
    TestIntersection();
    TestSubmodule1();
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
