//   Copyright (c)  2007  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingDistrMPolyClean.H"
#include "CoCoA/RingDistrMPolyInlFpPP.H"
#include "CoCoA/RingDistrMPolyInlPP.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/ideal.H"


#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

//----------------------------------------------------------------------
// Test for ideals in (commutative) SparsePolyRings.
// functions: [GBasis], NF, IsElem, IsContained, +, *
// environments: DMP (Fp, Q), DMPI (Q), DMPII
//----------------------------------------------------------------------
namespace CoCoA
{

  // test ideal ops: NF, IsElem, IsContained, IsZero, intersect, op==

  void TestRing(SparsePolyRing P)
  {
    const vector<RingElem>& x = indets(P);
    vector<RingElem> g;
    g.push_back(power(x[1],3) - power(x[3],3));
    g.push_back(power(x[1],2) - x[2]);
    g.push_back(power(x[3],4));
    g.push_back(power(x[0],6)*x[1] - power(x[0],5)*power(x[3],2));
    g.push_back(power(x[0],5)*power(x[1],2) - power(x[0],5)*power(x[3],3));

    vector<RingElem> f;
    f.push_back(g[0]);
    f.push_back(power(x[3],5));

    ideal I = ideal(g);
    ideal J = ideal(f);

    CoCoA_ASSERT_ALWAYS( IsDivisible(NF(power(x[1],2), I), x[2]) );
    CoCoA_ASSERT_ALWAYS( IsDivisible(NF(power(x[3],3), I), x[1]*x[2]) );

    CoCoA_ASSERT_ALWAYS( !IsElem(power(x[1],2), I) );
    CoCoA_ASSERT_ALWAYS( !IsElem(power(x[3],3), I) );
    CoCoA_ASSERT_ALWAYS(  IsElem(g[0], I) );

    CoCoA_ASSERT_ALWAYS( !IsContained(I, J) );
    CoCoA_ASSERT_ALWAYS( !IsContained(ideal(g), ideal(f)) );
    CoCoA_ASSERT_ALWAYS(  IsContained(J, I) ); // both with and without precomputed GB
    CoCoA_ASSERT_ALWAYS(  IsContained(ideal(f), ideal(g)) );
    CoCoA_ASSERT_ALWAYS( !IsOne(I) );
    f.push_back(-one(P));
    CoCoA_ASSERT_ALWAYS( IsOne(ideal(f)) );

    CoCoA_ASSERT_ALWAYS( ideal(x[0])+ideal(x[1]) == ideal(x[0], x[1]) );
    CoCoA_ASSERT_ALWAYS( ideal(x[0],x[1])*ideal(x[1],x[2])
                 == ideal(x[0]*x[1], x[0]*x[2], x[1]*x[1], x[1]*x[2]) );
    CoCoA_ASSERT_ALWAYS( power(ideal(x[0],x[1]), 3)
                 ==  ideal(x[0],x[1])*ideal(x[0],x[1])*ideal(x[0],x[1]));
    CoCoA_ASSERT_ALWAYS( IsZero(intersect(I, ideal(zero(P)))));
    CoCoA_ASSERT_ALWAYS( IsZero(intersect(ideal(zero(P)), J)));
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    const ring Fp = NewZZmod(101);
    const ring QQ = RingQQ();

    const SparsePolyRing Fpx = NewPolyRing(Fp, SymbolRange("x",0,3));
    TestRing(Fpx);

    const SparsePolyRing QQx = NewPolyRing(QQ, SymbolRange("x",0,3));
    TestRing(QQx);

    const SparsePolyRing FpxII = NewPolyRing_DMPII(Fp, SymbolRange("x",0,3));
    TestRing(FpxII);

    const SparsePolyRing QQxI = NewPolyRing_DMPI(QQ, SymbolRange("x",0,3));
    TestRing(QQxI);
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
