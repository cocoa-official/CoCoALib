//   Copyright (c)  2008  Anna Bigatti

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
#include "CoCoA/SparsePolyOps-ideal.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/symbol.H"
// only for testing:
#include "CoCoA/time.H"
#include "CoCoA/VectorOps.H"


#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

//----------------------------------------------------------------------
// Test for monomial ideal operations
// functions: GBasis, intersect // todo NF, IsElem, IsContained
// environments: DMP (Fp, Q), DMPI (Q), DMPII
//----------------------------------------------------------------------
namespace CoCoA
{

  void TestSparsePolyRing(SparsePolyRing P)
  {
    //  cout << "TESTING: " << P << endl << endl;

    RingElem x = indet(P,0);
    RingElem y = indet(P,1);
    RingElem z = indet(P,2);

    ideal I = ideal(x*y, y*z);
    ideal J = ideal(x*x, x*y*y);
    //  cout << I << endl;
    //  cout << J << endl;
    J = intersect(I, J);
    //  cout << J << endl;
    //  cout << TidyGens(J) << endl;
    CoCoA_ASSERT_ALWAYS(J == ideal(x*y*y, x*x*y));
    const vector<RingElem>& g = TidyGens(J);
    CoCoA_ASSERT_ALWAYS(g.size() == 2);
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    const QuotientRing Fp = NewZZmod(101);
    const SparsePolyRing QQxyz = NewPolyRing(RingQQ(), symbols("x,y,z"));
    const SparsePolyRing QQx = NewPolyRing_DMPI(RingQQ(), SymbolRange("x",1,4));
    const SparsePolyRing Fpy = NewPolyRing_DMPI(Fp, SymbolRange("y",1,3));
    const SparsePolyRing FpyII = NewPolyRing_DMPII(Fp, SymbolRange("y",1,4));

    TestSparsePolyRing(QQxyz);
    TestSparsePolyRing(Fpy);
    TestSparsePolyRing(FpyII);
    TestSparsePolyRing(QQx);

    {

      ring Fp = NewZZmod(32003);          // coefficient ring
      SparsePolyRing Fpx = NewPolyRing(Fp, SymbolRange("x", 0,7));
      SparsePolyRing P = Fpx;
      double t0;  // for CpuTime  
      bool IsPrintingMode = false;
      //  IsPrintingMode = true;

      const vector<RingElem>& x = indets(P);
      vector<RingElem> g;
      back_inserter(g) = power(x[2],2) * power(x[5],4);
      back_inserter(g) = power(x[1],3) * power(x[4],4);
      back_inserter(g) = power(x[1],3) * power(x[5],4);
      back_inserter(g) = power(x[3],3) * power(x[6],4);
      back_inserter(g) = power(x[3],4) * power(x[6],3);

      ideal J1(g);
      ideal J2(x[1]*x[2], x[2]*x[3]);
      if (IsPrintingMode) cout << "J1  = " << J1 << endl; 
      if (IsPrintingMode) cout << "J2  = " << J2 << endl << endl;

      CoCoA_ASSERT_ALWAYS(AreGensMonomial(J1));
      CoCoA_ASSERT_ALWAYS(!AreGensSqFreeMonomial(J1));

      CoCoA_ASSERT_ALWAYS(AreGensMonomial(J2));
      CoCoA_ASSERT_ALWAYS(AreGensSqFreeMonomial(J2));

      CoCoA_ASSERT_ALWAYS(!AreGensMonomial(ideal(x[0],x[1]-1)));
      CoCoA_ASSERT_ALWAYS(!AreGensSqFreeMonomial(ideal(x[0],x[1]-1)));

      CoCoA_ASSERT_ALWAYS(AreGensMonomial(ideal(x[0],x[1]*x[1])));
      CoCoA_ASSERT_ALWAYS(!AreGensSqFreeMonomial(ideal(x[0],x[1]*x[1])));

      t0 = CpuTime();
      ideal I = intersect(J1, J2);
      if (IsPrintingMode) cout << "Cpu Time = " << CpuTime()-t0 << endl;
      if (IsPrintingMode) cout << "intersect(J1, J2) = " << I << endl;
      if (IsPrintingMode) cout << endl;
      CoCoA_ASSERT_ALWAYS(I ==
                  ideal(power(x[2],2)*x[3]*power(x[5],4), x[1]*power(x[2],2)*power(x[5],4),
                        x[2]*power(x[3],3)*power(x[6],4), x[2]*power(x[3],4)*power(x[6],3)) +
                  ideal(power(x[1],3)*x[2]*power(x[5],4), power(x[1],3)*x[2]*power(x[4],4)));
  

      t0 = CpuTime();
      vector<ideal> PrimDec = PrimaryDecomposition(J2);
      if (IsPrintingMode) cout << "Cpu Time = " << CpuTime()-t0 << endl;
      if (IsPrintingMode) cout << "PrimaryDecomposition(J2) = " << PrimDec << endl;
      CoCoA_ASSERT_ALWAYS(PrimDec[0] == ideal(x[2]));
      CoCoA_ASSERT_ALWAYS(PrimDec[1] == ideal(x[1], x[3]));

      t0 = CpuTime();
      I = J1 * J2;
      if (IsPrintingMode) cout << "Cpu Time = " << CpuTime()-t0 << endl;
      if (IsPrintingMode) cout << "J1 * J2 = " << I << endl;
      if (IsPrintingMode) cout << endl;

      t0 = CpuTime();
      I = colon(J1, J2);
      if (IsPrintingMode) cout << "Cpu Time = " << CpuTime()-t0 << endl;
      if (IsPrintingMode) cout << "colon(J1, J2) = " << I << endl;
      if (IsPrintingMode) cout << endl;
      CoCoA_ASSERT_ALWAYS(I ==
                  ideal(x[2]*power(x[5],4), power(x[3],3)*power(x[6],4),
                        power(x[3],4)*power(x[6],3), power(x[1],3)*power(x[5],4))+
                  ideal(power(x[1],3)*power(x[4],4),
                        power(x[1],2)*power(x[3],2)*power(x[5],4)*power(x[6],4),
                        power(x[1],2)*power(x[3],2)*power(x[4],4)*power(x[6],4),
                        power(x[1],2)*power(x[3],3)*power(x[5],4)*power(x[6],3))+
                  ideal(power(x[1],2)*power(x[3],3)*power(x[4],4)*power(x[6],3)));

      t0 = CpuTime();
      I = elim(J1, {x[1],x[2]});
      if (IsPrintingMode) cout << "Cpu Time = " << CpuTime()-t0 << endl;
      if (IsPrintingMode) cout << "elim(ElimInds, J2) = " << I << endl;
      if (IsPrintingMode) cout << endl;
      CoCoA_ASSERT_ALWAYS(I ==
                  ideal(power(x[3],3)*power(x[6],4),
                        power(x[3],4)*power(x[6],3)));
    }

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
