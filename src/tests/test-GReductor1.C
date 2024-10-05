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
#include "CoCoA/FractionField.H"
#include "CoCoA/FreeModule.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/ModuleOrdering.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingDistrMPolyClean.H"
#include "CoCoA/RingDistrMPolyInlFpPP.H"
#include "CoCoA/RingDistrMPolyInlPP.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/TmpGOperations.H"
#include "CoCoA/TmpGPoly.H"
#include "CoCoA/VectorOps.H"
#include "CoCoA/submodule.H"
#include "CoCoA/symbol.H"

#include <algorithm>
using std::find_if;
#include <functional>
using std::not1;
using std::ptr_fun;
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::boolalpha;
#include <string>
#include <vector>
using std::vector;

//----------------------------------------------------------------------
// Test for Buchberger's Algorithm 
// functions: ComputeGBasis
// environments: DMP (Fp, Q), DMPI (Fp), DMPII
//----------------------------------------------------------------------
namespace CoCoA
{

  void PrintInfo(const std::string& s, const SparsePolyRing& P, const vector<RingElem>& InputPolys)
  {
    cout << "Executing test " << s << endl << P << "   ";
    cout << "GradingDim = " << GradingDim(P) << "   ";
    if (GradingDim(P)!=0)
      cout << "IsHomog(InputPolys) = " << IsHomog(InputPolys);
    cout << endl;
  }


  // expecting P to have 5 indeterminates
  void testC4_h(SparsePolyRing P)
  {
    RingElem x = indet(P,0),  y = indet(P,1),  z = indet(P,2), t=indet(P,3);
    RingElem h = indet(P,4);
    PolyList InputPolys,GB;
    InputPolys.push_back(x+y+z+t);
    InputPolys.push_back(x*y+y*z+z*t+t*x);
    InputPolys.push_back(x*y*z+y*z*t+z*t*x+t*x*y);
    InputPolys.push_back(x*y*z*t-power(h,4));
    PrintInfo("C4_h", P, InputPolys);
    PolyList MinGens;
    ComputeGBasis(GB, MinGens, InputPolys);
    CoCoA_ASSERT_ALWAYS(len(MinGens)==4);
    monic(GB);
    cout << "C4_h = " << GB << "\n"<<endl;
  }


  void testC4(SparsePolyRing P)
  {
    RingElem x = indet(P,0),  y = indet(P,1),  z = indet(P,2), t = indet(P,3);
    PolyList InputPolys,GB;
    InputPolys.push_back(x+y+z+t);
    InputPolys.push_back(x*y + y*z + z*t + t*x);
    InputPolys.push_back(x*y*z + y*z*t + z*t*x + t*x*y);
    InputPolys.push_back(x*y*z*t - 1);
    PrintInfo("C4", P, InputPolys);
    PolyList MinGens;
    ComputeGBasis(GB, MinGens, InputPolys);
    CoCoA_ASSERT_ALWAYS(len(MinGens)==0);
    monic(GB);
    cout << "C4 = " << GB << "\n"<<endl;
  }

  // expecting P to have 4 indeterminates
  void testCurve_h(SparsePolyRing P)
  {
    RingElem t = indet(P,0),  x = indet(P,1),  y = indet(P,2), z=indet(P,3);
    PolyList InputPolys,GB;
    InputPolys.push_back(power(t,3)-power(x,3));
    InputPolys.push_back(power(t,5)-power(y,5));
    InputPolys.push_back(power(t,7)-power(z,7));
    PrintInfo("Curve_h", P, InputPolys);
    PolyList MinGens;
    ComputeGBasis(GB, MinGens, InputPolys);
    CoCoA_ASSERT_ALWAYS(len(MinGens)==3);
    monic(GB);
    cout << "Curve_h = " << GB << "\n"<<endl;
  }

  // expecting P to have 5 indeterminates
  void testCurve(SparsePolyRing P)
  {
    const vector<RingElem>& x = indets(P);
    PolyList InputPolys, GB;
    InputPolys.push_back(power(x[0],3) - x[1]);
    InputPolys.push_back(power(x[0],5) - x[2]);
    InputPolys.push_back(power(x[0],7) - x[3]);
    PrintInfo("Curve", P, InputPolys);
    PolyList MinGens;
    ComputeGBasis(GB, MinGens, InputPolys);
    CoCoA_ASSERT_ALWAYS(len(MinGens)==0);
    monic(GB);
    cout << "Curve = " << GB << "\n" << endl;
  }

  void program()
  {
    GlobalManager CoCoAFoundations(UseNonNegResidues);
    cout << boolalpha;
  
    // This is a test of the correctness of GB computations for the three cases:
    // AllGraded, AllAffine and GradedRing, Inhomogenous Input for ideals only.
    // This does not check if the computations have been done in the right way,
    // only if the result is correct.

    // the results have been tested with cocoa 4
  
    ring ZZmod = NewZZmod(101);  
  
    // Cyclic 4
    const vector<symbol> X = SymbolRange("x", 0, 4);
    SparsePolyRing P(NewPolyRing_DMPI(ZZmod, X));
    SparsePolyRing Q(NewPolyRing_DMPII(ZZmod, X));
    SparsePolyRing R(NewPolyRing(RingQQ(), NewSymbols(5)));
    SparsePolyRing T(NewPolyRing_DMPI(RingQQ(), X));

    // Graded Ring, Homogeneous Input
    testC4_h(P);
    testC4_h(Q);
    testC4_h(R);
    testC4_h(T);
 
    // Graded Ring, Non-homogeneous Input
    testC4(P);
    testC4(Q);
    testC4(R);
    testC4(T);
  
    SparsePolyRing LexR = NewPolyRing(RingQQ(), X, lex);
    // Non-graded Ring, Non-homogeneous Input
    //testC4(LexR,"NewPolyRing Q");
    
    // Curve
    const vector<symbol> Y = SymbolRange("x", 0 ,3);
    SparsePolyRing P1(NewPolyRing_DMPI(ZZmod, Y));
    SparsePolyRing Q1(NewPolyRing_DMPII(ZZmod, Y));
    SparsePolyRing R1(NewPolyRing(RingQQ(), NewSymbols(4)));
    SparsePolyRing T1 (NewPolyRing_DMPI(RingQQ(), Y));

    // Graded Ring, Homogeneous Input
    testCurve_h(P1);
    testCurve_h(Q1);
    testCurve_h(R1);
    testCurve_h(T1);
  
    // Graded Ring, Non-homogeneous Input
    testCurve(P1);
    testCurve(Q1);
    testCurve(R1);
    testCurve(T1);
  
    SparsePolyRing LexR1 = NewPolyRing(RingQQ(), Y, lex);

    // Non-graded Ring, Non-homogeneous Input
    testCurve(LexR1);
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
