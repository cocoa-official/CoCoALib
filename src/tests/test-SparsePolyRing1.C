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


#include "CoCoA/BigIntOps.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/combinatorics.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/MatrixForOrdering.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingDistrMPolyClean.H"
#include "CoCoA/RingDistrMPolyInlFpPP.H"
#include "CoCoA/RingDistrMPolyInlPP.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/degree.H"
#include "CoCoA/matrix.H"
#include "CoCoA/symbol.H"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <vector>
using std::vector;



//----------------------------------------------------------------------
// Test for SparsePolyRing functions on RingElem
// functions: indets, *, wdeg, StdDeg, IsHomog, homog, gcd, CmpWDeg
// environments: DMP (Fp, ZZ), DMPI (Q), DMPII
//----------------------------------------------------------------------
namespace CoCoA
{

  void TestSparsePolyRing(SparsePolyRing P)
  {
    //    cout << "TESTING: " << P << endl << endl;
    cout << "TESTING:" << endl;
    P->myOutputSelfLong(cout);
    cout << endl << endl;
    cout << std::boolalpha; // prints true/false for bool

  
    const vector<RingElem>& x = indets(P);
    RingElem f1 = x[0]*x[0] + 3*x[1] -2,
      f2 = x[1]*x[1] - 5*x[0],
      f3 = x[0]-x[1];
    ConstRefRingElem h = x[NumIndets(P)-1];

    cout << "Given f1 = " << f1 << endl;
    cout << "  and f2 = " << f2 << endl << endl;
 
    cout << "f1*f2   gives  " << f1*f2 << endl;
    //  cout << "gcd(f1,f2)   gives  " << gcd(f1,f2) << endl;
    cout << "wdeg(f1)   gives  " << wdeg(f1) << endl;
    cout << "wdeg(f2)   gives  " << wdeg(f2) << endl;
    cout << "StdDeg(f1)   gives  " << StdDeg(f1) << endl;
    cout << "StdDeg(f2)   gives  " << StdDeg(f2) << endl;
    cout << "IsHomog(f1)   gives  " << IsHomog(f1) << endl;
    cout << "IsHomog(f2)   gives  " << IsHomog(f2) << endl;
    cout << "IsHomog(f3)   gives  " << IsHomog(f3) << endl;
    cout << "CmpWDeg(f1,f2)        gives  " << sign(CmpWDeg(f1,f2)) << endl; // CmpWDeg guarantees only the sign of result
    if (GradingDim(P)>0)
    {
      cout << "LF(f3)   gives  " << LF(f3) << endl;
      cout << "CmpWDegPartial(f1,f2) gives  "
           << sign(CmpWDegPartial(f1,f2,1)) << endl;  // CmpWDeg guarantees only the sign of result
      cout << "IsHomogPartial(f3,1)  gives  "
           << IsHomogPartial(f3,1) << endl;
    }
    if (GradingDim(P)<2)
    {
      cout << "  -- homogenizing with h = " << h << endl;
      cout << "homog(f1, h)   gives  " << homog(f1,h) << endl;
      cout << "homog(f2, h)   gives  " << homog(f2,h) << endl;
      cout << "homog(f3, h)   gives  " << homog(f3,h) << endl;
    }

    if ( IsField(CoeffRing(P)) || IsTrueGCDDomain(CoeffRing(P)))
    {
      CoCoA_ASSERT_ALWAYS( IsInvertible(gcd(f1*f2, f1*f3)/f1) );
      CoCoA_ASSERT_ALWAYS( IsInvertible(gcd(f1*f2, f2*f3)/f2) );
      CoCoA_ASSERT_ALWAYS( IsInvertible(gcd(f1*f3, f2*f3)/f3) );
    }
    if ( IsFractionFieldOfGCDDomain(CoeffRing(P)) )
    {
      CoCoA_ASSERT_ALWAYS( CommonDenom(x[0]/3+x[1]/2) ==  6);
      CoCoA_ASSERT_ALWAYS( ClearDenom(x[0]/3+x[1]/2) ==  2*x[0]+3*x[1]);
    }
  
    // Check some computations (without printing any results).
    //    cout << "one(P)  " << one(P) << endl;
    CoCoA_ASSERT_ALWAYS(IsOne(one(P)));
    CoCoA_ASSERT_ALWAYS(!IsOne(f1));
    CoCoA_ASSERT_ALWAYS(!IsOne(f2));

    RingElem f1f2 = f1*f2;
    RingElem f2f1 = f2*f1;
    CoCoA_ASSERT_ALWAYS(f2*2 == 2*x[1]*x[1] - 10*x[0]);
    CoCoA_ASSERT_ALWAYS(f1f2 == f2f1);
    CoCoA_ASSERT_ALWAYS(f1 != f2);
    CoCoA_ASSERT_ALWAYS(f1f2 != f1);
    CoCoA_ASSERT_ALWAYS(f1f2 != f2);
    if ( IsField(CoeffRing(P)) )
    {
      CoCoA_ASSERT_ALWAYS(IsDivisible(f1f2, f1));
      CoCoA_ASSERT_ALWAYS(IsDivisible(f1f2, f2));
      CoCoA_ASSERT_ALWAYS(!IsDivisible(f1, f2));
      CoCoA_ASSERT_ALWAYS(!IsDivisible(f2, f1));
      CoCoA_ASSERT_ALWAYS(f1f2/f1 == f2);
      CoCoA_ASSERT_ALWAYS(f1f2/f2 == f1);
      CoCoA_ASSERT_ALWAYS((power(f1f2,2))/(f1f2*f1f2) == 1);
    }
    CoCoA_ASSERT_ALWAYS(power(f1,2) == f1*f1);
    CoCoA_ASSERT_ALWAYS(power(f2,2) == f2*f2);    
    CoCoA_ASSERT_ALWAYS(deriv(f1f2, x[0]) == deriv(f1, x[0])*f2 + f1*deriv(f2, x[0]));
    CoCoA_ASSERT_ALWAYS(deriv(x[1]+2*x[0], x[0]) == 2);

    //  CoCoA_ASSERT_ALWAYS(gcd(f1*f1f2, f2*f1f2) == f1f2);

    CoCoA_ASSERT_ALWAYS(CmpWDeg(f1, one(P)) > 0);
    CoCoA_ASSERT_ALWAYS(CmpWDeg(f2, one(P)) > 0);
    CoCoA_ASSERT_ALWAYS(CmpWDeg(one(P), f1) < 0);
    CoCoA_ASSERT_ALWAYS(CmpWDeg(one(P), f2) < 0);
    CoCoA_ASSERT_ALWAYS(CmpWDeg(f1f2, f1) > 0);
    CoCoA_ASSERT_ALWAYS(CmpWDeg(f1f2, f2) > 0);
    CoCoA_ASSERT_ALWAYS(CmpWDeg(f1, f1f2) < 0);
    CoCoA_ASSERT_ALWAYS(CmpWDeg(f2, f1f2) < 0);

    CoCoA_ASSERT_ALWAYS(CmpWDeg(power(f1f2, 3), f1f2) > 0);

    f2f1 = 1;
    CoCoA_ASSERT_ALWAYS(IsOne(f2f1));

    CoCoA_ASSERT_ALWAYS(wdeg(f1) + wdeg(one(P)) == wdeg(f1));
    CoCoA_ASSERT_ALWAYS(wdeg(f2) + wdeg(one(P)) == wdeg(f2));
    CoCoA_ASSERT_ALWAYS(wdeg(f1) + wdeg(f2) == wdeg(f1f2));

    vector<long> exps0(NumIndets(P));  exps0[0] = 1;
    vector<long> exps1(NumIndets(P));  exps1[1] = 1;
    RingElem g0(x[0]), g1(x[1]);
    PushBack(g0, -one(CoeffRing(P)), exps1);
    PushFront(g1, one(CoeffRing(P)), exps0);
    CoCoA_ASSERT_ALWAYS(g0 == x[0]-x[1]);
    CoCoA_ASSERT_ALWAYS(g1 == x[0]+x[1]);

    cout << "------------------------------------------------" << endl << endl;
  }


  vector<long> intersection(const vector<long>& L1, const vector<long>& L2)
  {
    vector<long> ans;
    auto it1 = L1.begin();
    auto it2 = L2.begin();
    while (it1 != L1.end() && it2 != L2.end())
    {
      if (*it1 < *it2) { ++it1; continue; }
      if (*it1 > *it2) { ++it2; continue; }
      ans.push_back(*it1);
      ++it1; ++it2;
    }
    return ans;
  }

  void test_redmine1641()
  {
    SparsePolyRing P = NewPolyRing(NewZZmod(19), symbols("x,y,z"));
    RingElem f1(P, "x+1");
    RingElem f2(P, "y+2");
    RingElem f3(P, "z+3");
    RingElem f4(P, "x^2+y^2+4");
    RingElem f5(P, "x^2+z^2+5");
    RingElem f6(P, "y*z+6");
    RingElem f7(P, "x^2+y^2+z^2-x*y*z");
//     const vector<RingElem> L{f1,f2,f3,f4,f5,f6,f7}; // TOO SLOW!
    const vector<RingElem> L{f1,f2,f6,f7};
    for (auto it1=SubsetIter(len(L)); !IsEnded(it1); ++it1)
    {
      RingElem f = one(P);
      for (long j: *it1)  f *= L[j];
      for (auto it2=SubsetIter(len(L)); !IsEnded(it2); ++it2)
      {
        RingElem g = one(P);
        for (long j: *it2)  g *= L[j];
        RingElem h = gcd(f,g);
        RingElem k = one(P);
        for (long i: intersection(*it1, *it2))
          k *= L[i];
        CoCoA_ASSERT_ALWAYS(IsConstant(h/k));
        
      }
    }
  }


  void program()
  {
    GlobalManager CoCoAFoundations(UseSymmResidues);
    const vector<RingElem> v(3, one(RingZZ()));
    const ConstMatrix I = IdentityMat(RingZZ(), 3);
    PPOrdering Deg2(NewMatrixOrdering(MakeTermOrdMat(ConcatVer(RowMat(v), I)), 2));

    const QuotientRing Fp = NewZZmod(101);
    const QuotientRing R = NewZZmod(10);
    const SparsePolyRing ZZxy = NewPolyRing(RingZZ(), symbols("x,y"));
    const SparsePolyRing Rx_DMP = NewPolyRing_DMP(R, SymbolRange("x",0,2));
    const SparsePolyRing Rx_DMPI = NewPolyRing_DMPI(R, SymbolRange("x",0,2));
    //    const SparsePolyRing Rx_DMPII = NewPolyRing_DMPII(R, 3);
    const SparsePolyRing Fpy = NewPolyRing_DMPI(Fp, SymbolRange("y",0,2));
    const SparsePolyRing FpyII = NewPolyRing_DMPII(Fp, SymbolRange("y",0,3));
    const SparsePolyRing Qx = NewPolyRing_DMPI(RingQQ(), SymbolRange("x",1,4));
    const SparsePolyRing Q2x = NewPolyRing_DMPI(RingQQ(), symbols("x,y,z"), Deg2);

    TestSparsePolyRing(ZZxy);
    TestSparsePolyRing(Rx_DMP);
    TestSparsePolyRing(Rx_DMPI);
    //    TestSparsePolyRing(Rx_DMPII); // cannot have DMPII for ZZ/(10)
    TestSparsePolyRing(Fpy);
    TestSparsePolyRing(FpyII);
    TestSparsePolyRing(Qx);
    TestSparsePolyRing(Q2x);
    test_redmine1641();
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
