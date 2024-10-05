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
#include "CoCoA/GlobalManager.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/NumTheory-gcd.H"
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/SparsePolyOps-RealRadical.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/SparsePolyOps-SturmSeq.H"
#include "CoCoA/PolyFamilies.H"
#include "CoCoA/error.H"
#include "CoCoA/VectorOps.H"
#include "CoCoA/symbol.H"

#include <algorithm>
using std::min;
#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
#include <vector>
using std::vector;


namespace CoCoA
{

  void test_coeffs(ConstRefRingElem f, ConstRefRingElem x, ConstRefRingElem y)
  {
    // Ripped off from ex-PolyRing3 :-)

    cout << "In the following we have f = " << f << endl;

    // Accessing coeffs via SparsePolyIter:
    cout << "Using a SparsePolyIter we decompose f as follows:" << endl;
    for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
    {
      cout << "coeff = " << coeff(it) << "  in ring " << owner(coeff(it)) << endl;
      cout << "PP    = " << PP(it)    << "  in " << owner(PP(it)) << endl;
      cout << endl;
    }

    cout << endl;
    cout << "-------------------------------------------------------" << endl;
    cout << endl;

    // Regard f as a poly in just "x" or just "y" we obtain:
    cout << "Coefficients with respect to certain indeterminates:" << endl << endl;
    cout << "Case (1) considering f as a *univariate* polynomial in..." << endl;
    cout << "...the indeterminate x the coeffs are " << CoeffVecWRT(f, x) << endl;
    cout << "...the indeterminate y the coeffs are " << CoeffVecWRT(f, y) << endl;
    cout << endl;

    cout << "Case (2) considering f as a sparse multivariate polynomial in..." << endl;
    cout << "...the indet x its structure is " << CoefficientsWRT(f, x) << endl;
    cout << "...the indet y its structure is " << CoefficientsWRT(f, y) << endl;
    vector<long> XandY; XandY.push_back(0); XandY.push_back(1);
    cout << "...the indets x & y its structure is " << CoefficientsWRT(f, XandY) << endl;
    cout << endl;
  }


  void test_CoeffVecWRTSupport(ConstRefRingElem f)
  {
    const vector<RingElem> C = CoeffVecWRTSupport(f,f);
    // C is in "reverse order": C[0] is coeff of smallest PP in f, etc.
    long i=NumTerms(f);
    for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it)
    {
      --i;
      CoCoA_ASSERT_ALWAYS(C[i] == coeff(it));
    }
  }

  
  void test_content(ConstRefRingElem f, ConstRefRingElem x, ConstRefRingElem y)
  {
    cout << "-------------------------------------------------------" << endl;
    cout << "The content of a polynomial" << endl << endl;

    // Content of f
    RingElem ContF = content(f);
    cout << "The \"numerical\" content of f is " << ContF << "  -- an element of ring " << owner(ContF) << endl;
    cout << endl;
    RingElem ContWRTx = ContentWRT(f, x);
    cout << "Content WRT x is " << ContWRTx << "  -- element of " << owner(ContWRTx) << endl;

    RingElem ContWRTy = ContentWRT(f, y);
    cout << "Content WRT y is " << ContWRTy << "  -- element of " << owner(ContWRTy) << endl;

    vector<long> XandY; XandY.push_back(0); XandY.push_back(1);
    RingElem ContWRTxy = ContentWRT(f, XandY);
    cout << "Content WRT x & y is " << ContWRTxy << "  -- element of " << owner(ContWRTxy) << endl;
  }


  //----------------------------------------------------------------------
  // This test shows how to split a polynomial into its homogeneous components.
  // It uses the function CutLF to "cut off" the leading form of a polynomial.
  //----------------------------------------------------------------------
  void test_CutLF(ConstRefRingElem f)
  {
    const ring& P = owner(f);
    // Use CutLF to split g = f^4 apart ino its homog components:
    RingElem g = power(f,4);
    long LastDeg = deg(g)+1;
    RingElem SumLF = zero(P); // will contain successive homog parts
    while (!IsZero(g))
    {
      const RingElem LF = CutLF(g); // changes g!!!
      CoCoA_ASSERT_ALWAYS(deg(LF) < LastDeg);
      CoCoA_ASSERT_ALWAYS(IsHomog(LF));
      SumLF += LF;
      LastDeg = deg(LF);
    }
    CoCoA_ASSERT_ALWAYS(SumLF == power(f,4));
  }


  void test_FixedDivisor(ConstRefRingElem x)
  {
    RingElem f = x*x+x;
    CoCoA_ASSERT_ALWAYS(FixedDivisor(f) == 2);

    f = x*(x-1)*(x-2)*(x-3)*(x-4);
    for (int c=0; c < 120; ++c)
      CoCoA_ASSERT_ALWAYS(FixedDivisor(f+c) == gcd(c,120));
  }


  //----------------------------------------------------------------------
  // This test checks that radical of a polynomial produces the expected result.
  // Will be called with coeff fields QQ or ZZ/(2) or ZZ/(32003)
  //----------------------------------------------------------------------
  void test_radical(const ring& k)
  {
    ring P = NewPolyRing(k, symbols("x,y,z"));
    const RingElem f1 = RingElem(P, "x");
    const RingElem f2 = RingElem(P, "3*y");
    const RingElem f3 = RingElem(P, "5*x+7");
    const RingElem f4 = RingElem(P, "11*x+13*z");
    const RingElem f5 = RingElem(P, "17*x+19*y+23*z+29");

    // First test just on monomials
    for (int e1 = 0; e1 < 4; ++e1)    
    for (int e2 = 0; e2 < 4; ++e2)    
    {
      const RingElem f = power(f1,e1) * power(f2,e2);
      const RingElem expected = power(f1,min(e1,1)) * power(f2,min(e2,1));
      const RingElem radF = radical(f);
      CoCoA_ASSERT_ALWAYS(monic(radF) == monic(expected));
      CoCoA_ASSERT_ALWAYS(monic(radical(-f)) == monic(radF));
      CoCoA_ASSERT_ALWAYS(monic(radical(radF)) == monic(radF));
    }

    // More general test on polynomials
    for (int e1 = 0; e1 < 1; ++e1)    // effectively ignore f1 to make the test faster.
    for (int e2 = 0; e2 < 4; ++e2)    
    for (int e3 = 0; e3 < 4; ++e3)    
    for (int e4 = 0; e4 < 4; ++e4)    
    for (int e5 = 0; e5 < 3; ++e5)    
    {
      const RingElem f = power(f1,e1) * power(f2,e2) * power(f3,e3) * power(f4,e4) * power(f5,e5);
      const RingElem expected = power(f1,min(e1,1)) * power(f2,min(e2,1)) * power(f3,min(e3,1)) * power(f4,min(e4,1)) * power(f5,min(e5,1));
      const RingElem radF = radical(f);
      CoCoA_ASSERT_ALWAYS(monic(radF) == monic(expected));
      CoCoA_ASSERT_ALWAYS(monic(radical(-f)) == monic(radF));
      CoCoA_ASSERT_ALWAYS(monic(radical(radF)) == monic(radF));
    }
  }



  //----------------------------------------------------------------------
  // This proc test checks the function RealRadical in some simple cases.
  // It also checks HasRealRoots3.
  //----------------------------------------------------------------------
  void test_RealRadical()
  {
    ring P = NewPolyRing(RingQQ(), symbols("x,y,z"));
    RingElem f1(P, "(x^2+1)*(x^4+4)*(x^5-5)");
    CoCoA_ASSERT_ALWAYS(RealRadical(f1) == RingElem(P, "x^5-5"));
    CoCoA_ASSERT_ALWAYS(RealRadical(f1*f1) == RingElem(P, "x^5-5"));

    RingElem f2(P, "x^2+y^2+1");
    CoCoA_ASSERT_ALWAYS(RealRadical(f2) == one(P));
    CoCoA_ASSERT_ALWAYS(RealRadical(f2*f2) == one(P));
    CoCoA_ASSERT_ALWAYS(RealRadical(f1*f2) == RingElem(P, "x^5-5"));

    RingElem f3(P, "x^2+y^2");
    CoCoA_ASSERT_ALWAYS(RealRadical(f3) == f3);
    CoCoA_ASSERT_ALWAYS(RealRadical(f3*f3) == f3);
    CoCoA_ASSERT_ALWAYS(RealRadical(f2*f3) == f3);
    CoCoA_ASSERT_ALWAYS(RealRadical(f1*f3) == f3*RingElem(P,"x^5-5"));
    CoCoA_ASSERT_ALWAYS(RealRadical(f1*f2*f3) == f3*RingElem(P,"x^5-5"));

    RingElem g1(P, "x^2+y^4+z^6");
    CoCoA_ASSERT_ALWAYS(IsTrue3(HasRealRoot3(g1)));
    CoCoA_ASSERT_ALWAYS(IsTrue3(HasRealRoot3(g1-1)));
    CoCoA_ASSERT_ALWAYS(IsFalse3(HasRealRoot3(g1+1)));
  }


  //----------------------------------------------------------------------
  // This test checks that NumRealRoots gives the correct answer;
  // internally NumRealRoots calls SturmSeq, so hopefully this also
  // also checks that SturmSeq is correct.
  //----------------------------------------------------------------------
  void test_NumRealRoots()
  {
    const ring P = NewPolyRing(RingQQ(), symbols("x"));
    RingElem f(P, "x^4+x^3-x-1");
    CoCoA_ASSERT_ALWAYS(NumRealRoots(f) == 2);
    
    RingElem g = ChebyshevPoly(50, RingElem(P,"x"));
    CoCoA_ASSERT_ALWAYS(NumRealRoots(g) == 50);
  }
  


  void program()
  {
    GlobalManager CoCoAFoundations;

    ring ZZ = RingZZ();
    SparsePolyRing ZZxy = NewPolyRing(ZZ, symbols("x,y"));
    const RingElem x = indet(ZZxy,0);
    const RingElem y = indet(ZZxy,1);
    const RingElem f = 2*x*x*y - 4*y*y + 6*x*x + 36;

    test_coeffs(f,x,y);
    test_CoeffVecWRTSupport(f);
    test_CoeffVecWRTSupport(power(f,3));
    test_content(f,x,y);
    test_CutLF(RingElem(ZZxy, "4*x^4+2*x^2*y^2-3*y^3-1"));
    test_FixedDivisor(y);
    test_radical(RingQQ());
    test_radical(NewZZmod(2)); 
    test_radical(NewZZmod(32003));
    test_RealRadical();
    test_NumRealRoots();
  }

}  // end of namespace CoCoA


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
