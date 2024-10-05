//   Copyright (c)  2012  Anna Bigatti

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
#include "CoCoA/DenseUPolyRing.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingDenseUPolyClean.H"
#include "CoCoA/RingTwinFloat.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/VectorOps.H"
#include "CoCoA/symbol.H"


#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <vector>
using std::vector;


//----------------------------------------------------------------------
// Test for DenseUPolyRing functions on RingElem
// functions: +, -f, f-g, *, StdDeg, gcd, TidyGens, %, IsIrred
// environments: DenseUPolyClean (Fp, Z)
//----------------------------------------------------------------------
namespace CoCoA
{

  void TestDenseUPolyRing(DenseUPolyRing P)
  {
    //    cout << "TESTING: " << P << endl << endl;
    cout << "TESTING:" << endl;
    P->myOutputSelfLong(cout);
    cout  << endl << endl;

    if (!IsExact(P))
    {
      RingElem tf(NewRingTwinFloat(1024), 4);
      RingElem four(P);
      P->myRecvTwinFloat(raw(four), raw(tf));
      CoCoA_ASSERT_ALWAYS(four == 4);
    }

    const RingElem& x = indet(P, 0);
    RingElem f1 = 3*x +23,
      f2 = IndetPower(P, 0, 4) + x,  // x^4 + x
      f3 = x + 29,
      f4(P);

    P->myAssignNonZeroCoeff(raw(f1), raw(one(CoeffRing(P))), 2); // x^2

    cout << "Given f1 = " << f1 << endl;
    cout << "  and f2 = " << f2 << endl << endl;
    cout << "f1+f2   gives  " << f1+f2 << endl;
    cout << "f2+f1   gives  " << f2+f1 << endl;
    cout << "f2-f1   gives  " << f2-f1 << endl;
    cout << "-f1   gives  "   << -f1 << endl;
    cout << "deg((2*x+1)*(3*x+1)) gives  "   << deg((2*x+1)*(3*x+1)) << endl;
    //   //  cout << "gcd(f1,f2)   gives  " << gcd(f1,f2) << endl;
    cout << "StdDeg(f1)   gives  " << StdDeg(f1) << endl;
    cout << "StdDeg(f2)   gives  " << StdDeg(f2) << endl;

    // if (IsTrueGCDDomain(P)) CoCoA_ASSERT_ALWAYS(gcd(f1*f1f2, f2*f1f2) == f1f2);
    RingElem f1f2 = f1*f2;
    RingElem f2f1 = f2*f1;
    if (IsField(CoeffRing(P)) || IsTrueGCDDomain(CoeffRing(P)))
    {
      CoCoA_ASSERT_ALWAYS( gcd(f1,f2) == 1 );
      CoCoA_ASSERT_ALWAYS( IsInvertible(gcd(f1f2, f1*f3)/f1) );
      CoCoA_ASSERT_ALWAYS( IsInvertible(gcd(f1f2, f2*f3)/f2) );
      CoCoA_ASSERT_ALWAYS( IsInvertible(gcd(f1*f3, f2*f3)/f3) );
      CoCoA_ASSERT_ALWAYS( IsInvertible(gcd(f1*f1f2, f2f1*f2)/(f1f2)) );
      CoCoA_ASSERT_ALWAYS( IsDivisible(f1f2, f1) );
      CoCoA_ASSERT_ALWAYS( IsDivisible(f1f2, f2) );
      CoCoA_ASSERT_ALWAYS( !IsDivisible(f1, f2) );
      CoCoA_ASSERT_ALWAYS( !IsDivisible(f2, f1) );
      CoCoA_ASSERT_ALWAYS( (f1f2)/f1 == f2 );
      CoCoA_ASSERT_ALWAYS( (f1f2)/f2 == f1 );
      CoCoA_ASSERT_ALWAYS( deriv(f1,x) == 2*x+3 );
      try { std::cout << deriv(f1,x+1) << std::endl; }
      catch (const CoCoA::ErrorInfo& err) { if (err != ERR::ReqIndet) throw; }
    }
    if (IsField(CoeffRing(P)))
    {
      try
      {
        CoCoA_ASSERT_ALWAYS( IsIrred(f1) );
        CoCoA_ASSERT_ALWAYS( !IsIrred(f2) );
        CoCoA_ASSERT_ALWAYS( !IsIrred(f1f2) );
      }
      catch (const CoCoA::ErrorInfo& err) { if (err != ERR::NYI) throw; }
      RingElem g(P), acof(P), bcof(P);
      P->myExgcd(raw(g), raw(acof), raw(bcof), raw(f1), raw(f2));
      CoCoA_ASSERT_ALWAYS(IsInvertible(g) || (g == gcd(f1,f2)));
      CoCoA_ASSERT_ALWAYS(g == acof*f1 + bcof*f2);
    }
  
    P->myMulByXExp(raw(f1), 5);
    cout << "  P->myMulByXExp(raw(f1), 5) gives " << f1 << endl;
    P->myMulBy1MinusXExp(raw(f1), 2);
    cout << "  P->myMulBy1MinusXExp(raw(f1), 2) gives " << f1 << endl;
    P->myAddMulLM(raw(f1), raw(LC(x)), 2, raw(f2));
    cout << "  P->myAddMulLM(raw(f1), raw(LC(x)), 2, raw(f2)) gives " << f1 << endl;
  
    cout << "one(P) = " << one(P) << endl;

    f1f2 = f1*f2;  // previous lines changed f1
    f2f1 = f2*f1;  // previous lines changed f1
    CoCoA_ASSERT_ALWAYS(IsOne(one(P)));
    CoCoA_ASSERT_ALWAYS(!IsOne(f1));
    CoCoA_ASSERT_ALWAYS(!IsOne(f2));
    CoCoA_ASSERT_ALWAYS(P->myIsValid(raw(f1)));
    CoCoA_ASSERT_ALWAYS(P->myIsValid(raw(f2)));
    CoCoA_ASSERT_ALWAYS(P->myIsValid(raw(f3)));
    CoCoA_ASSERT_ALWAYS(P->myIsValid(raw(f4)));
    CoCoA_ASSERT_ALWAYS(P->myIsValid(raw(f1f2)));
    CoCoA_ASSERT_ALWAYS(P->myIsValid(raw(f2f1)));
    if (IsCommutative(P))  CoCoA_ASSERT_ALWAYS(f1f2 == f2f1);
    CoCoA_ASSERT_ALWAYS(f1 != f2);
    CoCoA_ASSERT_ALWAYS(f1+f2 == f2+f1);
    CoCoA_ASSERT_ALWAYS(f1+(-f1) == 0);
    CoCoA_ASSERT_ALWAYS(f1-f1 == 0);
    CoCoA_ASSERT_ALWAYS(f1f2 != f1);
    CoCoA_ASSERT_ALWAYS(f1f2 != f2);
  
    if (IsField(CoeffRing(P)))
    {
      ideal I = ideal(f1f2, f1*x, zero(P), f1*(x-1));
      cout << "TidyGens(I) = " << TidyGens(I) << endl;
      CoCoA_ASSERT_ALWAYS(!IsZero(I));
      CoCoA_ASSERT_ALWAYS(IsZero(f1 % I));
      CoCoA_ASSERT_ALWAYS(IsZero(f1f2 % I));
      CoCoA_ASSERT_ALWAYS(!IsZero((f1f2+x) % I));
    
      ideal J1 = ideal(f1f2);
      ideal J2 = ideal(f1*(3*x-7));
      CoCoA_ASSERT_ALWAYS(intersect(J1, J2) == ideal(f1f2*(3*x-7))); 
      CoCoA_ASSERT_ALWAYS(colon(J1, J2) == ideal(f2)); 
 
      vector<RingElem> h(3, zero(P));
      CoCoA_ASSERT_ALWAYS(TidyGens(ideal(h)).empty());
      CoCoA_ASSERT_ALWAYS(IsZero(ideal(h)));
    }
    if ( IsFractionFieldOfGCDDomain(CoeffRing(P)) )
    {
      CoCoA_ASSERT_ALWAYS( CommonDenom(x/3+one(P)/2) ==  6);
      CoCoA_ASSERT_ALWAYS( ClearDenom(x/3+one(P)/2) ==  2*x+3);
    }
    CoCoA_ASSERT_ALWAYS(IsZero(f4));
    f4 = 1;
    CoCoA_ASSERT_ALWAYS(IsOne(f4));
    f2f1 = 1;
    CoCoA_ASSERT_ALWAYS(P->myIsValid(raw(f2f1)));
    CoCoA_ASSERT_ALWAYS(IsOne(f2f1));

    cout << "------------------------------------------------" << endl;
  }


  void program()
  {
    GlobalManager CoCoAFoundations(UseNonNegResidues);

    const QuotientRing Fp = NewZZmod(101);
    const QuotientRing F6 = NewZZmod(6);
    const RingTwinFloat RR = NewRingTwinFloat(150);

    TestDenseUPolyRing(NewPolyRing_DUP(RingZZ()));
    TestDenseUPolyRing(NewPolyRing_DUP(Fp));
    TestDenseUPolyRing(NewPolyRing_DUP(F6, symbol("y")));
    TestDenseUPolyRing(NewPolyRing_DUP(RingQQ()));
    TestDenseUPolyRing(NewPolyRing_DUP(RR));
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
