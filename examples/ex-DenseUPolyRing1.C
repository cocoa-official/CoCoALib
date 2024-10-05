// Copyright (c) 2024  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "DenseUPolyRing functions on RingElem:  \n"
  "functions: +, -f, f-g, *, StdDeg, gcd, TidyGens.\n";

const string LongDescription =
  "This still has to be improved to become an interesting example (copied from test)"
  "Environments: DenseUPolyClean (Fp, Z)  \n";
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
      cout << "four == 4  gives " << (four == 4) << endl;
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
      cout << "gcd(f1,f2) == 1   gives  " << (gcd(f1,f2) == 1) << endl;
      cout << "IsInvertible(gcd(f1f2, f1*f3)/f1)   gives  " << IsInvertible(gcd(f1f2, f1*f3)/f1) << endl;
      cout << "IsDivisible(f1f2, f2)   gives  " << IsDivisible(f1f2, f2) << endl;
      cout << "!IsDivisible(f1, f2)   gives  " << !IsDivisible(f1, f2) << endl;
      cout << "(f1f2)/f1 == f2   gives  " << ((f1f2)/f1 == f2) << endl;
      cout << "deriv(f1,x)  gives  " << deriv(f1,x) << endl;
    }
    if (IsField(CoeffRing(P)))
    {
      try
      {
        cout << "IsIrred(f1)   gives  " << IsIrred(f1) << endl;
        cout << "IsIrred(f2)   gives  " << IsIrred(f2) << endl;
        cout << "IsIrred(f1f2)   gives  " << IsIrred(f1f2) << endl;
      }
      catch (const CoCoA::ErrorInfo& err) { if (err != ERR::NYI) throw; }
      RingElem g(P), acof(P), bcof(P);
      P->myExgcd(raw(g), raw(acof), raw(bcof), raw(f1), raw(f2));
      cout << "IsInvertible(g) || (g == gcd(f1,f2))  gives  " << (IsInvertible(g) || (g == gcd(f1,f2))) << endl;
      cout << "g == acof*f1 + bcof*f2  gives  " << (g == acof*f1 + bcof*f2) << endl;
    }
  
    P->myMulByXExp(raw(f1), 3);
    cout << "  P->myMulByXExp(raw(f1), 3) gives " << f1 << endl;
    P->myMulBy1MinusXExp(raw(f1), 2);
    cout << "  P->myMulBy1MinusXExp(raw(f1), 2) gives " << f1 << endl;
    P->myAddMulLM(raw(f1), raw(LC(x)), 2, raw(f2));
    cout << "  P->myAddMulLM(raw(f1), raw(LC(x)), 2, raw(f2)) gives " << f1 << endl;
  
    cout << "one(P) = " << one(P) << endl;

    f1f2 = f1*f2;  // previous lines changed f1
    f2f1 = f2*f1;  // previous lines changed f1
    cout << "IsOne(one(P))  gives  " << (IsOne(one(P))) << endl;
    cout << "IsOne(f1)  gives  " << (IsOne(f1)) << endl;
    cout << "P->myIsValid(raw(f1))  gives  " << (P->myIsValid(raw(f1))) << endl;
    if (IsCommutative(P))  cout << "f1f2 == f2f1  gives  " << (f1f2 == f2f1) << endl;
    cout << "f1 != f2  gives  " << (f1 != f2) << endl;
    cout << "f1+f2 == f2+f1  gives  " << (f1+f2 == f2+f1) << endl;
    cout << "f1+(-f1) == 0  gives  " << (f1+(-f1) == 0) << endl;
    cout << "f1-f1 == 0  gives  " << (f1-f1 == 0) << endl;
    cout << "f1f2 != f1  gives  " << (f1f2 != f1) << endl;
    cout << "f1f2 != f2  gives  " << (f1f2 != f2) << endl;
  
    if (IsField(CoeffRing(P)))
    {
      ideal I = ideal(f1f2, f1*x, zero(P), f1*(x-1));
      cout << "TidyGens(I) = " << TidyGens(I) << endl;
      cout << "!IsZero(I)  gives  " << (!IsZero(I)) << endl;
      cout << "IsZero(f1 % I)  gives  " << (IsZero(f1 % I)) << endl;
      cout << "IsZero(f1f2 % I)  gives  " << (IsZero(f1f2 % I)) << endl;
      cout << "!IsZero((f1f2+x) % I)  gives  " << (!IsZero((f1f2+x) % I)) << endl;
    
      ideal J1 = ideal(f1f2);
      ideal J2 = ideal(f1*(3*x-7));
      cout << "intersect(J1, J2) == ideal(f1f2*(3*x-7)));  gives  " << (intersect(J1, J2) == ideal(f1f2*(3*x-7))) << endl;                                                                                                             cout << "colon(J1, J2) == ideal(f2));  gives  " << (colon(J1, J2) == ideal(f2)) << endl;
 
      vector<RingElem> h(3, zero(P));
      cout << "TidyGens(ideal(h)).empty()  gives  " << (TidyGens(ideal(h)).empty()) << endl;
      cout << "IsZero(ideal(h))  gives  " << (IsZero(ideal(h))) << endl;
    }
    if ( IsFractionFieldOfGCDDomain(CoeffRing(P)) )
    {
      cout << "CommonDenom(x/3+one(P)/2) ==  6  gives  " << (CommonDenom(x/3+one(P)/2) ==  6) << endl;
      cout << "ClearDenom(x/3+one(P)/2) ==  2*x+3  gives  " << (ClearDenom(x/3+one(P)/2) ==  2*x+3) << endl;
    }
    cout << "IsZero(f4)  gives  " << (IsZero(f4)) << endl;
    f4 = 1;
    cout << "IsOne(f4)  gives  " << (IsOne(f4)) << endl;
    f2f1 = 1;
    cout << "P->myIsValid(raw(f2f1))  gives  " << (P->myIsValid(raw(f2f1))) << endl;
    cout << "IsOne(f2f1)  gives  " << (IsOne(f2f1)) << endl;

    RingHom phi = PolyAlgebraHom(P, CoeffRing(P), "1");
    cout << "phi is " << phi << endl;
    cout << "phi("<<f1<<")  gives  " << phi(f1) << endl;
    cout << "phi("<<f2<<")  gives  " << phi(f2) << endl;
    
    
    cout << "------------------------------------------------" << endl;
  }


  void program()
  {
    GlobalManager CoCoAFoundations(UseNonNegResidues);
    SignalWatcher MonitorInterrupt(SIGINT); // you must also call CheckForInterrupt every so often

    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false
    
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
