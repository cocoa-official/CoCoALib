// Copyright (c) 2005-2010  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Example showing how to create some simple polynomial rings.             \n"
  "It also shows some of the operations specific to elements of PolyRings. \n"
  "See also ex-ring2.C                                                     \n";

const string LongDescription =
  "This example program exhibits two things: various ways of creating \n"
  "polynomial rings, and several operations specific to elements of a \n"
  "polynomial ring.                                                   \n"
  "In the procedure \"program\" there are examples of creating        \n"
  "polynomial rings.                                                  \n"
  "In the procedure \"SomeComputations\" there are brief examples of  \n"
  "the operations specific to elements of PolyRings (e.g. deg).       \n"
  "See also ex-ring2.C.                                               \n";
//-----------------------------------------------------------------------

namespace CoCoA
{

  void SomeComputations(const PolyRing& P)
  {
    const long n = NumIndets(P);
    const vector<RingElem>& x = indets(P);

    // Put f = x[0] + x[1] + .. x[n-1]
    RingElem f(sum(indets(P)));

    // Put g = x[0] + 2*x[1] + ... + n*x[n-1] + 1234
    RingElem g(P, 1234);
    for (long i=0; i < n; ++i)  g += (i+1)*x[i];

    cout << "Some computations in the ring " << P << endl;
    cout << "Let f = " << f << endl
         << "and g = " << g << endl << endl;

    cout << "NumTerms(f) = " << NumTerms(f)
         << "   and NumTerms(g) = " << NumTerms(g) << endl;
    cout << "f+g = " << f+g << endl;
    cout << "f-g = " << f-g << endl;
    cout << "f*g = " << f*g << endl;
    if (IsDivisible(f, g))  cout << "f/g = " << f/g << endl;

    cout << "deg(f) = "  << deg(f) << endl;      // of type long
    cout << "StdDeg(f) = " << StdDeg(f) << endl; // same as deg(f)
    cout << "wdeg(f) = " << wdeg(f) << endl;     // of type CoCoA::degree
    cout << "deg(f, 0) = " << deg(f, 0) << endl; // of type long
    cout << "CmpWDeg(f, g) = " << CmpWDeg(f, g) << endl;
    cout << "CmpWDeg(power(f,2), g) = " << CmpWDeg(power(f,2), g) << endl;
    cout << "CmpWDeg(f, power(g,2)) = " << CmpWDeg(f, power(g,2)) << endl;
    if (!IsZero(f-g))
      cout << "LPP(f-g) = " << LPP(f-g) << "   in " << owner(LPP(f-g)) << endl;
    if (!IsZero(f+g))
      cout << "LC(f+g) = " << LC(f+g) << "   in " << owner(LC(f+g)) << endl;
    cout << "deriv(g, n-1) = " << deriv(g, n-1) << endl; // note second arg!
    if (GradingDim(owner(f)) > 0)
      cout << "IsHomog(f) = " << IsHomog(f) << endl;
    else
      cout << "IsHomog: works only if the grading dimension is positive!" << endl;
    cout << "--------------------------------------------" << endl << endl;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;
    cout << boolalpha; // so that bools print out as true/false

    cout << ShortDescription << endl;

    // Define some coefficient rings.
    ring QQ   = RingQQ();
    ring FF_7 = NewZZmod(7);

    // Now we show some different ways of creating polynomial rings.

    // [1] If we want, we can specify the names of the indets:
    PolyRing P1 = NewPolyRing(FF_7, symbols("beta,Psi")); // FF_7[beta,Psi]
    SomeComputations(P1);

    // [2] We can a range of names for the indets:
    PolyRing P2 = NewPolyRing(FF_7, SymbolRange("x", 1,3)); // FF_7[x[1],x[2],x[3]]
    SomeComputations(P2);

    // [3] We can let CoCoA choose "new" names for the indets:
    PolyRing P3 = NewPolyRing(QQ, NewSymbols(5));         // indets have (strange) "new" names
    SomeComputations(P3);
  }

} // end of namespace CoCoA

//----------------------------------------------------------------------
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

  CoCoA::BuildInfo::PrintAll(cerr);
  return 1;
}
