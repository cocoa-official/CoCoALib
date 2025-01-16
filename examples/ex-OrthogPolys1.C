// Copyright (c) 2017  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program illustrates calling the functions for making  \n"
  "orthogonal polynomials (Chebyshev, Hermite, Laguerre).     \n";

const string LongDescription =
  "This program illustrates calling the functions for making  \n"
  "orthogonal polynomials (Chebyshev, Hermite, Laguerre).     \n";


//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    const ring P = NewPolyRing(RingQQ(), symbols("x"));
    const RingElem& x = indet(P,0);

    cout << "ChebyshevPoly(3,x), 1st kind: " << ChebyshevPoly(3,x) << endl;
    cout << "ChebyshevPoly2(3,x), 2nd kind: " << ChebyshevPoly2(3,x) << endl;
    cout << endl;
    cout << "HermitePoly(3,x) as used in physics: " << HermitePoly(3,x) << endl;
    cout << "HermitePoly2(3,x) as used in probability: " << HermitePoly2(3,x) << endl;
    cout << endl;
    cout << "LaguerrePoly(3,x) (mult by factorial(3)): " << LaguerrePoly(3,x) << endl;

    // We can also evaluate the polynomials by giving a second arg which
    // is a number rather than an indet:

    cout << endl;
    const RingElem half(P, BigRat(1,2));
    cout << "Value of ChebyshevPoly(3,x) at x=1/2: " << ChebyshevPoly(3, half) << endl;
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
  // catch (const CoCoA::InterruptReceived& intr)
  // {
  //   cerr << endl
  //        << "------------------------------" << endl
  //        << ">>>  CoCoALib interrupted  <<<" << endl
  //        << "------------------------------" << endl
  //        << "-->>  " << intr << "  <<--" << endl;
  //   return 2;
  // }
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
