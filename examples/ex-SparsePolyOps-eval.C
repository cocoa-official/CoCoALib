// Copyright (c) 2023  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example illustrates evaluating a univariate polynomial  \n"
  "with integer coefficients at an integer or rational point.   \n";

const string LongDescription =
  "This example illustrates evaluating a univariate polynomial  \n"
  "with integer coefficients at an integer or rational point.   \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;
    SignalWatcher MonitorInterrupt(SIGINT); // you must also call CheckForInterrupt every so often

    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    ring P = NewPolyRing(RingQQ(), symbols("x,y"));
    RingElem f(P, "x^2+x+41");  // NB coeffs are integer

    // Repeated evaluation of a single poly via an EvalUPoly object:
    EvalUPoly F(f);
    cout << "f = " << f << endl;
    cout << "f(1) = " << F(1) << endl;
    cout << "numerator f(1/2) = " << F(1,2) << endl;

    // Direct evaluation without (explicitly) creating an EvalUPoly object:
    // (these direct functions are suitable when evaluating at only very few points)
    cout << "f(-1) = " << EvalAt(f, -1) << endl;
    cout << "numerator of f(-1/2) = " << EvalAt(f, -1,2) << endl;
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

  catch (const CoCoA::InterruptReceived& intr)
  {
    cerr << endl
         << "------------------------------" << endl
         << ">>>  CoCoALib interrupted  <<<" << endl
         << "------------------------------" << endl
         << "-->>  " << intr << "  <<--" << endl;
    return 2;
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
