// Copyright (c) 2017  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example shows simple use of RealRadical and HasRealRoot3.  \n";

const string LongDescription =
  "This example shows simple use of RealRadical and HasRealRoot3.  \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    ring P = NewPolyRing(RingQQ(), symbols("x,y,z"));
    const RingElem& x = indet(P,0);

    RingElem f1(P, "(x^2+1)*(x^4+4)*(x^5-5)");
    cout << "Let f1 = " << f1 << endl;
    cout << "(a multiple of) RealRadical of f1 is  " << RealRadical(f1) << endl;

    // HasRealRoot3 is reliable for univariate polynomials:
    cout << "Does f1 have a real root?  " << HasRealRoot3(f1) << endl;


    cout << endl;
    cout << "HasRealRoot3 is heuristic for multivariate polynomials:"  << endl;

    RingElem g(P, "x^2+y^2+1");
    cout << "Let g = " << g << endl;
    cout << "Does g have a real root?  " << HasRealRoot3(g)  << endl;
    cout << "Does g+x have a real root?     [no]  " << HasRealRoot3(g+x) << endl;
    cout << "Does g+2*x have a real root?  [yes]  " << HasRealRoot3(g+2*x) << endl;
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
