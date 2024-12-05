// Copyright (c) 2017  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program illustrates use of the RootBound function for      \n"
  "obtaining an upper bound on the abs value of the complex roots. \n";

const string LongDescription =
  "This program illustrates use of the RootBound function for      \n"
  "obtaining an upper bound on the abs value of the complex roots. \n"
  "The function has an optional second arg.  The second arg says   \n"
  "how many Graeffe iterations to perform; these lead to a better  \n"
  "bound, but can be slow for large polynomials (as the example    \n"
  "makes clear).                                                   \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    ring P = NewPolyRing(RingQQ(), symbols("x,y"));
    RingElem f= RingElem(P, "x^2-2");
    BigRat B = RootBound(f);

    //-------------------------------------------------------
    // Now we try a larger polynomial g = product((x+k*k*k), k=1,...,200)
    // For such a big poly "RootBound" does just a few Graeffe steps
    // but we can use the second arg to force it to make more Graeffe steps
    // (though this will take extra time).
    RingElem g = one(P);
    const RingElem& x = indet(P,0);
    for (int k=0; k <= 250; ++k)
      g *= (x+k*k*k);

    cout << "Let g = product((x+k*k*k), k=1,...,200)" << endl
         << "We know the true root bound for g is 15625000." << endl << endl;

    double start = CpuTime();
    BigRat Bg = RootBound(g);
    const double t = CpuTime() - start;
    cout << "RootBound uses a heuristic to decide how many Graffe iterations to use." << endl
         << "RootBound(g) = " << Bg << "   time to compute: " << t << endl;    

    cout << endl << "We can specify how many Graeffe iterations to perform as a second arg." << endl
         << "We obtain a trade-off between speed and tightness:" << endl;
    for (int NumGraeffe = 0; NumGraeffe < 5; ++NumGraeffe)
    {
      start = CpuTime();
      Bg = RootBound(g,NumGraeffe);
      const double t1 = CpuTime() - start;
      cout << "RootBound(g," << NumGraeffe << ") = " << Bg << "  time to compute: " << t1 << endl;
    }

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
