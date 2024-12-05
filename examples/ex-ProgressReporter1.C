// Copyright (c) 2014-2015  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is a simple example showing how to use a ProgressReporter. \n";

const string LongDescription =
  "An example of how to use a ProgressReporter to print occasional \n"
  "updates during a long iterative computation.                    \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  bool IsSquare(long n)
  {
    const long sqrtn = FloorSqrt(n);
    return (n == sqrtn*sqrtn);
  }

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    cout << "-------------------------------------------------------" << endl;
    cout << "First loop:" << endl;

    ProgressReporter ProgressLoop1(1.0); // print reports roughly every 1.0 seconds
    for (long n=1; n < 77777; ++n)
    {
      ProgressLoop1(); // print progress count (at specified intervals)
      if (n%4 == 0 || n%4 == 3) continue;
      long NumReprs = 0;
      for (long j=1; 2*j*j <= n; ++j)
        if (IsSquare(n-j*j)) ++NumReprs;
      if (NumReprs > 8)
        cout << n << " has " << NumReprs << " representations of the form A^2+B^2" << endl;
    }

    cout << endl;
    cout << "-------------------------------------------------------" << endl;
    cout << "Second loop:" << endl;

    // On my computer this loop takes less than 1 sec, so ProgressReporter prints nothing...
    ProgressReporter ProgressLoop2(1.0); // print reports roughly every 1.0 seconds
    for (long p=101; ; p = NextPrime(p))
    {
      ProgressLoop2(p); // printed progress report will also give value of p
      if (PrimitiveRoot(p) < 64) continue;
      cout << "Prime " << p << " has least primitive root " << PrimitiveRoot(p) << endl;
      break;
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
