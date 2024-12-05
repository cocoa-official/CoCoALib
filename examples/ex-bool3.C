// Copyright (c) 2012  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program shows a simple use of three-state booleans. \n";

const string LongDescription =
  "This program shows a simple use of three-state booleans.        \n"
  "We define a quick primality test which is guaranteed to be fast \n"
  "but which sometimes has to return a verdict of \"Don't know\" so\n"
  "it can keep its guarantee of speed.  We then see how often the  \n"
  "quick test gives a definite answer, and how often it has to say \n"
  "\"Don't know\".                                                 \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  // This fn is guaranteed to be fast; the price you must pay for
  // this speed is the possibility of a "Don't know" response.
  bool3 IsPrime3(long n)
  {
    if (n == 1) return false3;
    if (n == 2 || n == 3 || n == 5 || n == 7) return true3;
    if (n%2 == 0 || n%3 == 0 || n%5 == 0 || n%7 == 0) return false3;
    if (n < 121) return true3;
    return uncertain3;
  }

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    // We shall count how many numbers are classified as "prime",
    // "NotPrime" and "undecided" in these 3 variables:
    long yes = 0;
    long no = 0;
    long DontKnow = 0;

    for (long n=1; n <= 200; ++n)
    {
      const bool3 b = IsPrime3(n);
      if (IsTrue3(b)) ++yes;
      if (IsFalse3(b)) ++no;
      if (IsUncertain3(b)) ++DontKnow;
    }

    // Print out the results of the survey:
    cout << "Quick analysis of the primeness of the numbers 1..200:\n";
    cout << "  How many are surely prime?      " << yes      << endl;
    cout << "  How many are surely composite?  " << no       << endl;
    cout << "  How many need more analysis?    " << DontKnow << endl;
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
