// Copyright (c) 2016  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

#include <csignal>
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This short example shows an easy way of making CoCoALib programs  \n"
  "handle signals (or interrupts), by converting them to exceptions. \n";

const string LongDescription =
  "This short example shows an easy way of making CoCoALib programs  \n"
  "handle signals (or interrupts), by converting them to exceptions. \n"
  "There are two crucial parts: create a SignalWatcher to say which  \n"
  "signals to watch for, and call CheckForInterrupt when it is       \n"
  "convenient to act upon the interrupting signal.                   \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  //////////////////////////////////////////////////////////////////
  // A long computation which checks occasionally for interrupts.
  // Relevant details: main loop calls CheckForInterrupt frequently.
  void LongComputation()
  {
    BigRat prod(1,1);
    constexpr int MaxPrime = 700000; // so it takes about 4s on my computer
    for (int p=2; p < MaxPrime; p = NextPrime(p))
    {
      CheckForInterrupt("LongComputation"); // arg gives context info
      prod *= BigRat(p-1,p);
    }
    cout << "Product of (1-1/p) for primes p up to " << MaxPrime
         << " is about " << FloatStr(prod) << endl;
  }


  void program()
  {
    cout << ShortDescription << endl;
    GlobalManager CoCoAFoundations;

    SignalWatcher MonitorSIGINT(SIGINT); // RAII object
    // Now do a long computation which checks for interrupts...
    // Call it inside a try..catch block so any InterruptedBySignal
    // exception can be handled appropriately.
    try
    {
      LongComputation();
    }
    catch (const InterruptedBySignal& intr)
    {
      // LongComputation was interrupted by a signal.
      // Here we "handle" it by just printing out a message.
      PrintInFrame(cout, intr);
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
