// Copyright (c) 2023  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

#include <csignal>
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example shows how to make a CoCoALib program watch for   \n"
  "signals, and if one is detected, it prints extra info.        \n";


const string LongDescription =
  "This example shows how to make a CoCoALib program watch for   \n"
  "signals, and if one is detected, it prints extra info.  Look  \n"
  "at ex-interrupt1.C and ex-interrupt2.C before looking at this.\n";

//----------------------------------------------------------------------

namespace CoCoA
{

  ///////////////////////////////////////////////////////
  // A long computation which checks occasionally for interrupts
  // If an interrupt was received, a progress message is printed, and
  // then an InterruptReceived exception is thrown.
  BigRat LongComputation()
  {
    BigRat prod(1,1);
    constexpr int MaxPrime = 1000000; // so it takes about 2.5s on my computer
    for (int p=2; p < MaxPrime; p = NextPrime(p))
    {
      // Instead of calling CheckForInterrupt(...)
      if (GetAndResetSignalReceived())
      { // an interrupt was received; so we print out some extra info then throw an exception
        clog << "\nDetected interrupt at p=" << p << endl;
        ThrowException(InterruptReceived("LongComputation"));
      }
      prod *= BigRat(p-1,p);
    }
    return prod;
  }


  void program()
  {
    cout << ShortDescription << endl;
    GlobalManager CoCoAFoundations;

    SignalWatcher MonitorSIGINT(SIGINT); // RAII object, IMMEDIATELY STARTS "watching" for SIGINT
    // -------------------------------------------------------
    // Now do a long computation which checks for interrupts...
    cout << "Starting interruptible computation..." << endl;
    const BigRat ans1 = LongComputation();
    cout << "Computation finished: ans = " << FloatStr(ans1) << endl;
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
    PrintInFrame(cerr, intr);
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
