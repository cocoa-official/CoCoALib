// Copyright (c) 2017  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

#include <csignal>
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example shows how to make a CoCoALib programs watch for   \n"
  "signals only during a section of code.  Look at ex-interrupt1  \n"
  "before looking at this example!  Compare what happens when you \n"
  "interrupt the first computation, and when you interrupt the    \n"
  "second (identical) computation.                                \n";

const string LongDescription =
  "This example shows how to make a CoCoALib programs watch for   \n"
  "signals only during a section of code.  Look at ex-interrupt1  \n"
  "before looking at this example!  Compare what happens when you \n"
  "interrupt the first computation, and when you interrupt the    \n"
  "second (identical) computation.                                \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  ///////////////////////////////////////////////////////
  // A long computation which checks occasionally for interrupts
  BigRat LongComputation()
  {
    BigRat prod(1,1);
    constexpr int MaxPrime = 800000; // so it takes about 2.5s on my computer
    for (int p=2; p < MaxPrime; p = NextPrime(p))
    {
      CheckForInterrupt("LongComputation"); // arg gives context info
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
    cout << "Starting first computation...[CoCoALib interrupt handler]" << endl;
    const BigRat ans1 = LongComputation();
    cout << "First computation finished: ans = " << FloatStr(ans1) << endl;

    deactivate(MonitorSIGINT); // now STOP "watching" for SIGINT
    cout << endl
         << "The SignalWatcher is now: " << MonitorSIGINT << endl
         << endl;

    // -------------------------------------------------------
    // Repeat the long computation, but without CoCoALib watching for signals.

    // Since the SignalWatcher is inactive, it does not set the global "flag"
    // when a signal is received, so the calls to CheckForInterrupt inside
    // the function LongComputation no longer "see" an interrupt; instead
    // the system interrupt handler will react to the signal.
    cout << "Starting second computation...[system interrupt handler]" << endl;
    const BigRat ans2 = LongComputation();
    cout << "Second computation finished: ans = " << FloatStr(ans2) << endl;
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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-interrupt2.C,v 1.8 2024/03/08 19:59:58 abbott Exp $
// $Log: ex-interrupt2.C,v $
// Revision 1.8  2024/03/08 19:59:58  abbott
// Summary: Used constexpr (redmine 1511)
//
// Revision 1.7  2024/01/21 19:52:07  abbott
// Summary: Added clarification
//
// Revision 1.6  2022/02/13 09:57:01  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.5  2017/08/08 13:48:48  abbott
// Summary: Improved var name and a comment
//
// Revision 1.4  2017/07/22 16:10:53  abbott
// Summary: doc for SignalWatcher
//
// Revision 1.3  2017/07/22 12:53:11  abbott
// Summary: Minor improvements
//
// Revision 1.2  2017/07/14 14:02:24  abbott
// Summary: Removed cruft
//
// Revision 1.1  2017/07/08 19:07:48  abbott
// Summary: updated example for interrupt; added new example too.
//
//
