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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-interrupt1.C,v 1.10 2024/03/08 19:59:58 abbott Exp $
// $Log: ex-interrupt1.C,v $
// Revision 1.10  2024/03/08 19:59:58  abbott
// Summary: Used constexpr (redmine 1511)
//
// Revision 1.9  2022/02/13 09:57:00  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.8  2017/08/08 13:48:08  abbott
// Summary: Improved comment, and var name
//
// Revision 1.7  2017/07/22 12:52:40  abbott
// Summary: Updated exception types
//
// Revision 1.6  2017/07/08 19:07:48  abbott
// Summary: updated example for interrupt; added new example too.
//
// Revision 1.5  2016/11/18 18:06:06  abbott
// Summary: Made simpler by using new CoCoALib built-in "signal monitor"
//
// Revision 1.4  2015/11/04 10:12:41  abbott
// Summary: Made it faster
//
// Revision 1.3  2015/06/26 14:55:20  abbott
// Summary: Improved readability
// Author: JAA
//
// Revision 1.2  2015/06/25 14:49:51  abbott
// Summary: Added string arg to CheckForInterrupt
// Author: JAA
//
// Revision 1.1  2015/05/20 16:52:43  abbott
// Summary: New example illustrating checking and handling interrupts
// Author: JAA
//
//
