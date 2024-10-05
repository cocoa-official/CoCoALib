// Copyright (c) 2017  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example shows how to use a CpuTimeLimit object to limit the   \n"
  "CPU time used in a section of code, and what myReset   \n"
  "does.  Compare this example with ex-CpuTimeLimit1.C.";

const string LongDescription =
  "This example shows how to use a CpuTimeLimit object to limit the   \n"
  "CPU time used in a section of code.  In particular, it shows use of\n"
  "myReset between two loops (with differing costs for a  \n"
  "single iteration).  Compare this example with ex-CpuTimeLimit1.C.";

//----------------------------------------------------------------------

namespace CoCoA
{

  void AnotherLongComputation()
  {
    if (FloorLog2(factorial(5867281)) != 123456790)
      CoCoA_THROW_ERROR("Wrong answer!", "AnotherLongComputation");
  }

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    CpuTimeLimit CheckForTimeout(1.0);  // <--- start checking CPU usage from here
    // Calling CheckForTimeout() will throw TimeoutException when timeout occurs
    // this exception is then caught and handled in main (below).

    // First loop
    int CountPrimes = 0;
    for (int k=1; k < 4000; ++k)
    {
      CheckForTimeout("Loop 1");
      if (IsPrime(k)) ++CountPrimes;
    }
    cout << "Loop 1: CountPrimes=" << CountPrimes << endl;

    // Second loop
    CheckForTimeout.myReset(); // IMPORTANT to call myReset; try running without this call!
    CountPrimes = 0;
    for (int k=1; k < 4000; ++k)
    {
      CheckForTimeout("Loop 2");
      if (IsProbPrime((power(6,k)-1)/5)) ++CountPrimes;
    }
    cout << "Loop 2: CountPrimes=" << CountPrimes << endl;
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
  catch (const CoCoA::TimeoutException& exc)
  {
    // For this example we do not consider time-out as an error
    cout << endl
         << "-------------------------------" << endl
         << ">>>  Computation timed out  <<<" << endl
         << "-------------------------------" << endl;
      cout << endl << "PS timeout occurred in " << context(exc) << endl;
    return 0; // do not consider time-out as an error
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-CpuTimeLimit2.C,v 1.16 2024/06/09 20:10:19 abbott Exp $
// $Log: ex-CpuTimeLimit2.C,v $
// Revision 1.16  2024/06/09 20:10:19  abbott
// Summary: Renamed myPrepareForNewLoop to myReset (ProgressReporter & CpuTimeLimit); remine 1705
//
// Revision 1.15  2023/07/04 09:14:41  abbott
// Summary: Changed TimeOut to Timeout
//
// Revision 1.14  2022/02/13 09:56:56  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.13  2020/06/17 15:49:19  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.12  2019/09/16 14:38:14  abbott
// Summary: Added using namespace std
//
// Revision 1.11  2018/05/25 09:24:46  abbott
// Summary: Major redesign of CpuTimeLimit (many consequences)
//
// Revision 1.10  2018/03/12 12:42:26  abbott
// Summary: Now prints out on cout rather than clog when timeout occurs
//
// Revision 1.9  2017/10/16 15:25:32  abbott
// Summary: Time-out now announced on cout (instead of clog)
//
// Revision 1.8  2017/09/13 08:45:34  abbott
// Summary: Now exit status is 0 even wen time-out occurs.
//
// Revision 1.7  2017/09/06 13:55:44  abbott
// Summary: Added a comment
//
// Revision 1.6  2017/08/08 13:47:26  abbott
// Summary: Changed var name (limit --> timeout)
//
// Revision 1.5  2017/07/22 12:55:15  abbott
// Summary: Used PrintInFrame
//
// Revision 1.4  2017/07/21 15:05:32  abbott
// Summary: Added helpful(?) print statement
//
// Revision 1.3  2017/07/21 14:27:23  bigatti
// -- increased time limit ;-)
//
// Revision 1.2  2017/07/21 13:17:57  abbott
// Summary: Minor improvements
//
// Revision 1.1  2017/07/15 15:15:44  abbott
// Summary: Added examples for CpuTimeLimit
//
//
