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
