// Copyright (c) 2017  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example shows how to use a CpuTimeLimit object to limit        \n"
  "the CPU time used in a function.  Compare this example with that in \n"
  "ex-CpuTimeLimit2.C.";

const string LongDescription =
  "This example shows how to use a CpuTimeLimit object to limit the    \n"
  "CPU time used in a function containing a loop: just call the memfn  \n"
  "operator(), which will throw TimeoutException if timeout occurs.    \n"
  "Note that timeout may occur a bit later than requested.             \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  BigRat SumOfReciprocals(int n, const CpuTimeLimit& CheckForTimeout)
  {
    BigRat ans;
    for (int i=1; i <= n; ++i)
    {
      CheckForTimeout("SumOfReciprocals");
      ans += BigRat(1,i);
    }
    return ans;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    // Now call fn defined above: if TimeoutException is thrown, we
    // will catch it (and just print a message).
    try
    {
      const BigRat sum = SumOfReciprocals(50000, CpuTimeLimit(1.0, IterationVariability::low));
      cout << "sum = " << FloatStr(sum) << endl;
    }
    catch (const CoCoA::TimeoutException&)
    {
      // For this example, we just print a message saying it timed out.
      cout << "-------------------------------" << endl
           << ">>>  Computation timed out  <<<" << endl
           << "-------------------------------" << endl;
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
