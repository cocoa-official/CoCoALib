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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-CpuTimeLimit1.C,v 1.14 2023/07/04 09:14:41 abbott Exp $
// $Log: ex-CpuTimeLimit1.C,v $
// Revision 1.14  2023/07/04 09:14:41  abbott
// Summary: Changed TimeOut to Timeout
//
// Revision 1.13  2022/02/13 09:56:56  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.12  2021/02/22 20:08:11  abbott
// Summary: Updated after changing interface to CpuTimeLimit
//
// Revision 1.11  2021/01/12 13:25:56  abbott
// Summary: Added "variability" to one ctor call
//
// Revision 1.10  2019/09/16 14:37:53  abbott
// Summary: Added using namespace std
//
// Revision 1.9  2018/05/25 09:24:45  abbott
// Summary: Major redesign of CpuTimeLimit (many consequences)
//
// Revision 1.8  2017/10/16 15:25:32  abbott
// Summary: Time-out now announced on cout (instead of clog)
//
// Revision 1.7  2017/09/13 08:45:33  abbott
// Summary: Now exit status is 0 even wen time-out occurs.
//
// Revision 1.6  2017/09/06 13:55:24  abbott
// Summary: Updated comment: exception is now called TimeoutException
//
// Revision 1.5  2017/08/08 13:46:54  abbott
// Summary: Added comments
//
// Revision 1.4  2017/07/23 15:27:24  abbott
// Summary: Renamed demo fn to avoid name clash
//
// Revision 1.3  2017/07/22 12:55:15  abbott
// Summary: Used PrintInFrame
//
// Revision 1.2  2017/07/21 13:17:57  abbott
// Summary: Minor improvements
//
// Revision 1.1  2017/07/15 15:15:44  abbott
// Summary: Added examples for CpuTimeLimit
//
//
