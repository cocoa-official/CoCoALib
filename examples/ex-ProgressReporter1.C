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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-ProgressReporter1.C,v 1.9 2022/02/13 09:56:58 abbott Exp $
// $Log: ex-ProgressReporter1.C,v $
// Revision 1.9  2022/02/13 09:56:58  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.8  2015/11/24 12:46:04  abbott
// Summary: Renamed "isqrt" --> "FloorSqrt"
//
// Revision 1.7  2015/11/04 10:13:19  abbott
// Summary: Made example faster
//
// Revision 1.6  2015/06/29 15:42:11  bigatti
// *** empty log message ***
//
// Revision 1.5  2015/06/29 13:06:53  bigatti
// -- code in namespace CoCoA
//
// Revision 1.4  2015/06/25 14:18:29  abbott
// Summary: Improved comments/readability
// Author: JAA
//
// Revision 1.3  2014/07/28 14:43:55  abbott
// Summary: Reduced a parameter to make the example faster (was too slow on my netbook)
// Author: JAA
//
// Revision 1.2  2014/04/22 14:52:23  abbott
// Summary: Made slightly clearer
// Author: JAA
//
// Revision 1.1  2014/04/22 14:09:12  abbott
// Summary: Example for ProgressReporter
// Author: JAA
//
// Revision 1.7  2013/05/28 07:07:04  bigatti
// -- added "cout << boolalpha": useful for testing
//
// Revision 1.6  2012/11/30 14:04:55  abbott
// Increased visibility of comment saying "put your code here".
//
// Revision 1.5  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.4  2008/10/07 12:12:54  abbott
// Removed useless commented out #include.
//
// Revision 1.3  2007/05/31 16:06:16  bigatti
// -- removed previous unwanted checked-in version
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.9  2007/03/07 11:51:40  bigatti
// -- improved test alignment
//
// Revision 1.8  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.7  2007/03/02 17:46:40  bigatti
// -- unique RingZ and RingQ
// -- requires foundations.H ;  foundations blah;  (thik of a better name)
//
// Revision 1.6  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.5  2007/03/01 13:52:59  bigatti
// -- minor: fixed typo
//
// Revision 1.4  2007/02/28 15:15:56  bigatti
// -- minor: removed quotes in description
//
// Revision 1.3  2007/02/12 16:27:43  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.2  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.1  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
