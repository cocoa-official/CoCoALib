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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-bool3.C,v 1.7 2022/02/13 09:57:00 abbott Exp $
// $Log: ex-bool3.C,v $
// Revision 1.7  2022/02/13 09:57:00  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.6  2018/08/28 12:32:21  abbott
// Summary: Added useful comments; minor improvement to readability
//
// Revision 1.5  2015/06/29 15:51:44  bigatti
// -- code in namespace CoCoA
//
// Revision 1.4  2012/05/30 14:13:38  abbott
// Mildly simplified.  Approved by Alessandra!
//
// Revision 1.3  2012/05/30 12:39:03  bigatti
// -- changed IsPrimeFast --> IsPrime3
//
// Revision 1.2  2012/05/29 16:35:43  abbott
// Improved long description.
//
// Revision 1.1  2012/05/29 16:25:53  abbott
// New example prog for bool3.
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
