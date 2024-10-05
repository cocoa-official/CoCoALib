// Copyright (c) 2005,2012  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "When computing over a finite field normally it is best to use the\n"
  "function NewZZmod to create the field.  However, for the curious \n"
  "it is possible to create small prime finite fields stating the   \n"
  "particular implemention method (there are 3 possibilities).      \n";

const string LongDescription =
  "This program compares the speeds of computing sums and products in   \n"
  "the three different implementations of small prime finite fields.    \n"
  "It performs the timing tests for different sizes of prime, and       \n"
  "illustrates that different implementations have different upper      \n"
  "limits for the characteristic -- these limits depend on the platform.\n"
  "It also gives the performance of the finite field created by NewZZmod\n";

//----------------------------------------------------------------------

namespace CoCoA
{

  // Measure how long it takes to build an NxN addition table.
  double TimeAdd(ring R, long N)
  {
    vector<RingElem> value;
    value.reserve(N); // not necessary but good style: reserves memory (C++ STL)
    for (long i=1; i <= N; ++i)
      value.push_back(RingElem(R, i));

    RingElem tmp(R);
    const double StartTime = CpuTime();
    for (long i=0; i < N; ++i)
      for (long j=0; j < N; ++j)
        tmp = value[i] + value[j];
    return CpuTime()-StartTime;
  }

  // Measure how long it takes to build an NxN multiplication table.
  double TimeMult(ring R, long N)
  {
    vector<RingElem> value; value.reserve(N);
    for (long i=1; i <= N; ++i)
      value.push_back(RingElem(R, i));

    RingElem tmp(R);
    const double StartTime = CpuTime();
    for (long i=0; i < N; ++i)
      for (long j=0; j < N; ++j)
        tmp = value[i] * value[j];
    return CpuTime()-StartTime;
  }


  void test(long p)
  {
    cout << "Timing test for characteristic " << p << endl;

    constexpr long N=1000;

    cout << "Time to compute addition table of size " << N << endl;

    cout << "   ZZmod:        " << TimeAdd(NewZZmod(p), N) << endl;

    cout << "   RingFp:       ";
    try { cout << TimeAdd(NewRingFp(p), N) << endl; }
    catch (...) { cout << "CHAR TOO BIG\n"; }

    cout << "   RingFpLog:    ";
    try { cout << TimeAdd(NewRingFpLog(p), N) << endl; }
    catch (...) { cout << "CHAR TOO BIG\n"; }

    cout << "   RingFpDouble: ";
    try { cout << TimeAdd(NewRingFpDouble(p), N) << endl; }
    catch (...) { cout << "CHAR TOO BIG\n"; }

    cout << endl;

    cout << "Time to compute multiplication table of size " << N << endl;

    cout << "   ZZmod:        " << TimeMult(NewZZmod(p), N) << endl;

    cout << "   RingFp:       ";
    try { cout << TimeMult(NewRingFp(p), N) << endl; }
    catch (...) { cout << "CHAR TOO BIG\n"; }

    cout << "   RingFpLog:    ";
    try { cout << TimeMult(NewRingFpLog(p), N) << endl; }
    catch (...) { cout << "CHAR TOO BIG\n"; }

    cout << "   RingFpDouble: ";
    try { cout << TimeMult(NewRingFpDouble(p), N) << endl; }
    catch (...) { cout << "CHAR TOO BIG\n"; }

    cout << "--------------------------------------------------" << endl;

  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    for (long P = 256; P < power(2,30); P *= 16)
    {
      test(PrevPrime(P));
      if (P > numeric_limits<long>::max()/16) break;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-RingFp2.C,v 1.15 2024/03/08 19:59:58 abbott Exp $
// $Log: ex-RingFp2.C,v $
// Revision 1.15  2024/03/08 19:59:58  abbott
// Summary: Used constexpr (redmine 1511)
//
// Revision 1.14  2022/02/13 09:56:59  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.13  2015/07/01 16:31:35  abbott
// Removed superfluous "using namespace CoCoA"
//
// Revision 1.12  2015/06/29 15:47:58  bigatti
// -- code in namespace CoCoA
//
// Revision 1.11  2014/04/03 15:42:35  abbott
// Summary: Reduced parameter N (o/w too slow with debugging on some machines)
// Author: JAA
//
// Revision 1.10  2012/05/30 12:53:01  bigatti
// -- added comment for reserve
//
// Revision 1.9  2012/05/04 16:57:00  abbott
// Shortened some fn names.
// Corrected code (moved a print outside a try block), & improved readability.
//
// Revision 1.8  2012/05/02 17:22:39  abbott
// Reduced iteration count (for slower machines).
// Added check for overflow (on 32-bitters).
//
// Revision 1.7  2012/04/20 08:55:56  abbott
// Major change because the previous version was not usefully illustrative.
//
// Revision 1.6  2012/02/08 17:43:14  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.5  2011/08/23 12:04:04  bigatti
// -- updated after renaming ZZ --> BigInt
//
// Revision 1.4  2010/12/26 13:02:03  abbott
// Finished changing "GlobalXXXput()" into the corresponding standard C++ stream.
//
// Revision 1.3  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.2  2008/10/07 15:45:22  abbott
// Changed ErrorInfo objects so they include the name of their own error ID.
// Changed catch statements to catch const objects.
// Removed calls to the member fn which accessed the error ID member of an
// ErrorInfo; now you simply compare directly with the error ID (makes the
// code easier to read).
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.5  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.4  2007/02/26 17:41:50  bigatti
// -- getting ready for unique ring Z: using NewZmod(N), NewRingQ()
//
// Revision 1.3  2007/02/12 16:11:12  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.2  2007/02/10 18:44:04  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/10/17 10:46:53  cocoa
// Imported files
//
// Revision 1.1  2005/10/14 15:25:07  cocoa
// Major tidying and cleaning to small prime finite fields.
// Several consequential changes.  Improved their documentation.
//
// Added Makefile and script to include/CoCoA/ directory to
// keep library.H up to date.
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.3  2005/04/27 16:14:56  cocoa
// Cleaned up example programs -- added "free use" permit.
// Changed a couple of ErrorInfo object names, and added
// ERR::NotTrueGCDDomain.
//
// Revision 1.2  2005/04/21 15:12:19  cocoa
// Revised NewPolyRing as Dag Arneson suggested (perhaps just an interim
// measure).
// Brought example programs up to date (new name for CoCoA error
// information objects).
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.2  2004/12/09 15:08:42  cocoa
// -- added log info
//
