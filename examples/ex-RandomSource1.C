// Copyright (c) 2010  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program shows the easiest way to produce                   \n"
  "random booleans, machine integers, and big integers.            \n";

const string LongDescription =
  "This program uses GlobalRandomSource() to generate uniformly   \n"
  "distributed random bits/booleans, machine integers in a given  \n"
  "range, and big integers in a given range.  For the booleans and\n"
  "machine integers it produces a histogram.                      \n"
  "Use of globals could be dangerous in multi-threaded programs.  \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void HistogramRandomBool(const long NumTrials)
  {
    cout << "*** Random booleans ***\n";
    cout << "*** RandomBool() ***\n";
    int TrueCount = 0;
    for (long i=0; i < NumTrials; ++i)
      if (RandomBool())
        ++TrueCount;
    cout << "Histogram after " << NumTrials << " trials:\n"
         << "  True  = " << TrueCount << "\n"
         << "  False = " << NumTrials - TrueCount << "\n"
         << endl;
  }


  void HistogramRandomLong(const long NumTrials, const long lo, const long hi)
  {
    cout << "*** Random Machine Integers ***\n";
    cout << "*** RandomLong(" << lo << ", " << hi << ") ***\n";
    const long NumCases = hi-lo+1;
    const long TotTrials = NumTrials*NumCases;
    vector<long> hist(NumCases);
    for (long i=0; i < TotTrials; ++i)
      ++hist[-lo + RandomLong(lo, hi)];
    cout << "Histogram after " << NumTrials << "*" << NumCases << " trials:\n";
    for (int k=0; k < NumCases; ++k)
      cout << "  " << k+lo << "\t " << hist[k] << "\n";
    cout << endl;
  }


  void ExRandomBigInt(const long NumTrials, const long NumBits)
  {
    cout << "*** Random Large Integers ***\n";
    cout << "*** RandomBigInt(-upb, upb) ***\n";
    const BigInt upb = power(2, NumBits)-1;
    cout << "Here are " << NumTrials << " signed random numbers of length "
         << NumBits << " bits.\n";
    for (long i=0; i < NumTrials; ++i)
      cout << RandomBigInt(-upb, upb) << "\n";
    cout << endl;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    HistogramRandomBool(2000000);
    HistogramRandomLong(1000, -2, 7);
    ExRandomBigInt(5, 100);
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-RandomSource1.C,v 1.6 2022/02/13 09:56:59 abbott Exp $
// $Log: ex-RandomSource1.C,v $
// Revision 1.6  2022/02/13 09:56:59  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.5  2015/06/29 15:43:36  bigatti
// *** empty log message ***
//
// Revision 1.4  2015/06/29 13:06:53  bigatti
// -- code in namespace CoCoA
//
// Revision 1.3  2012/02/03 10:10:02  bigatti
// -- reorganized and simplified: implicit use of GlobalRandomSource
//
// Revision 1.2  2011/08/23 12:04:04  bigatti
// -- updated after renaming ZZ --> BigInt
//
// Revision 1.1  2010/12/17 16:01:34  abbott
// New example for RandomSource.
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
