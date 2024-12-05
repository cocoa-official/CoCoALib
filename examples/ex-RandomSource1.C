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
