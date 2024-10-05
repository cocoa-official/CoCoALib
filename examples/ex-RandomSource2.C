// Copyright (c) 2010  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program shows basic use of a RandomSource object to produce\n"
  "random booleans, machine integers, and big integers.            \n"
  "Also shows how to seed and reseed a RandomSource.               \n"
;

const string LongDescription =
  "This program uses explicitly a RandomSource object to generate \n"
  "distributed random bits/booleans, machine integers in a given  \n"
  "range, and big integers in a given range.                      \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    RandomSource src;

    cout << "*** Random booleans ***\n";
    constexpr long NumTrials = 15;
    cout << NumTrials << " times RandomBool(src); -->  ";
    for (long i=0; i < NumTrials; ++i)
      cout << RandomBool(src) << " ";
    cout << endl << endl;

    cout << "*** Random Machine Integers ***\n";
    RandomSource src1(123);
    cout << "(seed 123)   6 times RandomLong(src1, 10, 90); -->  "
         << RandomLong(src1, 10, 90) << " "
         << RandomLong(src1, 10, 90) << " "
         << RandomLong(src1, 10, 90) << " "
         << RandomLong(src1, 10, 90) << " "
         << RandomLong(src1, 10, 90) << " "
         << RandomLong(src1, 10, 90);
    cout << endl;
    
    reseed(src1, 123);
    cout << "(reseed 123) 6 times RandomLong(src1, 10, 90); -->  "
         << RandomLong(src1, 10, 90) << " "
         << RandomLong(src1, 10, 90) << " "
         << RandomLong(src1, 10, 90) << " "
         << RandomLong(src1, 10, 90) << " "
         << RandomLong(src1, 10, 90) << " "
         << RandomLong(src1, 10, 90);
    cout << endl << endl;

    cout << "*** Random Large Integers ***\n";
    const BigInt lo = BigIntFromString("111111111111111111");
    const BigInt hi = BigIntFromString("999999999999999999");
    cout << "3 times RandomLong(src, lo, hi); -->  " << "\n"
         << RandomBigInt(src, lo, hi) << "\n"
         << RandomBigInt(src, lo, hi) << "\n"
         << RandomBigInt(src, lo, hi) << "\n";
    cout << endl;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-RandomSource2.C,v 1.7 2024/03/08 19:59:58 abbott Exp $
// $Log: ex-RandomSource2.C,v $
// Revision 1.7  2024/03/08 19:59:58  abbott
// Summary: Used constexpr (redmine 1511)
//
// Revision 1.6  2022/02/13 09:56:59  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.5  2018/04/20 18:51:24  abbott
// Summary: Changed ctors for BigInt/BigRat from string or from MPZ/MPQ
//
// Revision 1.4  2015/06/29 15:43:36  bigatti
// *** empty log message ***
//
// Revision 1.3  2015/06/29 13:06:53  bigatti
// -- code in namespace CoCoA
//
// Revision 1.2  2012/12/04 19:55:45  abbott
// Corrected typo in a description string.
//
// Revision 1.1  2012/02/03 10:34:53  bigatti
// first import
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
