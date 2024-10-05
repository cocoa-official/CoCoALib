// Copyright (c) 2010  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program illustrates the pseudo-random number generator of CoCoALib: \n"
  "  RandomSeqLong.                                                         \n"
  "The numbers are independent and uniformly distributed in the given range;\n"
  "both ends of the range are reachable.                                    \n"
  "See RandomSeqBool if you want random bools,                              \n"
  "and RandomSeqBigInt if you want random large integers.                   \n"
  "See also RandomSource for a general random generator.                    \n";

const string LongDescription =
  "CoCoALib offers a way to make uniform pseudo-random number generators.   \n"
  "When creating the generator you must specify the (inclusive)             \n"
  "upper and lower bounds for the random numbers which will be generated.   \n"
  "When creating a generator you may specify a `seed';                      \n"
  "this allows different pseudo-random sequences to be produced,            \n"
  "though the sequence is completely determined by its initial seed value.  \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    cout << "Here are 20 random integers in the range [-10,+10]\n";
    RandomSeqLong RndLong(-10,10);
    for (int i=0; i < 20; ++i)
      cout << NextValue(RndLong) << " ";
    cout << endl << endl;

    // If you prefer you can use RndLong as an (endless) input iterator...
    cout << "Here are 20 more random integers in the range [0,99]\n";
    RandomSeqLong RndLong2(0,99);
    for (int i=0; i < 20; ++i)
    {
      ++RndLong2;
      cout << *RndLong2 << " ";
    }
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-RandomLong1.C,v 1.5 2022/02/13 09:56:59 abbott Exp $
// $Log: ex-RandomLong1.C,v $
// Revision 1.5  2022/02/13 09:56:59  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.4  2017/03/01 17:16:43  bigatti
// -- resorted ShortDescription and LongDescription
//
// Revision 1.3  2015/06/29 15:43:36  bigatti
// *** empty log message ***
//
// Revision 1.2  2015/06/29 13:06:53  bigatti
// -- code in namespace CoCoA
//
// Revision 1.1  2012/12/05 11:01:02  abbott
// Renamed existing example programs.
//
// Revision 1.7  2012/12/04 19:56:21  abbott
// Replaced calls to "sample" by calls to "NextValue".
//
// Revision 1.6  2012/12/04 09:59:28  abbott
// Improved doc and examples for pseudo-random generators:
// e.g. better var names in the examples.
//
// Revision 1.5  2012/11/30 15:19:34  abbott
// Added reference to RandomSource in short description.
// Improved a comment.
//
// Revision 1.4  2011/05/03 09:43:53  abbott
// Renamed RandomBitStream into RandomBoolStream (in cross-references).
//
// Revision 1.3  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.2  2010/06/29 15:14:13  abbott
// Improved descriptions: added reference to new RandomZZStream class.
//
// Revision 1.1  2010/02/16 10:19:29  abbott
// Added new class RandomLongStream with example and test.
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
