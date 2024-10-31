// Copyright (c) 2007  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program illustrates use of the pseudo-random bit generator of CoCoALib.  \n"
  "The bits are independent, identically distributed; each with equal probability\n"
  "of being true or false.  It is also possible to generate biased bits.         \n"
  "See RandomSeqLong & RandomSeqBigInt if you want random integers.              \n"
  "See also RandomSource for a general random generator.                         \n";

const string LongDescription =
  "CoCoALib offers a pseudorandom bit generator.  The generator can be   \n"
  "seeded when it is created; this allows different pseudo-random        \n"
  "sequences to be produced, though the sequence is completely determined\n"
  "by the initial seed value.  The `NextBiasedBool' function filters a   \n"
  "random bit sequence to produce `true' with a specified probability.   \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    RandomSeqBool RndBool;
    cout << "This is a brand new random bit generator: " << RndBool << endl << endl;

    constexpr int NumBits = 20;
    cout << "The first " << NumBits << " random bits are:";
    for (int i=0; i < NumBits; ++i)
      cout << " " << NextValue(RndBool);
    cout << endl;

    cout << endl
         << "The ""prob"" function simulates a biased coin toss." << endl
         << "Here we do 1000000 tosses of a 0.001 probability coin." << endl;
    constexpr double LowProb = 0.001;
    const size_t StartIndex = RndBool.myIndex();
    constexpr int NumIters = 1000000;
    int count = 0;
    for (int i=0; i < NumIters; ++i)
      if (NextBiasedBool(RndBool, LowProb)) ++count;
    cout << "The coin came up heads " << count << " times -- this count should be about " << LowProb*NumIters << endl
         << endl
         << "The ""prob"" function uses on average about two random bits per call." << endl
         << "The actual number of random bits used for this trial is " << RndBool.myIndex() - StartIndex << endl
         << "We see that this value is indeed not far from " << 2*NumIters << endl
         << endl;

    cout << "The final state of the generator is " << RndBool << endl;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-RandomBool1.C,v 1.6 2024/03/08 19:59:58 abbott Exp $
// $Log: ex-RandomBool1.C,v $
// Revision 1.6  2024/03/08 19:59:58  abbott
// Summary: Used constexpr (redmine 1511)
//
// Revision 1.5  2022/02/13 09:56:58  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.4  2015/06/29 15:42:11  bigatti
// *** empty log message ***
//
// Revision 1.3  2015/06/29 13:06:53  bigatti
// -- code in namespace CoCoA
//
// Revision 1.2  2013/02/15 16:30:40  abbott
// Consequential change: prob --> NextBiasedBool.
//
// Revision 1.1  2012/12/05 11:01:02  abbott
// Renamed existing example programs.
//
// Revision 1.4  2012/12/04 19:56:21  abbott
// Replaced calls to "sample" by calls to "NextValue".
//
// Revision 1.3  2012/12/04 09:59:28  abbott
// Improved doc and examples for pseudo-random generators:
// e.g. better var names in the examples.
//
// Revision 1.2  2011/08/23 12:04:04  bigatti
// -- updated after renaming ZZ --> BigInt
//
// Revision 1.1  2011/05/03 09:37:08  abbott
// Renamed RandomBitStream into RandomBoolStream.
//
// Revision 1.6  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.5  2010/06/29 15:14:13  abbott
// Improved descriptions: added reference to new RandomBigIntStream class.
//
// Revision 1.4  2010/02/16 10:20:58  abbott
// Added new fn sample; SampleBool & SampleLong are obsolescent.
//
// Revision 1.3  2010/02/01 22:35:40  abbott
// Very minor changes to some printed strings.
//
// Revision 1.2  2009/01/23 16:08:08  abbott
// Added missing newline inside the LongDescription string.
//
// Revision 1.1  2007/06/06 15:16:48  abbott
// Added new RandomBoolStream class (now based upon GMP's random generator).
// Consequential changes to Makefiles, etc.  There's even doc and an example!
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
