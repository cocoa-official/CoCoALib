// Copyright (c) 2010  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program illustrates use of the pseudo-random number generator of CoCoALib.\n"
  "The numbers are independent and uniformly distributed in the given range; both \n"
  "ends of the range are reachable.                                               \n"
  "See RandomSeqBool if you want random bools,                                    \n"
  "  & RandomSeqLong if you want random machine integers.                         \n"
  "See also RandomSource for a general random generator.                          \n";

const string LongDescription =
  "CoCoALib offers a way to make uniform pseudo-random number generators.      \n"
  "When creating the generator you must specify the (inclusive) upper and lower\n"
  "bounds for the random numbers which will be generated.  When creating a     \n"
  "generator you may specify a `seed'; this allows different pseudo-random     \n"
  "sequences to be produced, though the sequence is completely determined by   \n"
  "its initial seed value.                                                     \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    cout << "Here are 20 random samples from the uniform distribution on [-10^10,+10^10]\n";
    const BigInt N = power(10,10);
    RandomSeqBigInt RndBigInt(-N,N);
    for (int i=0; i < 20; ++i)
      cout << NextValue(RndBigInt) << " ";
    cout << endl << endl;

    // If you prefer you can use a random sequence as an (endless) input iterator
    cout << "Here are 20 more random samples; this time from uniform distr on [0,10^99]\n";
    RandomSeqBigInt RndBigInt2(0,power(10,99));
    for (int i=0; i < 20; ++i)
    {
      ++RndBigInt2;
      cout << *RndBigInt2 << "\n";
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-RandomBigInt1.C,v 1.4 2022/02/13 09:56:58 abbott Exp $
// $Log: ex-RandomBigInt1.C,v $
// Revision 1.4  2022/02/13 09:56:58  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.3  2015/06/29 15:42:11  bigatti
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
// Revision 1.5  2012/03/30 10:40:27  abbott
// Changed variable names -- previously they were abbreviations of the original
// type name (which is about to change...again!)
//
// Revision 1.4  2011/08/23 12:04:04  bigatti
// -- updated after renaming ZZ --> BigInt
//
// Revision 1.3  2011/05/03 09:43:53  abbott
// Renamed RandomBitStream into RandomBoolStream (in cross-references).
//
// Revision 1.2  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.1  2010/06/29 15:17:02  abbott
// New class RandomZZStream
//
//
