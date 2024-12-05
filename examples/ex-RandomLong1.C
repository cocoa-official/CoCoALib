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
