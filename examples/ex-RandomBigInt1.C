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
