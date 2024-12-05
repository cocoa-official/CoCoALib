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
