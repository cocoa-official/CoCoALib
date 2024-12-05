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
