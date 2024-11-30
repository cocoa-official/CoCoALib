// Copyright (c) 2023  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Example of using a SumBigRat object to sum many rationals;        \n"
  "and an example of using a SumBigInt object to sum many integers.  \n";

const string LongDescription =
  "Example of using a SumBigRat object to sum many rationals.     \n"
  "SumBigRat can be rather faster than directly summing.          \n"
  "Similarly SumBigInt is suitable for summing many big integers. \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;
    SignalWatcher MonitorInterrupt(SIGINT); // you must also call CheckForInterrupt every so often

    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    // Sum some BigInts: normal summation would produce many carries.
    SumBigInt ZZsum;
    ZZsum += power(2,10000);
    for (int i=0; i < 1000000; ++i)
    {
      ZZsum += -i;
      ZZsum += i;
    }
    if (ZZsum.myTotal() != power(2,10000))
      CoCoA_THROW_ERROR("ZZsum is wrong", "ex-SumBigRat");

    // Sum the reciprocals of the first few primes
    SumBigRat SUM;
    for (PrimeSeq Pseq; *Pseq < 1000000; ++Pseq)
      SUM += BigRat(1,*Pseq);
    cout << FloatStr(SUM.myTotal()) << endl;
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

  catch (const CoCoA::InterruptReceived& intr)
  {
    cerr << endl
         << "------------------------------" << endl
         << ">>>  CoCoALib interrupted  <<<" << endl
         << "------------------------------" << endl
         << "-->>  " << intr << "  <<--" << endl;
    return 2;
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
