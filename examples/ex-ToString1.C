// Copyright (c) 2011  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program illustrates use of the decimal string conversion functions.\n";

const string LongDescription =
  "Example of use FloatStr, ScientificStr and DecimalStr for printing \n"
  "rationals in a comprehensible/decimal format.                      \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    // Compute an approximation to exp(5) using the power series.
    BigRat ans(1,1);
    BigRat term(1,1);
    BigRat x(5,1);

    cout << "Successive sums of the series for exp(5):" << endl;
    for (int i=1; i < 16; ++i)
    {
      term = (term*x)/i;
      ans += term;
      cout << FloatStr(ans) << "\t  "
           << DecimalStr(ans) << "\t  "
           << ScientificStr(ans) << endl;
    }
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
