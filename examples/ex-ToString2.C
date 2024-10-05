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

    BigInt N = factorial(27);
    string DecimalDigits = ToString(N);
    cout << "N = factorial(27) = " << N << endl;
    cout << "N as a string: " << DecimalDigits << endl;
    cout << "FloatStr(N): " << FloatStr(N) << endl;
    cout << endl;

    BigRat q(fibonacci(101),fibonacci(100));
    DecimalDigits = ToString(q);
    cout << "q = fib(101)/fib(100) = " << q << endl;
    cout << "q as a string: " << DecimalDigits << endl;
    cout << "FloatStr(q,10): " << FloatStr(q,10) << endl;
    cout << endl;

    cout << "FloatStr aids comprehension of large numbers: e.g." << endl;
    cout << "factorial(1000000) = " << FloatStr(factorial(1000000)) << endl;
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
