// Copyright (c) 2014  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program shows how to use some of the \"template utilities\" \n"
  "in CoCoALib: e.g. sum, product, LexCmp3.                         \n";

const string LongDescription =
  "This program illustrates use of product, sum, and LexCmp3.      \n"
  "The functions are fairly general, but we present just a simple  \n"
  "case.  They can also be used with lists instead of vectors      \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    const vector<long> OneToTen = LongRange(1,10);
    cout << "Factorial of 10 is " << product(OneToTen) << endl;
    cout << "1+2+3+...+10 = " << sum(OneToTen) << endl;

    vector<long> v1230(4);
    v1230[0] = 1; v1230[1] = 2; v1230[2] = 3; v1230[3] = 0;
    cout << "Doing Lex comparison of the following two vectors:\n"
         << OneToTen << endl
         << v1230 <<endl;
    const int result = LexCmp3(OneToTen.begin(),OneToTen.end(), v1230.begin(),v1230.end());
    cout << "Result is " << result << endl;
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
