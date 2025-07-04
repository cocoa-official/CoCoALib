// Copyright (c) 2017  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is an example showing some basic arithmetic ops in C++.   \n";

const string LongDescription =
  "This is an example showing some basic arithmetic ops in C++. \n"
  "With a strong caution about division and computing powers.   \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;
    cout << ShortDescription << endl;

    const long a = 3;
    const long b = 4;

    // Most arithmetic operations work "as expected"
    cout << a << "+" << b << " = " << a+b << endl; // sum
    cout << a << "-" << b << " = " << a-b << endl; // difference
    cout << a << "*" << b << " = " << a*b << endl; // product
    cout << a << "/" << b << " = " << a/b << "  <--- INTEGER DIVISION" << endl;

    // ***** TWO IMPORTANT WARNINGS BELOW *****

    // ==========================================================
    // ***NO RATIONALS***  ***NO RATIONALS***  ***NO RATIONALS***
    cout << "This is NOT one half: 1/2 gives " << 1/2 << endl << endl;
    // ***NO RATIONALS***  ***NO RATIONALS***  ***NO RATIONALS***
    // ==========================================================
    
    // =================================================
    // ***NOT POWER***  ***NOT POWER***  ***NOT POWER***
    cout << a << "^" << b << " = " << (a^b) << "  <--- ***NOT POWER***" << endl;
    // ***NOT POWER***  ***NOT POWER***  ***NOT POWER***
    // =================================================
  }

} // end of namespace CoCoA

// IGNORE THE STUFF BELOW (at least for now)

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
