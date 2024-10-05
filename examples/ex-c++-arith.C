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

    long a = 3;
    long b = 4;

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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-c++-arith.C,v 1.3 2022/02/13 09:57:00 abbott Exp $
// $Log: ex-c++-arith.C,v $
// Revision 1.3  2022/02/13 09:57:00  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.2  2021/09/03 12:58:16  abbott
// Summary: Minor improvements
//
// Revision 1.1  2017/02/10 16:31:25  abbott
// Summary: Added new examples for C++: ex-c++-basic, ex-c++-arith, ex-c++-for-loop
//
//
