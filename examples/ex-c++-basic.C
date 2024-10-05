// Copyright (c) 2017  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is an example showing some very basic C++:  comments, \n"
  "variable declaration (and init) for basic types and how to print. \n";

const string LongDescription =
  "This is an example showing some very basic C++:  comments, \n"
  "variable declaration (and init) for basic types and        \n"
  "how to print.  See also ex-c++-bool.C                      \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  // Here is a very basic intro to C++.
  // You should consult a good introductory book.
  // The "Effective" books by Scott Meyers contain many useful
  // design hints (but they assume you know basic C++ already).

  void program()
  {
    GlobalManager CoCoAFoundations;
    cout << ShortDescription << endl;

    // This is a comment... up to the end of the line --->

    /* This is a comment...
       ... up to the end-of-comment marker --> */

    // Some BASIC TYPES in C++
    // Note: = denotes initialization or assignment; == is equality test.
    int a = 1;  // an integer value (of limited range)
    long b = 2; // an integer value (of greater range)
    string str = "a succession of characters"; // string is part of STL
    double pi = 3.14159; // a floating-point value
    bool flag = true; // boolean value: either "true" or "false"
    
    //-- PRINTING -------------------------------------------------
    // The ostream `cout' is the computer screen (usually).
    // The operator `<<' means "send to" or "print on".
    // The manipulator `endl' ends the current line; following
    // output will start on the next line.

    cout << "Hello, World!" << endl;

    cout << "The string variable `str' has value " << str << endl;

    cout << "pi is approximately " << pi << endl;
    cout << endl;  // just produces an empty line
    
    cout << a << "+" << b << " = " << a+b << endl;

    cout << "The boolean variable `flag' has value " << flag << endl;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-c++-basic.C,v 1.8 2022/03/17 14:38:39 abbott Exp $
// $Log: ex-c++-basic.C,v $
// Revision 1.8  2022/03/17 14:38:39  abbott
// Summary: Improved comments
//
// Revision 1.7  2022/02/13 09:57:00  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.6  2021/10/07 14:27:45  abbott
// Summary: Added comment
//
// Revision 1.5  2021/09/03 12:58:41  abbott
// Summary: Improved comments
//
// Revision 1.4  2017/02/24 16:58:00  abbott
// Summary: Removed an empty line
//
// Revision 1.3  2017/02/15 12:21:53  abbott
// Summary: Added bool as a basic type
//
// Revision 1.2  2017/02/10 16:40:35  abbott
// Summary: Minor improvement
//
// Revision 1.1  2017/02/10 16:31:25  abbott
// Summary: Added new examples for C++: ex-c++-basic, ex-c++-arith, ex-c++-for-loop
//
//
