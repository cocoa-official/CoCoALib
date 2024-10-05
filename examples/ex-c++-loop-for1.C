// Copyright (c) 2017  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is an example showing some basic \"for\" loops in C++. \n"
  "It also shows \"continue\" and \"break\" inside a loop.     \n";

const string LongDescription =
  "This is an example showing some basic \"for\" loops in C++.    \n"
  "It also shows \"continue\" and \"break\" inside a loop.        \n"
  "We restrict to simple integer \"for\" loops.  See also examples\n"
  "for vectors and iterators.                                     \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  // SYNTAX of classical "for" loop:
  //
  // for  (<init>;  <keep-going-condition>;  <incr>)
  // {
  //   <commands>
  // }
  // In reality the "for" loop is more flexible than this.

  void program()
  {
    GlobalManager CoCoAFoundations;
    cout << ShortDescription << endl;

    // Print out the numbers from 1 to 10
    // Read "++i" as "increment i".
    for (int i=1; i <= 10; ++i)
    {
      cout << i << "  ";
    }
    cout << endl;

    // For DECREMENTING use "--i".  Print numbers 10 down to 1.
    for (int i=10; i >= 1; --i)
    {
      cout << i << "  ";
    }
    cout << endl;

    // "CONTINUE" command: skips to next iteration
    // Loop below prints: 1  2  3  8  9  10
    for (int i=1; i <= 10; ++i)
    {
      if (i > 3 && i < 8) continue;
      cout << i << "  ";
    }
    cout << endl;
    
    // "BREAK" command: exit entire loop immediately
    // Loop below prints: 1  2  3
    for (int i=1; i <= 10; ++i)
    {
      if (i > 3) break;
      cout << i << "  ";
    }
    cout << endl;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-c++-loop-for1.C,v 1.2 2022/02/13 09:57:00 abbott Exp $
// $Log: ex-c++-loop-for1.C,v $
// Revision 1.2  2022/02/13 09:57:00  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.1  2021/09/03 13:01:24  abbott
// Summary: Renamed from ex-c++-for-loop
//
// Revision 1.3  2019/03/27 14:57:37  bigatti
// -- better layout for example of C++ syntax
//
// Revision 1.2  2017/02/15 12:22:24  abbott
// Summary: Added countdown loop; several minor changes
//
// Revision 1.1  2017/02/10 16:31:25  abbott
// Summary: Added new examples for C++: ex-c++-basic, ex-c++-arith, ex-c++-for-loop
//
//
