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
