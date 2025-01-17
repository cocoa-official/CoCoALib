// Copyright (c) 2021  John Abbott,  Anna M. Bigatti
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

  // SYNTAX of "while" loop: -- may do 0 iterations
  //
  // while  (<keep-going-condition>)
  // {
  //   <commands>
  // }

  // SYNTAX of "do..while" loop: -- always does at least 1 iter
  //
  // do
  // {
  //   <commands>
  // }
  // while  (<keep-going-condition>);


  void program()
  {
    GlobalManager CoCoAFoundations;
    cout << ShortDescription << endl;

    //--------------------------------------------
    // Print out the primes up to 20
    long p = NextPrime(1);
    while (p <= 20)
    {
      cout << p << "  ";
      p = NextPrime(p);
    }
    cout << endl;

    //--------------------------------------------
    // Print out the primes whose least prim root is >= 5
    // until one is found with least prim root >= 10
    p = NextPrime(1);
    while (true) // exit condition must appear inside loop body
    {
      long r = PrimitiveRoot(p);
      if (r >= 5) { cout << p << "  "; }
      if (r >= 10) break; // exit condition calls "break"
      p = NextPrime(p);
    }
    cout << endl;

    //--------------------------------------------
    // Alternative for of search in loop above:
    long r = 0;
    do
    {
      p = NextPrime(p);
      r = PrimitiveRoot(p);
    }
    while (r < 10);
    cout << "Smallest prime with first PrimRoot >= 10 is " << p << endl;
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
