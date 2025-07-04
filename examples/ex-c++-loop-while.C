// Copyright (c) 2021  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This is an example showing some basic \"while\" loops in C++.     \n"
  "See ex-c++-loop-for1.C for the \"continue\" and \"break\" commands. \n";

const string LongDescription =
  "This is an example showing some basic \"while\" loops in C++.     \n"
  "See ex-c++-loop-for1.C for the \"continue\" and \"break\" commands. \n";

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
    cout << "All primes up to 20:" << endl;
    long p = NextPrime(1);
    while (p <= 20)
    {
      cout << p << "  ";
      p = NextPrime(p);
    }
    cout << endl << endl;
    // At this point the variable p contains the smallest prime greater than 20.
    

    //--------------------------------------------
    // Print out the primes whose least primitive root is >= 5
    // until one is found with least primitive root >= 10.
    cout << "First few primes with least primitive root >= 5:" << endl;
    p = NextPrime(1);
    while (true) // exit condition must appear inside loop body
    {
      const long r = PrimitiveRoot(p);
      if (r >= 5)  { cout << p << "  "; }
      if (r >= 10)  break; // exit loop using "break" (see also ex-c++-loop-for1.C)
      p = NextPrime(p);
    }
    cout << endl << endl;
    // Here the variable p contains the smallest prime with primitive root >= 10

    // NOTE: it is common to write  while (true) { ... if (exit_cond) break; ... }
    //       when the exit condition is known only part way through an iteration.

    //--------------------------------------------
    // Alternative form of search in loop above:
    p = 1;  // reset p to start from 1
    long r = 0; // now outside the loop body; initial value is unimportant.
    do
    {
      p = NextPrime(p);
      r = PrimitiveRoot(p);
    }
    while (r < 10);  // logical negation of the exit condition
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
