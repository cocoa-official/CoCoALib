// Copyright (c) 2005  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

#include <cstdlib>
// using exit

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Inefficient program to compute sqrt(2) modulo a given number.    \n"
  "Simple example using finite fields or integers modulo N.\n";

const string LongDescription =
  "The program asks the user for the value of N, it creates the     \n"
  "ring of integers mod N, and finally uses \"brute force\" to find \n"
  "all square-roots of 2 modulo N.                                  \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    cout << "Enter a small number: ";
    int p;
    cin >> p;
    if ( !cin )
    {
      cerr << "*** ERROR *** Input must be a positive integer" << endl;
      exit(1);
    }
    if (p == 0)
    {
      cout << "There is no integer square-root of 2 modulo 0." << endl;
      return;
    }

    ring Fp = NewZZmod(p);
    if (!IsField(Fp))
      cout << "The number you entered, p=" << p << ", is not prime:" << endl
           << "so the integers mod p are implemented as a general QuotientRing."
           << endl;
    cout << "Fp is " << Fp << endl;

    // Just blindly try all elements of Fp (except -1)
    for (RingElem x(Fp); !IsZero(x+1); x += 1)
      if (x*x-2 == 0)
        cout << "A square-root of 2 modulo " << p
             << " is " << x << endl;
    cout << "Search finished" << endl;
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
