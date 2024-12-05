// Copyright (c) 2004-2007,2009  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Program to calculate the modular inverse of a number.  \n"
  "You must give as input the modulus and the residue;    \n"
  "the residue need not be reduced.                       \n";

const string LongDescription =
  "This program illustrates that BigInt values can be used much like normal C++ \n"
  "ints except that there is a almost no limit on the magnitude of the values.  \n"
  "NB If you need extreme efficiency then use the GMP library directly.         \n"
  "Contrast this example with ex-RingZZ1.                                       \n";
//-----------------------------------------------------------------------------


namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    cout << "Input the modulus: ";
    BigInt modulus;
    cin >> modulus;
    if (!cin || modulus < 2) { cerr << "The modulus must be an integer >= 2\n"; return; }

    cout << "Input the residue: ";
    BigInt residue;
    cin >> residue;
    if (!cin) { cerr << "The residue must be an integer.\n"; return; }

    if (gcd(residue, modulus) != 1)
    {
      cout << "No inverse exists because the residue and modulus have a" << endl
           << "non-trivial common factor: " << gcd(residue, modulus) << endl;
      return;
    }

    const BigInt inverse = InvMod(residue, modulus); // defined in NumTheory
    cout << "The inverse of " << residue << " modulo " << modulus
         << " is " << inverse << endl;

    // We confirm that the result is right.
    if ((inverse*residue)%modulus != 1)
      cout << "THIS SHOULD NEVER BE PRINTED!" << endl;
  }

} // end of namespace CoCoA


// We write main() like this so we can handle uncaught CoCoA errors in
// a sensible way (i.e. by announcing them).
int main()
{
  try
  {
    CoCoA::program();
    return 0;
  }
  catch (const CoCoA::ErrorInfo& err)
  {
    cerr << "***ERROR***  UNCAUGHT CoCoA error" << endl;
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
