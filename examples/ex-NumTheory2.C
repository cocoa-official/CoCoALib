// Copyright (c) 2015  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program illustrates use of CRTMill to build a large integer \n"
  "from its residues modulo various different primes.               \n";

const string LongDescription =
  "This program illustrates use of CRTMill to build a large integer \n"
  "from its residues modulo various different primes.               \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    // Daft example: we shall construct N from its modular images.
    // We pretend we know only the UPB, and have a way to compute
    // N modulo p for any prime p (in this case we just compute N%p).
    const BigInt N = power(10,100);
    const BigInt UPB = 2*N+1;

    CRTMill crt;
    int p = 101;
    while (true)
    {
      p = NextPrime(p);
      crt.myAddInfo(N%p, p); // tell crt the new residue-modulus pair
      if (CombinedModulus(crt) >= UPB) break;
    }

    // Since we already know the answer, we can check it is correct.
    if (CombinedResidue(crt) != N)
      CoCoA_THROW_ERROR("Wrong answer", "CoCoA::Program");
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
