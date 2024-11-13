// Copyright (c) 2004  John Abbott,  Anna M. Bigatti
// Orig authors: 2004 Daniele Venzano; modified by John Abbott.

#include "CoCoA/library.H"

using namespace std;

//-----------------------------------------------------------------------------
const string ShortDescription =
  "Program to find a (probable) prime with a specified number of bits.\n";


const string LongDescription =
  "This program shows that BigInts can be used much like normal C++ ints\n"
  "with the advantage that there is almost no limit on the magnitude of \n"
  "the values.  Here we generate random BigInts and test them for being \n"
  "probable primes -- stopping as soon as we find a likely prime.       \n"
  "NB If you need extreme efficiency then use the GMP library directly. \n";
//-----------------------------------------------------------------------------


namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    cout << "Enter the number of bits the prime should have: ";
    long nbits;
    cin >> nbits;
    if (!cin) CoCoA_THROW_ERROR(ERR::InputFail, "ex-BigInt3: input was not a (small) integer");
    if (nbits < 2)
    {
      cout << "There are no primes with fewer than 2 bits." << endl;
      return;
    }

    cout << endl << "Starting search for a prime with " << nbits << " bits..." << endl;
    if (nbits > 9999) cout << "WARNING: the computation will take a very long time!" << endl;
    else if (nbits > 999) cout << "NOTE: this computation may take some time." << endl;

    const BigInt min = power(2, nbits-1);
    const BigInt max = 2*min-1;

    BigInt candidate = RandomBigInt(min, max);
    while (!IsProbPrime(candidate))
    {
      candidate = RandomBigInt(min, max);
    }

    cout << candidate << " seems to be a prime! (" << nbits << " bits)" << endl;

    cout << "And further test with more iterations ";
    if (IsProbPrime(candidate, 50))
      cout << "confirms";
    else
      cout << "rejects";
    cout << " this result." << endl;
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
