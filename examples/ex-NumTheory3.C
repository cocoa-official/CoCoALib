// Copyright (c) 2014  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This example uses eratosthenes to build a sieve for testing quickly   \n"
  "whether a number is prime.  We compute many Goldbach representations. \n";

const string LongDescription =
  "This program tests how hard it is to find a \"Goldbach\" representation \n"
  "of an integer; i.e. a sum of two primes.  It has to do many primality  \n"
  "tests, so it is faster to use a table than call IsPrime repeatedly.    \n";
//----------------------------------------------------------------------

namespace CoCoA
{

  // See how hard it is to find a Goldbach representation for even
  // numbers up to NMAX; print out successive "maxima".
  void SimpleSearch(long NMAX)
  {
    long BiggestSmallPrime = 3;
    const double t0 = CpuTime();
    for (long n = 6; n <= NMAX; n += 2)
    {
      for (long p=3;; p = NextPrime(p))
        if (IsPrime(n-p))
        {
          if (p > BiggestSmallPrime)
          {
            BiggestSmallPrime = p;
            cout << "For n=" << n << " simplest repr is " << p << "+" << n-p << endl;
          }
          break;
        }
    }
    const double t1 = CpuTime();
    cout << "Time for SimpleSearch (using IsPrime): " << t1-t0 << endl;
  }


  // Does the same as SimpleSearch, but uses a table to store the first
  // few primes, and a boolean table to test for primality.
  // Harder to understand than SimpleSearch, but much faster!
  void TblSearch(long NMAX)
  {
    const double t0 = CpuTime();
    vector<long> prime;
    for (int p=3; p < 1000; p = NextPrime(p))
      prime.push_back(p);
    const int nprimes = len(prime);
    const vector<bool> IsPrimeTbl = eratosthenes(NMAX);
    long BiggestSmallPrime = 3;
    for (long n = 6; n <= NMAX; n += 2)
    {
      for (int i=0; i < nprimes; ++i)
      {
        if (IsPrimeTbl[(n-prime[i])/2])
        {
          if (prime[i] > BiggestSmallPrime)
          {
            BiggestSmallPrime = prime[i];
            cout << "For n=" << n << " simplest repr is " << prime[i] << "+" << n-prime[i] << endl;
          }
          break;
        }
      }
      
    }
    const double t1 = CpuTime();
    cout << "Time for TblSearch (using eratosthenes): " << t1-t0 << endl;
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;
    cout << boolalpha; // so that bools print out as true/false

    const long NMAX = 555555;
    SimpleSearch(NMAX);
    cout << endl;
    TblSearch(NMAX);
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
