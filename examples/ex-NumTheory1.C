// Copyright (c) 2006,2012  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "This program illustrates the use of some basic number theoretic functions.\n";

const string LongDescription =
  "This programs show how to use some of the basic number theoretic functions.\n"
  "Many of the examples use machine integers for convenience, but all the     \n"
  "functions also work with BigInt values (except NextPrime and PrevPrime).   \n";

//----------------------------------------------------------------------

namespace CoCoA
{

  void PrettyPrint(const factorization<long>& FacInfo)
  {
    const vector<long>& facs = FacInfo.myFactors();
    const vector<long>& mults = FacInfo.myMultiplicities();
    const long NumFacs = len(facs);
    for (long i=0; i < NumFacs; ++i)
      cout << facs[i] << "^" << mults[i] << " * ";
    cout << FacInfo.myRemainingFactor() << endl << endl;
  }

  void program()
  {
    GlobalManager CoCoAFoundations;

    cout << ShortDescription << endl;

    const int n = 123456789;
    const int m = 987654321;
    cout << "Some example computations with" << endl
         << "  m=" << m << "  and" << endl
         << "  n=" << n << endl
         << endl;

    cout << "GCD computation  (NB result is always non-negative)" << endl;
    cout << "  gcd(m,n) = " << gcd(m,n) << endl;
    cout << "  gcd(-m,-n) = " << gcd(-m,-n) << "  -- result is positive!" << endl;
    cout << endl;
    // Compute the cofactors for m and n: note that a & b must be of type long!
    long a,b;
    ExtGcd(a,b,m,n);
    cout << "Cofactors for gcd(m,n) are: a = " << a << "    b = " << b << endl
         << endl;
    cout << "To compute lcm in this case we must use big integers because" << endl
         << "the result is too big to fit into a 32-bit machine integer" << endl
         << "  lcm(m,n) = " << lcm(BigInt(m), n) << endl
         << endl;

    cout << "Partial factorization (ex: looking for factors <= 99)" << endl;
    cout << "  factor_TrialDiv(m,99) = " << endl
         << "    " << factor_TrialDiv(m,99) << endl;
    cout << "  which means  m = ";  PrettyPrint(factor_TrialDiv(m,99));

    factorization<long> nfactors = factor_TrialDiv(n,99);
    cout << "  factor_TrialDiv(n,99) = " << endl
         << "    " << nfactors << endl;
    cout << "  which means  n = ";  PrettyPrint(nfactors);

    cout << endl
         << "Complete factorization" << endl;
    nfactors = factor(n);
    cout << "  factor(n) = " << endl
         << "    " << nfactors << endl;
    cout << "  which means  n = "; PrettyPrint(nfactors);

    // Testing primality & generating primes
    cout << endl
         << "Testing primality & generating primes." << endl;

    // Example showing that IsPrime can be much slower than IsProbPrime.
    const BigInt N = power(2,64);
    const BigInt P = NextProbPrime(N+300); // NB NextPrime is *only* for machine integers!
    cout << "Comparing IsPrime and IsProbPrime on the (probable) prime P = 2^64 + " << P-N << endl;
    double StartTime = CpuTime();
    IsProbPrime(P);
    cout << "Time for IsProbPrime(P): " << CpuTime() - StartTime << "     [relatively quick]" << endl;
    StartTime = CpuTime();
    IsPrime(P);
    cout << "Time for IsPrime(P):     " << CpuTime() - StartTime << "     [relatively slow]" << endl;
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
