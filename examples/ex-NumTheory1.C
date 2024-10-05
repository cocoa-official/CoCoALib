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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-NumTheory1.C,v 1.13 2022/02/13 09:56:57 abbott Exp $
// $Log: ex-NumTheory1.C,v $
// Revision 1.13  2022/02/13 09:56:57  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.12  2022/02/02 09:19:06  abbott
// Summary: Renamed SmoothFactor to factor_TrialDiv (redmine 950)
//
// Revision 1.11  2015/06/29 15:27:25  bigatti
// -- minor fix
//
// Revision 1.10  2015/06/29 12:20:43  bigatti
// -- moved in cocoa namespace
//
// Revision 1.9  2014/03/24 12:09:20  abbott
// Summary: Major revision to public interface of factorization template class
// Author: JAA
//
// Revision 1.8  2013/02/26 11:28:35  abbott
// Corrected minor mistake in a printed message.
//
// Revision 1.7  2012/10/05 09:29:32  abbott
// Changed myExponents into myMultiplicities.
//
// Revision 1.6  2012/05/30 16:21:03  abbott
// Cleaned and improved the leggibility of program output.
// Approved by Alessandra!
//
// Revision 1.5  2012/04/27 15:05:52  bigatti
// -- clarified
//
// Revision 1.4  2011/08/23 12:04:04  bigatti
// -- updated after renaming ZZ --> BigInt
//
// Revision 1.3  2011/05/24 09:58:54  abbott
// Slight change to one test to make program faster on some platforms.
//
// Revision 1.2  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.1  2010/03/03 10:41:50  abbott
// Added example for basic number theory functions.
//
// Revision 1.4  2008/10/07 12:12:54  abbott
// Removed useless commented out #include.
//
// Revision 1.3  2007/05/31 16:06:16  bigatti
// -- removed previous unwanted checked-in version
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.9  2007/03/07 11:51:40  bigatti
// -- improved test alignment
//
// Revision 1.8  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.7  2007/03/02 17:46:40  bigatti
// -- unique RingZ and RingQ
// -- requires foundations.H ;  foundations blah;  (thik of a better name)
//
// Revision 1.6  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.5  2007/03/01 13:52:59  bigatti
// -- minor: fixed typo
//
// Revision 1.4  2007/02/28 15:15:56  bigatti
// -- minor: removed quotes in description
//
// Revision 1.3  2007/02/12 16:27:43  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.2  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.1  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
