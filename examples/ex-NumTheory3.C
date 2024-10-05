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
  "This program tests how hard it is to find a \"Goldbach\" represetation \n"
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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-NumTheory3.C,v 1.5 2022/02/13 09:56:57 abbott Exp $
// $Log: ex-NumTheory3.C,v $
// Revision 1.5  2022/02/13 09:56:57  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.4  2019/09/16 14:40:15  abbott
// Summary: Improved readability
//
// Revision 1.3  2015/06/29 15:27:25  bigatti
// -- minor fix
//
// Revision 1.2  2015/06/29 12:20:43  bigatti
// -- moved in cocoa namespace
//
// Revision 1.1  2014/09/16 10:41:41  abbott
// Summary: Added new fn eratosthenes (with doc, example, test)
// Author: JAA
//
// Revision 1.7  2013/05/28 07:07:04  bigatti
// -- added "cout << boolalpha": useful for testing
//
// Revision 1.6  2012/11/30 14:04:55  abbott
// Increased visibility of comment saying "put your code here".
//
// Revision 1.5  2010/12/17 16:07:54  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
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
