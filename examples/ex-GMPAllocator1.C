// Copyright (c) 2005,2010,2012  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Program to find lengths of 3n+1 sequences -- using CoCoA::GMPAllocator  \n";

const string LongDescription =
  "This example shows how to specify that a CoCoALib custom allocator be used\n"
  "for managing the memory for small GMP values.  If you do computations    \n"
  "involving many small big-integers (e.g. up to about 40 decimal digits) then\n"
  "the custom allocator may give slightly better run-time performance (with \n"
  "debugging turned off!)  The argument to the constructor for GlobalManager\n"
  "indicates that the CoCoALib allocator is to be used for GMP values.      \n";
//----------------------------------------------------------------------

namespace CoCoA
{

// This function counts the number of iterations of the 3n+1 sequence
// are needed to reach 1 (the first time) starting from N.  We chose
// this example because it performs many operations on smallish values.
  long NumIters(BigInt N)
  {
    N = abs(N); // ignore sign of N
    long iters = 0;
    while (N > 1)
    {
      if (IsEven(N)) N /= 2;  // the defining relation for the sequence
      else N = 3*N+1;         //
      ++iters;
    }
    return iters;
  }


  void program()
  {
    GlobalManager CoCoAFoundations(UseGMPAllocator);

    cout << ShortDescription << endl;

    BigInt Nmax;
    cout << "Enter highest starting value to try (e.g. 1000000): ";
    cin >> Nmax;
    long MaxIters = 0;

    for (BigInt N(3); N <= Nmax; N += 2)
    {
      long iters = NumIters(N);
      if (iters > MaxIters)
      {
        MaxIters = iters;
        cout << "The sequence starting from " << N << " has length " << MaxIters << endl;
      }
    }
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
