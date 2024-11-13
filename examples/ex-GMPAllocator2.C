// Copyright (c) 2005,2010,2012  John Abbott,  Anna M. Bigatti
// This file is part of the CoCoALib suite of examples.
// You are free to use any part of this example in your own programs.

#include "CoCoA/library.H"

using namespace std;

//----------------------------------------------------------------------
const string ShortDescription =
  "Program comparing memory allocators for GMP values; compares speed.  \n";

const string LongDescription =
  "This example shows the various ways of specifying the memory manager to  \n"
  "be used by the GMP library.  The choice of memory manager is indicated as\n"
  "an argument to the GlobalManager constructor.  Here we illustrate four   \n"
  "choices: the CoCoALib default, the system allocator, the specialized     \n"
  "GMPAllocator (with and without explicit indication of the slice size).   \n"
  "The program measures the speed of computation with these various choices.\n";
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


  void computation()
  {
    // Assumes GlobalManager has already been created.

    constexpr int Nmax = 123456;
    long MaxIters = 0;
    for (BigInt N(3); N <= Nmax; N += 2)
    {
      const long iters = NumIters(N);
      if (iters > MaxIters)
      {
        MaxIters = iters;
      }
    }
//  cout << "The longest sequence had length " << MaxIters << endl;
  }

  void program()
  {
    cout << ShortDescription << endl;

    {
      // Computation using CoCoALib's default choice of mem mgr for GMP values...
      GlobalManager CoCoAFoundations;
      cout << "Using CoCoALib's default mem mgr for GMP ..." << endl;
      const double t0 = CpuTime();
      computation();
      cout << "... time was " << CpuTime()-t0 << endl;
    }

    {
      // Computation using system allocator as mem mgr for GMP values...
      GlobalManager CoCoAFoundations(UseSystemAllocatorForGMP);
      cout << "Using standard system allocator mem mgr for GMP ..." << endl;
      const double t0 = CpuTime();
      computation();
      cout << "... time was " << CpuTime()-t0 << endl;
    }

    {
      // Computation using GMPAllocator as mem mgr for GMP values...
      GlobalManager CoCoAFoundations(UseGMPAllocator);
      cout << "Using GMPAllocator as mem mgr for GMP ..." << endl;
      const double t0 = CpuTime();
      computation();
      cout << "... time was " << CpuTime()-t0 << endl;
    }

    {
      // Computation using GMPAllocator(slice_size) as mem mgr for GMP values...
      GlobalManager CoCoAFoundations(UseGMPAllocator(256));
      cout << "Using GMPAllocator (with slice = 256 bytes) as mem mgr for GMP ..." << endl;
      const double t0 = CpuTime();
      computation();
      cout << "... time was " << CpuTime()-t0 << endl;
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
