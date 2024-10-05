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

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/ex-GMPAllocator1.C,v 1.7 2022/02/13 09:56:56 abbott Exp $
// $Log: ex-GMPAllocator1.C,v $
// Revision 1.7  2022/02/13 09:56:56  abbott
// Summary: Updated copyright (John & Anna in almost all cases, redmine 855)
//
// Revision 1.6  2015/06/26 15:34:48  abbott
// Summary: Moved code into namespace CoCoA (see redmine 739)
// Author: JAA
//
// Revision 1.5  2012/10/15 12:32:55  abbott
// Replaced size_t by long; improved LongDescription.
//
// Revision 1.4  2011/08/23 12:04:04  bigatti
// -- updated after renaming ZZ --> BigInt
//
// Revision 1.3  2010/12/26 13:02:03  abbott
// Finished changing "GlobalXXXput()" into the corresponding standard C++ stream.
//
// Revision 1.2  2010/12/17 16:07:55  abbott
// Ensured that all i/o in examples is on standard C++ streams
// (rather than GlobalInput(), etc).
//
// Revision 1.1  2010/10/22 09:14:47  abbott
// Renamed ex-GMPAllocator to ex-GMPAllocator1.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.5  2007/03/07 15:45:42  cocoa
// GMPAllocator must be created before GlobalManager, otherwise bus error.
//
// Revision 1.4  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.3  2007/02/12 15:29:07  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.2  2007/02/10 18:44:04  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/10/17 10:46:53  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.1  2005/04/29 15:42:02  cocoa
// Improved documentation for GMPAllocator.
// Added example program for GMPAllocator.
// Added example program for simple ops on polynomials.
// Added two new ctors for (principal) ideals (from long, and from BigInt).
// Added (crude) printing for PPMonoids.
// Updated library.H (#included GMPAllocator.H).
//
// Revision 1.3  2005/04/27 16:14:56  cocoa
// Cleaned up example programs -- added "free use" permit.
// Changed a couple of ErrorInfo object names, and added
// ERR::NotTrueGCDDomain.
//
// Revision 1.2  2005/04/21 15:12:19  cocoa
// Revised NewPolyRing as Dag Arneson suggested (perhaps just an interim
// measure).
// Brought example programs up to date (new name for CoCoA error
// information objects).
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.3  2004/12/09 15:08:42  cocoa
// -- added log info
//
