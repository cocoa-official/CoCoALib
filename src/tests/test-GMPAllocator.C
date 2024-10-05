//   Copyright (c)  2005  John Abbott,  Anna M. Bigatti

//   This file is part of the source of CoCoALib, the CoCoA Library.

//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.

//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.

//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.


#include "CoCoA/BigIntOps.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/error.H"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

namespace CoCoA
{

  // Test GMPAllocator by computing many "3n+1" sequences, and printing the successive maxima.

  // This function counts the number of iterations of the 3n+1 sequence
  // are needed to reach 1 (the first time) starting from N
  long NumIters(BigInt N)
  {
    N = abs(N); // ignore sign of N
    long iters = 0;
    while (N > 1)
    {
      if (IsEven(N)) N /= 2;
      else N = 3*N+1;
      ++iters;
    }
    return iters;
  }


  void program()
  {
    GlobalManager CoCoAFoundations(UseGMPAllocator);

    BigInt Nmax;
    Nmax = 999;
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
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/tests/test-GMPAllocator.C,v 1.11 2022/02/11 09:53:56 abbott Exp $
// $Log: test-GMPAllocator.C,v $
// Revision 1.11  2022/02/11 09:53:56  abbott
// Summary: Updated copyright (redmine 855)
//
// Revision 1.10  2018/05/18 12:29:29  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.9  2014/10/30 15:12:56  abbott
// Summary: Moved code into namespace CoCoA (redmine #642); improved test-IdealOfPoints, test-IsInteger
// Author: JAA
//
// Revision 1.8  2014/04/30 16:30:56  abbott
// Summary: Replaced size_t by long
// Author: JAA
//
// Revision 1.7  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.6  2011/08/23 06:40:31  bigatti
// -- fixed with new name "BigInt" for old "ZZ"
//
// Revision 1.5  2010/12/17 16:06:26  abbott
// Ensured that all i/o is on standard C++ streams (instead of GlobalInput, etc)
//
// Revision 1.4  2010/11/17 15:53:27  abbott
// Removed out-of-date include of GMPAllocator.H
//
// Revision 1.3  2010/10/22 14:06:56  abbott
// Revised following change of syntax for specifying GMPAllocator.
//
// Revision 1.2  2007/10/30 17:14:06  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:12  abbott
// Imported files
//
// Revision 1.2  2007/03/07 11:34:44  bigatti
// -- minimized #include's
//
// Revision 1.1  2007/03/06 12:34:54  bigatti
// -- added test for GMPAllocator (problems before/after GlobalManager)
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
