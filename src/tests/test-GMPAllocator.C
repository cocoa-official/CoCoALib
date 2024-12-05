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
