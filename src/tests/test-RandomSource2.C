//   Copyright (c)  2010  John Abbott,  Anna M. Bigatti

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


#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/error.H"
#include "CoCoA/random.H"

#include <iostream>
using std::cerr;
using std::endl;
#include <limits>
using std::numeric_limits;
#include <vector>
using std::vector;


namespace CoCoA
{

  void program()
  {
    // This test checks the behaviour of the RandomLong function which
    // generates pseudo-random machine integer values from a RandomSource.
    GlobalManager CoCoAFoundations;

    RandomSource RndSrc;

    // Check that the arg checks are working, and that generated values
    // are within the specified range.
    for (int i=-9; i < 10; ++i)
      for (int j=-9; j < 10; ++j)
        for (int k=0; k < 100; ++k)
        {
          try
          {
            const long rnd = RandomLong(RndSrc, i, j);
            CoCoA_ASSERT_ALWAYS(i <= j && i <= rnd && rnd <= j);
          }
          catch (const ErrorInfo& err) { CoCoA_ASSERT_ALWAYS(err == ERR::ReqNonEmpty); }
        }


    {
      // Check that we can generate the full range of longs
      constexpr int NumIters = 2147000;
      constexpr long lwb = numeric_limits<long>::min();
      constexpr long upb = numeric_limits<long>::max();
      long smallest = 0;
      long largest = 0;
      for (int i=0; i < NumIters; ++i)
      {
        const long rnd = RandomLong(RndSrc, lwb, upb);
        if (rnd > largest)  largest = rnd;
        else if (rnd < smallest)  smallest = rnd;
      }
      CoCoA_ASSERT_ALWAYS(smallest < lwb + 14*(upb/NumIters)); // 99.9% chance of being true
      CoCoA_ASSERT_ALWAYS(largest > upb - 14*(upb/NumIters));  // 99.9% chance of being true
    }


    // Generate some histograms, and check that they are fairly uniform.
    constexpr long NumCases = 10;
    constexpr long ExpectedFreq = 10000;
    constexpr long NumTrials = ExpectedFreq*NumCases;
    for (int lwb=-12345; lwb < 23456; lwb += 321)
    {
      vector<long> hist(NumCases);
      for (long i=0; i < NumTrials; ++i)
        ++hist[RandomLong(RndSrc, lwb, lwb+NumCases-1) - lwb];

      int min = ExpectedFreq;
      int max = 0;
      for (int i=0; i < NumCases; ++i)
      {
        if (hist[i] > max)  max = hist[i];
        else if (hist[i] < min)  min = hist[i];
      }
      CoCoA_ASSERT_ALWAYS(min > 0.8*max);
    }

  }

} // end of namespace CoCoA


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
    cerr << "***ERROR***  UNCAUGHT CoCoA Error";
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
