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


#include "CoCoA/BigIntOps.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/error.H"
#include "CoCoA/random.H"

#include <iostream>
using std::cerr;
using std::endl;


namespace CoCoA
{

  void program()
  {
    // This test checks the generation of pseudo-random large integers from a
    // RandomSource object.  The test is really very mild; I'm not even sure
    // how to write a thorough test.
    GlobalManager CoCoAFoundations;

    RandomSource RndSrc;

    for (int i=-9; i < 10; ++i)
      for (int j=-9; j < 10; ++j)
        for (int k=0; k < 100; ++k)
        {
          try
          {
            const BigInt rnd = RandomBigInt(RndSrc, i, j);
            CoCoA_ASSERT_ALWAYS(i <= j && i <= rnd && rnd <= j);
          }
          catch (const ErrorInfo& err) { CoCoA_ASSERT_ALWAYS(err == ERR::BadArg); }
        }

    // Now choose a random (negative) range, and check that even numbers
    // are chosen about half the time, and that about half the values
    // produced lie in the central interval.
    const BigInt upb = RandomBigInt(RndSrc, -power(2,100), -1);
    const BigInt lwb = RandomBigInt(RndSrc, 2*upb, upb);
    const BigInt width = upb-lwb; // we hope that this will not be too small (e.g. width > 100)
    const BigInt left = lwb+width/4;
    const BigInt right = upb-width/4;

    constexpr int NumTrials = 2000; // if you change this, then change also 0.45 & 0.55 in the asserts after the loop!!
    int CountEven = 0;
    int CountCentral = 0;
    for  (int i=0; i < NumTrials; ++i)
    {
      const BigInt rnd = RandomBigInt(RndSrc, lwb, upb);
      CoCoA_ASSERT_ALWAYS(lwb <= rnd && rnd <= upb);
      if (IsEven(rnd)) ++CountEven;
      if (left <= rnd && rnd <= right) ++CountCentral;
    }
    // These assertions are wrong with prob < 10^(-5)
    CoCoA_ASSERT_ALWAYS( 0.45*NumTrials < CountEven && CountEven <= 0.55*NumTrials );
    CoCoA_ASSERT_ALWAYS( 0.45*NumTrials < CountCentral && CountCentral <= 0.55*NumTrials );
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
