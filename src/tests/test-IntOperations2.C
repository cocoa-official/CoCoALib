//   Copyright (c)  2014  John Abbott,  Anna M. Bigatti

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
using std::cerr;
using std::endl;
#include <limits>
using std::numeric_limits;



namespace CoCoA
{

  // tests for power(zero,zero), IsPowerOf2
  void program()
  {
    // Sundry tests for edge cases for several integer functions.
    GlobalManager CoCoAFoundations;

    const int zero = 0;
    const BigInt ZERO;
    CoCoA_ASSERT_ALWAYS(IsOne(power(zero,zero)));
    CoCoA_ASSERT_ALWAYS(IsOne(power(ZERO,zero)));
    CoCoA_ASSERT_ALWAYS(IsOne(power(zero,ZERO)));
    CoCoA_ASSERT_ALWAYS(IsOne(power(ZERO,ZERO)));

    CoCoA_ASSERT_ALWAYS(SmallPower(zero,zero) == 1);

    const int IntMax = numeric_limits<int>::max();
    const int IntMin = numeric_limits<int>::min();
    const long LongMax = numeric_limits<long>::max();
    const long LongMin = numeric_limits<long>::min();
    const unsigned long ULongMax = numeric_limits<unsigned long>::max();
    const unsigned long ULongMin = numeric_limits<unsigned long>::min();

    CoCoA_ASSERT_ALWAYS(IsPowerOf2(1));
    CoCoA_ASSERT_ALWAYS(IsPowerOf2(2));
    CoCoA_ASSERT_ALWAYS(!IsPowerOf2(3));
    CoCoA_ASSERT_ALWAYS(IsPowerOf2(4));

    CoCoA_ASSERT_ALWAYS(!IsPowerOf2(-1));
    CoCoA_ASSERT_ALWAYS(!IsPowerOf2(-2));
    CoCoA_ASSERT_ALWAYS(!IsPowerOf2(3));
    CoCoA_ASSERT_ALWAYS(!IsPowerOf2(-4));

    CoCoA_ASSERT_ALWAYS(!IsPowerOf2(IntMax));
    CoCoA_ASSERT_ALWAYS(!IsPowerOf2(IntMin));
    CoCoA_ASSERT_ALWAYS(!IsPowerOf2(LongMax));
    CoCoA_ASSERT_ALWAYS(!IsPowerOf2(LongMin));
    CoCoA_ASSERT_ALWAYS(!IsPowerOf2(ULongMax));
    CoCoA_ASSERT_ALWAYS(!IsPowerOf2(ULongMin));

    CoCoA_ASSERT_ALWAYS(!IsPowerOf2(BigInt(0)));
    for (int n=0; n <= 1024; ++n)
    {
      CoCoA_ASSERT_ALWAYS(IsPowerOf2(power(2,n)));
      CoCoA_ASSERT_ALWAYS(!IsPowerOf2(-power(2,n)));
      CoCoA_ASSERT_ALWAYS(n==0 || !IsPowerOf2(1+power(2,n)));
      CoCoA_ASSERT_ALWAYS(n==1 || !IsPowerOf2(-1+power(2,n)));
      CoCoA_ASSERT_ALWAYS(!IsPowerOf2(3*power(2,n)));
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
