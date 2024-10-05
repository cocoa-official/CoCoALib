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


#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/NumTheory-prime.H"
#include "CoCoA/SmallFpImpl.H"
#include "CoCoA/error.H"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;


namespace CoCoA
{

  // Check that NonRedValue cannot overflow (within limit of myMaxIters)
  void OverFlowTest(long p)
  {
    SmallFpImpl ModP(p);

    SmallFpImpl::value p1=ModP.myReduce(p-1);
    const long limit = ModP.myMaxIters();
    for (long InitVal = 0; InitVal < p; InitVal += 3) // stepsize=3 to be faster
    {
      SmallFpImpl::NonRedValue ans=ModP.myReduce(InitVal);
      long IterCount=0;
      for (long i=0; i<p; ++i)
      {
        ans += p1*p1; // p1 is auto converted to NonRedValue here!
        ++IterCount;
        if (IterCount < limit) continue;
        ans = ModP.myHalfNormalize(ans);
        IterCount = 0;
      }
      CoCoA_ASSERT_ALWAYS(ModP.myNormalize(ans) == ModP.myReduce(InitVal));
    }
  }

  void program()
  {
    GlobalManager CoCoAFoundations;

    for (int p=NextPrime(4000); p < 33000; p = NextPrime(p+7100))
      OverFlowTest(p);
  }

}  // end of namespace CoCoA

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
