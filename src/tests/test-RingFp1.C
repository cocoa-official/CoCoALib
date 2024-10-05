//   Copyright (c)  2007  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/NumTheory-prime.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/RingFpDouble.H"
#include "CoCoA/RingFpLog.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/convert.H"
#include "CoCoA/ring.H"
#include "CoCoA/utils.H"

#include <limits>
using std::numeric_limits;
#include <iostream>
using std::cerr;
using std::endl;
#include<vector>
using std::vector;

namespace CoCoA
{

  void test(const QuotientRing& Fp)
  {
    unsigned long p;
    CoCoA_ASSERT_ALWAYS(IsConvertible(p, characteristic(Fp)));
    CoCoA_ASSERT_ALWAYS(p > 0);
    RingElem x(Fp);
    CoCoA_ASSERT_ALWAYS(x == p);

    // We shall do some arithmetic spot checks; use Fibonacci numbers as "random" numbers.
    constexpr long NumTrials = 100;
    vector<BigInt> fibonacci(NumTrials);
    fibonacci[0] = 0;
    fibonacci[1] = 1;
    for (long n=2; n < NumTrials; ++n)
      fibonacci[n] = fibonacci[n-1] + fibonacci[n-2];

    // This loop checks Fermat's little theorem: x^(p-1) = 1 mod p.
    // For speed check only for x values in the array fibonacci.
    const BigInt BIG(numeric_limits<unsigned long>::max());
    for (long i=1; i < NumTrials; ++i)
    {
      RingElem x(Fp, fibonacci[i]);
      if (IsZero(x)) continue;
      CoCoA_ASSERT_ALWAYS(power(x, p-1) == 1);
      if (p==2) continue;
      RingElem powx = power(x, (p-1)/2);
      CoCoA_ASSERT_ALWAYS(powx == 1 || powx == -1);
      powx = power(x, BIG*(p-1)); // a large power, also multiple of p-1
      CoCoA_ASSERT_ALWAYS(powx == 1);
    }

    // This variable will be handy in the loops below.
    BigInt preimage;
    // Do some spot checks on multiplication.
    for (long i=0; i < NumTrials; ++i)
    {
      for (long j=0; j < NumTrials; ++j)
      {
        const BigInt prod = fibonacci[i]*fibonacci[j];
        RingElem ModProd = RingElem(Fp,fibonacci[i]) * RingElem(Fp,fibonacci[j]);
        CoCoA_ASSERT_ALWAYS(IsInteger(preimage, ModProd));
        CoCoA_ASSERT_ALWAYS(IsZero((prod-preimage)%p));
      }
    }

    // Do some spot checks on division.
    for (long i=0; i < NumTrials; ++i)
    {
      for (long j=0; j < NumTrials; ++j)
      {
        RingElem denom(Fp, fibonacci[j]);
        if (IsZero(denom)) continue;
        RingElem ModQuot = RingElem(Fp,fibonacci[i])/denom;
        CoCoA_ASSERT_ALWAYS(IsInteger(preimage, ModQuot));
        CoCoA_ASSERT_ALWAYS(IsZero((fibonacci[i] - fibonacci[j]*preimage)%p));
      }
    }
  }


  void program()
  {
    using std::floor;
    GlobalManager CoCoAFoundations;

    const ring ZZ = RingZZ();

    vector<unsigned long> p;
    p.push_back(2);

    constexpr double MAX = static_cast<double>(numeric_limits<unsigned long>::max());
    for (double n = 16.0; n <= MAX; n *= 2)
    {
      const long l = ConvertTo<long>(floor(sqrt(n))-1); // will surely succeed
      if (IsPrime(l)) p.push_back(l);
      else p.push_back(PrevPrime(l));
      if (NextPrime(l) == 0) break;
      p.push_back(NextPrime(l));
    }

    for (long i=0; i < len(p); ++i)
    {
      if (IsGoodForRingFp(p[i])) test(NewRingFp(p[i]));
      if (IsGoodForRingFpLog(p[i])) test(NewRingFpLog(p[i]));
      if (IsGoodForRingFpDouble(p[i])) test(NewRingFpDouble(p[i]));
      test(NewQuotientRing(ZZ, ideal(RingElem(ZZ,p[i]))));
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
