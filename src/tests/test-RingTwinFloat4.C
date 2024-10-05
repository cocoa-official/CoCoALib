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


#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/RingTwinFloat.H"


#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <algorithm>
using std::swap;



namespace CoCoA
{

  RingElem SquareRoot(RingElem x)
  {
    if (!IsRingTwinFloat(owner(x)))
      CoCoA_THROW_ERROR("Argument must be element of RingTwinFloat", "SquareRoot");
    ring RR = owner(x);
    RingElem approx(RR);
    RingElem approx2(RR, 1);

//  while (approx != approx2)
    while (!IsPracticallyEqual(approx, approx2))
    {
      approx = (approx2 + x/approx2)/2;
      swap(approx, approx2);
    }
    return (approx2 + x/approx2)/2;
  }


  void trial(long BitPrec)
  {
    const RingTwinFloat RR = NewRingTwinFloat(BitPrec);
    CoCoA_ASSERT_ALWAYS(PrecisionBits(RR) >= BitPrec);

    // cout << "First trial: loop comparing x+3 with (x^2-9)/(x-3) for\n";
    // cout << "x = 3 + eps where  eps = 1/(2^n) for  n=1,2,...200" << endl;
    const RingElem three(RR, 3);
    RingElem eps(RR, 1);
    try
    {
      for (long i = 1; i < 200; ++i)
      {
        eps /= 2;
        const RingElem x = three + eps;
        const RingElem v1 = x + 3;
        const RingElem v2 = (x*x-9)/(x-3);
        CoCoA_ASSERT_ALWAYS(v1 == v2);
      }
      CoCoA_ASSERT_ALWAYS(BitPrec >= 200);
    }
    catch (const RingTwinFloat::InsufficientPrecision&)
    {
      CoCoA_ASSERT_ALWAYS(BitPrec <= 200);
    }


    // cout << "Second trial: see whether (sqrt(2)-1)^6 = 99-70*sqrt(2)\n";
    try
    {
      const RingElem sqrt2 = SquareRoot(RingElem(RR,2));
      CoCoA_ASSERT_ALWAYS(sqrt2*sqrt2 == 2);
      CoCoA_ASSERT_ALWAYS(power(sqrt2-1, 6) == 99-70*sqrt2);
    }
    catch (const RingTwinFloat::InsufficientPrecision&)
    {
      // We will never get here unless the approx sqrt(2) ~ 99/70 is replaced by
      // a much finer one.
    }


    // cout << "Third trial: an almost equality between square roots:\n";
    // cout << "sqrt(176)+sqrt(195)+sqrt(2025) =?= sqrt(190)+sqrt(398)+sqrt(1482)" << endl;
    const RingElem sum1 = SquareRoot(RingElem(RR,176))+SquareRoot(RingElem(RR,195))+SquareRoot(RingElem(RR,2025));
    const RingElem sum2 = SquareRoot(RingElem(RR,190))+SquareRoot(RingElem(RR,398))+SquareRoot(RingElem(RR,1482));

    try
    {
      // The comparison  sum1 != sum2  can produce 3 possible outcomes:
      // with low precision (<= 27 bits) the values may be seen as equal
      // with high precision (>= 35 bits) the values may be seen as unequal
      // while from 23 to 35 bits (incl), the comparison may fail.
      if (sum1 == sum2) CoCoA_ASSERT_ALWAYS(BitPrec <= 27);
      else CoCoA_ASSERT_ALWAYS(BitPrec >= 31);
    }
    catch (const RingTwinFloat::InsufficientPrecision&)
    {
      CoCoA_ASSERT_ALWAYS((BitPrec >= 22) && (BitPrec <= 35));
    }
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    for (long prec=10; prec <= 220; ++prec)
      trial(prec);
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
