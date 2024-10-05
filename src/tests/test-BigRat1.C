//   Copyright (c)  2009  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/BigRatOps.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/error.H"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

namespace CoCoA
{

  // test BigRat ctor, basic arith, comparisons, num & den
  void program()
  {
    GlobalManager CoCoAFoundations;

    constexpr int a=13;
    constexpr int b=21;

    cout << 101 << " = " << BigRat(101) << endl;
    cout << power(101, 9) << " = " << BigRat(power(101,9)) << endl;

    // Test the ctors.
    BigRat frac(a,b);
    CoCoA_ASSERT_ALWAYS(frac == BigRat(a,BigInt(b)));
    CoCoA_ASSERT_ALWAYS(frac == BigRat(BigInt(a),b));
    CoCoA_ASSERT_ALWAYS(frac == BigRat(BigInt(a),BigInt(b)));

    // Generators are negative
    CoCoA_ASSERT_ALWAYS(frac == BigRat(-a,-b));
    CoCoA_ASSERT_ALWAYS(frac == BigRat(-a,BigInt(-b)));
    CoCoA_ASSERT_ALWAYS(frac == BigRat(BigInt(-a),-b));
    CoCoA_ASSERT_ALWAYS(frac == BigRat(BigInt(-a),BigInt(-b)));

    // Generators are not coprime, but positive.
    CoCoA_ASSERT_ALWAYS(frac == BigRat(a*b, b*b));
    CoCoA_ASSERT_ALWAYS(frac == BigRat(a*b, BigInt(b*b)));
    CoCoA_ASSERT_ALWAYS(frac == BigRat(BigInt(a*b), b*b));
    CoCoA_ASSERT_ALWAYS(frac == BigRat(BigInt(a*b), BigInt(b*b)));

    // Generators are not coprime, and negative.
    CoCoA_ASSERT_ALWAYS(frac == BigRat(-a*b, -b*b));
    CoCoA_ASSERT_ALWAYS(frac == BigRat(-a*b, BigInt(-b*b)));
    CoCoA_ASSERT_ALWAYS(frac == BigRat(BigInt(-a*b), -b*b));
    CoCoA_ASSERT_ALWAYS(frac == BigRat(BigInt(-a*b), BigInt(-b*b)));

    // Now check that a/b - a/b gives zero.
    CoCoA_ASSERT_ALWAYS(IsZero(frac-frac));
    CoCoA_ASSERT_ALWAYS(IsZero(frac + BigRat(-a,b)));
    CoCoA_ASSERT_ALWAYS(IsZero(frac + BigRat(a,-b)));
    CoCoA_ASSERT_ALWAYS(IsZero(frac + BigRat(-a,BigInt(b))));
    CoCoA_ASSERT_ALWAYS(IsZero(frac + BigRat(a,BigInt(-b))));
    CoCoA_ASSERT_ALWAYS(IsZero(frac + BigRat(BigInt(-a),b)));
    CoCoA_ASSERT_ALWAYS(IsZero(frac + BigRat(BigInt(a),-b)));
    CoCoA_ASSERT_ALWAYS(IsZero(frac + BigRat(BigInt(-a),BigInt(b))));
    CoCoA_ASSERT_ALWAYS(IsZero(frac + BigRat(BigInt(a),BigInt(-b))));

    // Non trivial GCD between the ctor args.
    CoCoA_ASSERT_ALWAYS(IsZero(frac + BigRat(-a*b,b*b)));
    CoCoA_ASSERT_ALWAYS(IsZero(frac + BigRat(a*b,-b*b)));
    CoCoA_ASSERT_ALWAYS(IsZero(frac + BigRat(-a*b,BigInt(b*b))));
    CoCoA_ASSERT_ALWAYS(IsZero(frac + BigRat(a*b,BigInt(-b*b))));
    CoCoA_ASSERT_ALWAYS(IsZero(frac + BigRat(BigInt(-a*b),b*b)));
    CoCoA_ASSERT_ALWAYS(IsZero(frac + BigRat(BigInt(a*b),-b*b)));
    CoCoA_ASSERT_ALWAYS(IsZero(frac + BigRat(BigInt(-a*b),BigInt(b*b))));
    CoCoA_ASSERT_ALWAYS(IsZero(frac + BigRat(BigInt(a*b),BigInt(-b*b))));


    BigRat frac2;
    CoCoA_ASSERT_ALWAYS(IsZero(frac2));
    frac2 += frac;
    CoCoA_ASSERT_ALWAYS(frac2 == frac);
    frac2 -= frac;
    CoCoA_ASSERT_ALWAYS(IsZero(frac2));
    frac2 = frac;
    frac2 *= frac;
    CoCoA_ASSERT_ALWAYS(frac2 == BigRat(a*a, b*b));
    frac2 /= frac2 +0; // "+0" to avoid compiler warning about self-assignment
    CoCoA_ASSERT_ALWAYS(IsOne(frac2));

    frac2 += BigInt(99);
    CoCoA_ASSERT_ALWAYS(frac2 == 100);
    frac2 -= BigInt(-100);
    CoCoA_ASSERT_ALWAYS(frac2 == 200);
    frac2 *= BigInt(-5);
    CoCoA_ASSERT_ALWAYS(frac2 == -1000);
    frac2 /= BigInt(1000);
    CoCoA_ASSERT_ALWAYS(IsMinusOne(frac2));

    frac2 += 11;
    CoCoA_ASSERT_ALWAYS(frac2 == 10);
    frac2 -= -6;
    CoCoA_ASSERT_ALWAYS(frac2 == 16);
    frac2 *= -5;
    CoCoA_ASSERT_ALWAYS(frac2 == -80);
    frac2 /= 9;
    CoCoA_ASSERT_ALWAYS(frac2 == BigRat(80,-9));

    CoCoA_ASSERT_ALWAYS(BigRat(-71,9) == ++frac2);
    CoCoA_ASSERT_ALWAYS(frac2++ == BigRat(71,-9));
    CoCoA_ASSERT_ALWAYS(frac2 == BigRat(62,-9));
    CoCoA_ASSERT_ALWAYS(--frac2 == BigRat(-71,9));
    CoCoA_ASSERT_ALWAYS(BigRat(71,-9) == frac2--);
    CoCoA_ASSERT_ALWAYS(frac2 == BigRat(-80,9));

    swap(frac, frac2);
    CoCoA_ASSERT_ALWAYS(frac == BigRat(-80,9));
    CoCoA_ASSERT_ALWAYS(frac2 == BigRat(a,b));

    CoCoA_ASSERT_ALWAYS(num(frac) == -80);
    CoCoA_ASSERT_ALWAYS(den(frac) == 9);
    CoCoA_ASSERT_ALWAYS(num(abs(frac)) == 80);
    CoCoA_ASSERT_ALWAYS(den(abs(frac)) == 9);
    CoCoA_ASSERT_ALWAYS(num(-frac) == 80);
    CoCoA_ASSERT_ALWAYS(den(-frac) == 9);

    CoCoA_ASSERT_ALWAYS(frac+frac2 == BigRat(9*a-80*b,9*b));
    CoCoA_ASSERT_ALWAYS(frac-frac2 == BigRat(-9*a-80*b,9*b));
    CoCoA_ASSERT_ALWAYS(frac*frac2 == BigRat(-80*a,9*b));
    CoCoA_ASSERT_ALWAYS(frac/frac2 == BigRat(-80*b,9*a));

    CoCoA_ASSERT_ALWAYS(frac+BigInt(10) == BigRat(10,9));
    CoCoA_ASSERT_ALWAYS(frac-BigInt(10) == BigRat(-170,9));
    CoCoA_ASSERT_ALWAYS(frac*BigInt(10) == BigRat(-800,9));
    CoCoA_ASSERT_ALWAYS(frac/BigInt(10) == BigRat(-8,9));

    CoCoA_ASSERT_ALWAYS(BigInt(-10)+frac == BigRat(-170,9));
    CoCoA_ASSERT_ALWAYS(BigInt(-10)-frac == BigRat(-10,9));
    CoCoA_ASSERT_ALWAYS(BigInt(-10)*frac == BigRat(800,9));
    CoCoA_ASSERT_ALWAYS(BigInt(-10)/frac == BigRat(9,8));

    CoCoA_ASSERT_ALWAYS(frac+10 == BigRat(10,9));
    CoCoA_ASSERT_ALWAYS(frac-10 == BigRat(-170,9));
    CoCoA_ASSERT_ALWAYS(frac*10 == BigRat(-800,9));
    CoCoA_ASSERT_ALWAYS(frac/10 == BigRat(-8,9));

    CoCoA_ASSERT_ALWAYS(10+frac == BigRat(10,9));
    CoCoA_ASSERT_ALWAYS(10-frac == BigRat(170,9));
    CoCoA_ASSERT_ALWAYS(10*frac == BigRat(-800,9));
    CoCoA_ASSERT_ALWAYS(10/frac == BigRat(-9,8));

    BigRat frac3 = power(frac2, 17);
    CoCoA_ASSERT_ALWAYS(frac3 == BigRat(power(a,17), power(b,17)));

    frac3 = power(frac2, BigInt(255));
    CoCoA_ASSERT_ALWAYS(frac3 == BigRat(power(a,255), power(b,255)));

    CoCoA_ASSERT_ALWAYS(!IsZero(frac));
    CoCoA_ASSERT_ALWAYS(!IsOne(frac));
    CoCoA_ASSERT_ALWAYS(IsOne(frac/frac));
    CoCoA_ASSERT_ALWAYS(!IsMinusOne(frac/frac));
    CoCoA_ASSERT_ALWAYS(IsMinusOne(frac*((-1)/frac)));
    CoCoA_ASSERT_ALWAYS(sign(frac) == -1 && sign(frac2) == 1);

    // Comparisons
    CoCoA_ASSERT_ALWAYS(cmp(frac, frac) == 0);
    CoCoA_ASSERT_ALWAYS(cmp(frac, -frac) < 0);
    CoCoA_ASSERT_ALWAYS(cmp(-frac, frac) > 0);

    CoCoA_ASSERT_ALWAYS(cmp(frac, -BigInt(9)) > 0);
    CoCoA_ASSERT_ALWAYS(cmp(frac, -BigInt(8)) < 0);
    CoCoA_ASSERT_ALWAYS(cmp(BigInt(9), -frac) > 0);
    CoCoA_ASSERT_ALWAYS(cmp(BigInt(8), -frac) < 0);

    CoCoA_ASSERT_ALWAYS(cmp(frac, -9) > 0);
    CoCoA_ASSERT_ALWAYS(cmp(frac, -8) < 0);
    CoCoA_ASSERT_ALWAYS(cmp(9, -frac) > 0);
    CoCoA_ASSERT_ALWAYS(cmp(8, -frac) < 0);

    CoCoA_ASSERT_ALWAYS(frac == frac);
    CoCoA_ASSERT_ALWAYS(!(frac != frac));
    CoCoA_ASSERT_ALWAYS(!(frac == frac2));
    CoCoA_ASSERT_ALWAYS(frac != frac2);
    CoCoA_ASSERT_ALWAYS(!(frac < frac));
    CoCoA_ASSERT_ALWAYS(frac < frac2);
    CoCoA_ASSERT_ALWAYS(frac <= frac);
    CoCoA_ASSERT_ALWAYS(frac <= frac2);
    CoCoA_ASSERT_ALWAYS(!(frac > frac));
    CoCoA_ASSERT_ALWAYS(frac2 > frac);
    CoCoA_ASSERT_ALWAYS(frac >= frac);
    CoCoA_ASSERT_ALWAYS(frac2 >= frac);

    CoCoA_ASSERT_ALWAYS(!(frac == BigInt(-9)));
    CoCoA_ASSERT_ALWAYS(frac != BigInt(-9));
    CoCoA_ASSERT_ALWAYS(frac > BigInt(-9));
    CoCoA_ASSERT_ALWAYS(frac >= BigInt(-9));
    CoCoA_ASSERT_ALWAYS(frac < BigInt(-8));
    CoCoA_ASSERT_ALWAYS(frac <= BigInt(-8));

    CoCoA_ASSERT_ALWAYS(!(BigInt(-9) == frac));
    CoCoA_ASSERT_ALWAYS(BigInt(-9) != frac);
    CoCoA_ASSERT_ALWAYS(BigInt(-8) > frac);
    CoCoA_ASSERT_ALWAYS(BigInt(-8) >= frac);
    CoCoA_ASSERT_ALWAYS(-BigInt(9) < frac);
    CoCoA_ASSERT_ALWAYS(-BigInt(9) <= frac);

    CoCoA_ASSERT_ALWAYS(!(frac == -9));
    CoCoA_ASSERT_ALWAYS(frac != -9);
    CoCoA_ASSERT_ALWAYS(frac > -9);
    CoCoA_ASSERT_ALWAYS(frac >= -9);
    CoCoA_ASSERT_ALWAYS(frac < -8);
    CoCoA_ASSERT_ALWAYS(frac <= -8);

    CoCoA_ASSERT_ALWAYS(!(-9 == frac));
    CoCoA_ASSERT_ALWAYS(-9 != frac);
    CoCoA_ASSERT_ALWAYS(-8 > frac);
    CoCoA_ASSERT_ALWAYS(-8 >= frac);
    CoCoA_ASSERT_ALWAYS(-9 < frac);
    CoCoA_ASSERT_ALWAYS(-9 <= frac);

    CoCoA_ASSERT_ALWAYS(round(frac) == -9);
    CoCoA_ASSERT_ALWAYS(round(BigRat(1,2)) == RoundDiv(1,2));
    CoCoA_ASSERT_ALWAYS(round(BigRat(-1,2)) == RoundDiv(-1,2));
  }

}


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
