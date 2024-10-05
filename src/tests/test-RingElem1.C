//   Copyright (c)  2012  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/PolyRing.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/error.H"
#include "CoCoA/ring.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"


#include <iostream>
using std::cerr;
using std::endl;



namespace CoCoA
{

  void TestGCD(ring R)
  {
    if (!IsTrueGCDDomain(R)) return;
    CoCoA_ASSERT_ALWAYS(IsZero(gcd(zero(R),zero(R))));
    CoCoA_ASSERT_ALWAYS(IsOne(gcd(zero(R),one(R))));
    CoCoA_ASSERT_ALWAYS(IsOne(gcd(0,one(R))));
    CoCoA_ASSERT_ALWAYS(IsOne(gcd(zero(R),1)));
    CoCoA_ASSERT_ALWAYS(IsOne(gcd(one(R),zero(R))));
    CoCoA_ASSERT_ALWAYS(IsOne(gcd(1,zero(R))));
    CoCoA_ASSERT_ALWAYS(IsOne(gcd(one(R),0)));
    CoCoA_ASSERT_ALWAYS(IsOne(gcd(one(R),one(R))));
    CoCoA_ASSERT_ALWAYS(IsOne(gcd(1,one(R))));
    CoCoA_ASSERT_ALWAYS(IsOne(gcd(one(R),1)));

    RingElem two(R,2);
    if (!IsZero(two))
    {
      CoCoA_ASSERT_ALWAYS(IsInvertible(gcd(zero(R),two)/two));
      CoCoA_ASSERT_ALWAYS(IsInvertible(gcd(0,two)/two));
      CoCoA_ASSERT_ALWAYS(IsInvertible(gcd(zero(R),2)/two));
      CoCoA_ASSERT_ALWAYS(IsInvertible(gcd(two,two)/two));
      CoCoA_ASSERT_ALWAYS(IsInvertible(gcd(2,two)/two));
      CoCoA_ASSERT_ALWAYS(IsInvertible(gcd(two,2)/two));
      CoCoA_ASSERT_ALWAYS(IsInvertible(gcd(two*two,two)/two));
      CoCoA_ASSERT_ALWAYS(IsInvertible(gcd(4,two)/two));
      CoCoA_ASSERT_ALWAYS(IsInvertible(gcd(two*two,2)/two));
    }

    if (!IsPolyRing(R)) return;
    const RingElem x = indet(R,0);
    CoCoA_ASSERT_ALWAYS(IsInvertible(gcd(x,x)/x));
    CoCoA_ASSERT_ALWAYS(IsInvertible(gcd(0,x)/x));
    CoCoA_ASSERT_ALWAYS(IsInvertible(gcd(x,0)/x));
  }

  void TestLCM(ring R)
  {
    if (!IsTrueGCDDomain(R)) return;
    CoCoA_ASSERT_ALWAYS(IsZero(lcm(zero(R),zero(R))));
    CoCoA_ASSERT_ALWAYS(IsZero(lcm(zero(R),one(R))));
    CoCoA_ASSERT_ALWAYS(IsZero(lcm(0,one(R))));
    CoCoA_ASSERT_ALWAYS(IsZero(lcm(zero(R),1)));
    CoCoA_ASSERT_ALWAYS(IsZero(lcm(one(R),zero(R))));
    CoCoA_ASSERT_ALWAYS(IsZero(lcm(1,zero(R))));
    CoCoA_ASSERT_ALWAYS(IsZero(lcm(one(R),0)));
    CoCoA_ASSERT_ALWAYS(IsOne(lcm(one(R),one(R))));
    CoCoA_ASSERT_ALWAYS(IsOne(lcm(1,one(R))));
    CoCoA_ASSERT_ALWAYS(IsOne(lcm(one(R),1)));

    RingElem two(R,2);
    if (!IsZero(two))
    {
      CoCoA_ASSERT_ALWAYS(IsZero(lcm(zero(R),two)));
      CoCoA_ASSERT_ALWAYS(IsZero(lcm(0,two)));
      CoCoA_ASSERT_ALWAYS(IsZero(lcm(zero(R),2)));
      CoCoA_ASSERT_ALWAYS(IsInvertible(lcm(two,two)/two));
      CoCoA_ASSERT_ALWAYS(IsInvertible(lcm(2,two)/two));
      CoCoA_ASSERT_ALWAYS(IsInvertible(lcm(two,2)/two));
      CoCoA_ASSERT_ALWAYS(IsInvertible(lcm(two*two,two)/4));
      CoCoA_ASSERT_ALWAYS(IsInvertible(lcm(4,two)/4));
      CoCoA_ASSERT_ALWAYS(IsInvertible(lcm(two*two,2)/4));
    }

    if (!IsPolyRing(R)) return;
    const RingElem x = indet(R,0);
    CoCoA_ASSERT_ALWAYS(IsInvertible(lcm(x,x)/x));
    CoCoA_ASSERT_ALWAYS(IsZero(lcm(0,x)));
    CoCoA_ASSERT_ALWAYS(IsZero(lcm(x,0)));
  }


  void TestIsCoprime(ring R)
  {
    if (!IsTrueGCDDomain(R)) return;
    if (IsPolyRing(R))  CoCoA_ASSERT_ALWAYS(IsField(CoeffRing(R)));
    CoCoA_ASSERT_ALWAYS(IsCoprime(zero(R),one(R)));
    CoCoA_ASSERT_ALWAYS(IsCoprime(0,one(R)));
    CoCoA_ASSERT_ALWAYS(IsCoprime(zero(R),1));
    CoCoA_ASSERT_ALWAYS(IsCoprime(one(R),zero(R)));
    CoCoA_ASSERT_ALWAYS(IsCoprime(1,zero(R)));
    CoCoA_ASSERT_ALWAYS(IsCoprime(one(R),0));
    CoCoA_ASSERT_ALWAYS(IsCoprime(one(R),one(R)));
    CoCoA_ASSERT_ALWAYS(IsCoprime(one(R),-one(R)));
    CoCoA_ASSERT_ALWAYS(IsCoprime(-one(R),one(R)));
    CoCoA_ASSERT_ALWAYS(IsCoprime(-one(R),-one(R)));
    CoCoA_ASSERT_ALWAYS(IsCoprime(1,one(R)));
    CoCoA_ASSERT_ALWAYS(IsCoprime(1,-one(R)));
    CoCoA_ASSERT_ALWAYS(IsCoprime(-1,one(R)));
    CoCoA_ASSERT_ALWAYS(IsCoprime(-1,-one(R)));
    CoCoA_ASSERT_ALWAYS(IsCoprime(one(R),1));
    CoCoA_ASSERT_ALWAYS(IsCoprime(one(R),-1));
    CoCoA_ASSERT_ALWAYS(IsCoprime(-one(R),1));
    CoCoA_ASSERT_ALWAYS(IsCoprime(-one(R),-1));

    RingElem two(R,2);
    if (IsInvertible(two))
    {
      CoCoA_ASSERT_ALWAYS(IsCoprime(zero(R),two));
      CoCoA_ASSERT_ALWAYS(IsCoprime(two,zero(R)));
      CoCoA_ASSERT_ALWAYS(IsCoprime(zero(R),-two));
      CoCoA_ASSERT_ALWAYS(IsCoprime(-two,zero(R)));
      CoCoA_ASSERT_ALWAYS(IsCoprime(zero(R),2));
      CoCoA_ASSERT_ALWAYS(IsCoprime(2,zero(R)));
      CoCoA_ASSERT_ALWAYS(IsCoprime(zero(R),-2));
      CoCoA_ASSERT_ALWAYS(IsCoprime(-2,zero(R)));
    }

    if (!IsPolyRing(R)) return;
    const RingElem x = indet(R,0);
    CoCoA_ASSERT_ALWAYS(IsCoprime(x, one(R)));
    CoCoA_ASSERT_ALWAYS(IsCoprime(1,x));
    CoCoA_ASSERT_ALWAYS(IsCoprime(x,1));
  }

  void program()
  {
    // Test GCD and LCM in various rings.
    GlobalManager CoCoAFoundations;

    ring R1 = RingZZ();
    TestGCD(R1);   TestLCM(R1);  TestIsCoprime(R1);

    ring R2 = RingQQ();
    TestGCD(R2);   TestLCM(R2);

    ring R3 = NewPolyRing(RingQQ(), symbols("x"));
    TestGCD(R3);   TestLCM(R3);  TestIsCoprime(R3);

    ring R4 = NewPolyRing(RingQQ(), symbols("x,y"));
    TestGCD(R4);   TestLCM(R4);  TestIsCoprime(R4);

    ring R5 = NewPolyRing(RingZZ(), symbols("a"));
    TestGCD(R5);   TestLCM(R5);

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
