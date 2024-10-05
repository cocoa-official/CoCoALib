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
#include "CoCoA/RingZZ.H"
#include "CoCoA/error.H"
#include "CoCoA/ring.H"

#include <iostream>
using std::cerr;
using std::endl;



namespace CoCoA
{

  void test(const ring& ZZmodN)
  {

    CoCoA_ASSERT_ALWAYS(IsQuotientRing(ZZmodN));
    const QuotientRing QR = ZZmodN;
    CoCoA_ASSERT_ALWAYS(IsZZ(BaseRing(QR)));

    const BigInt p = characteristic(ZZmodN);
    CoCoA_ASSERT_ALWAYS(IsProbPrime(p) == IsField(ZZmodN));
    CoCoA_ASSERT_ALWAYS(IsZero(RingElem(ZZmodN, p)));

    RingElem n(ZZmodN);

    CoCoA_ASSERT_ALWAYS(IsZero(n));
    CoCoA_ASSERT_ALWAYS(!IsInvertible(n));
    CoCoA_ASSERT_ALWAYS(n == 0);
    CoCoA_ASSERT_ALWAYS(0 == n);
    CoCoA_ASSERT_ALWAYS(n == zero(ZZmodN));

    // We cannot do inequalities: (not even in ZZ/(0))
    try
    {
      const bool junk = (n > 0);
      CoCoA_ASSERT_ALWAYS(!"NEVER GET HERE!");
      cerr << junk; // never executed, just to keep compiler quiet!
    }
    catch (const CoCoA::ErrorInfo& err) { CoCoA_ASSERT_ALWAYS(err == ERR::ReqOrdDom); }

    const BigInt N = BigInt(1);
    n = N;
    CoCoA_ASSERT_ALWAYS(IsInvertible(n));
    CoCoA_ASSERT_ALWAYS(IsOne(n));
    CoCoA_ASSERT_ALWAYS(IsMinusOne(-n));
    CoCoA_ASSERT_ALWAYS(n == 1);
    CoCoA_ASSERT_ALWAYS(1 == n);
    CoCoA_ASSERT_ALWAYS(n == N);
    CoCoA_ASSERT_ALWAYS(N == n);
    n = N;
    CoCoA_ASSERT_ALWAYS(IsOne(n));

    // Some simple arithmetic
    n = n+1; CoCoA_ASSERT_ALWAYS(n == 2);
    n = n-1; CoCoA_ASSERT_ALWAYS(n == 1);
    n = n*1; CoCoA_ASSERT_ALWAYS(n == 1);
    n = n/1; CoCoA_ASSERT_ALWAYS(n == 1);
    n = 1+n; CoCoA_ASSERT_ALWAYS(n == 2);
    n = 1-n; CoCoA_ASSERT_ALWAYS(n == -1);
    n = 1*n; CoCoA_ASSERT_ALWAYS(n == -1);
    n = 1/n; CoCoA_ASSERT_ALWAYS(n == -1);
    n = n+n; CoCoA_ASSERT_ALWAYS(n == -2);
    n = n-n; CoCoA_ASSERT_ALWAYS(n == 0);
    n = n*n; CoCoA_ASSERT_ALWAYS(n == 0);
    n = n/(n+1); CoCoA_ASSERT_ALWAYS(n == 0);

    // assignment arithmetic
    n = p-1;
    const RingElem m = n;

    n += m;
    n -= m;
    n *= m;
    n /= m;

    // 2023-03-16 clang++ warns about self assignment (sigh)
    /*TEST SELF ASSIGN*/ n = m; n += n; CoCoA_ASSERT_ALWAYS(n == 2*m);
    /*TEST SELF ASSIGN*/ n = m; n -= n; CoCoA_ASSERT_ALWAYS(IsZero(n));
    /*TEST SELF ASSIGN*/ n = m; n *= n; CoCoA_ASSERT_ALWAYS(n == m*m);
    /*TEST SELF ASSIGN*/ n = m; n /= n; CoCoA_ASSERT_ALWAYS(IsOne(n));
  }




  void program()
  {
    GlobalManager CoCoAFoundations;

    test(NewZZmod(2));
    test(NewZZmod(-2));
    test(NewZZmod(3));
    test(NewZZmod(-3));
    test(NewZZmod(29641));
    test(NewZZmod(-29641));
    test(NewZZmod(10000019));
    test(NewZZmod(-10000019));
    test(NewZZmod(32768)); // not a field

    const ring ZZ = RingZZ();
    ideal I(RingElem(ZZ, 100000007));
    test(NewQuotientRing(ZZ, I));  // field

    I = ideal(power(RingElem(ZZ, 10), 20));
    test(NewQuotientRing(ZZ, I));  // not a field

    // Now check the three forbidden cases.

    // Should this really be forbidden?
    try { test(NewZZmod(0)); CoCoA_ASSERT_ALWAYS(!"NEVER GET HERE!"); }
    catch (const CoCoA::ErrorInfo& err) { CoCoA_ASSERT_ALWAYS(err == ERR::BadQuotRing); }

    try { test(NewZZmod(1)); CoCoA_ASSERT_ALWAYS(!"NEVER GET HERE!"); }
    catch (const CoCoA::ErrorInfo& err) { CoCoA_ASSERT_ALWAYS(err == ERR::BadQuotRing); }

    try { test(NewZZmod(-1)); CoCoA_ASSERT_ALWAYS(!"NEVER GET HERE!"); }
    catch (const CoCoA::ErrorInfo& err) { CoCoA_ASSERT_ALWAYS(err == ERR::BadQuotRing); }
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
