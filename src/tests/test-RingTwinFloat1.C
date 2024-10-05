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
#include "CoCoA/BigRat.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/RingTwinFloat.H"
#include "CoCoA/ring.H"


#include <iostream>
using std::cerr;
using std::endl;



namespace CoCoA
{

// ASSUMES argument is a RingElem with integer value.
// Makes some trivial check that the value is equal to itself.
  void tautology(const RingElem& n)
  {
    CoCoA_ASSERT_ALWAYS(IsZero(n) ^ IsInvertible(n)); // RingTwinFloat is supposed to be a field.
    BigInt N;
    CoCoA_ASSERT_ALWAYS(IsInteger(N, n));
    CoCoA_ASSERT_ALWAYS(n == N);
    CoCoA_ASSERT_ALWAYS(N == n);
    RingElem n2(owner(n), N);
    CoCoA_ASSERT_ALWAYS(n == n2);
    CoCoA_ASSERT_ALWAYS(n2 == n);
  }


// ASSUMES that n > m.  Checks various inequalities.
  void greater(const RingElem& n, int m)
  {
    CoCoA_ASSERT_ALWAYS(!(n == m));
    CoCoA_ASSERT_ALWAYS(n != m);
    CoCoA_ASSERT_ALWAYS(n >= m);
    CoCoA_ASSERT_ALWAYS(n > m);
    CoCoA_ASSERT_ALWAYS(!(n <= m));
    CoCoA_ASSERT_ALWAYS(!(n < m));
    CoCoA_ASSERT_ALWAYS(!(m == n));
    CoCoA_ASSERT_ALWAYS(m != n);
    CoCoA_ASSERT_ALWAYS(!(m >= n));
    CoCoA_ASSERT_ALWAYS(!(m > n));
    CoCoA_ASSERT_ALWAYS(m <= n);
    CoCoA_ASSERT_ALWAYS(m < n);
  }

// ASSUMES that n == m.  Checks various inequalities.
  void equal(const RingElem& n, int m)
  {
    CoCoA_ASSERT_ALWAYS(n == m);
    CoCoA_ASSERT_ALWAYS(!(n != m));
    CoCoA_ASSERT_ALWAYS(n >= m);
    CoCoA_ASSERT_ALWAYS(!(n > m));
    CoCoA_ASSERT_ALWAYS(n <= m);
    CoCoA_ASSERT_ALWAYS(!(n < m));
    CoCoA_ASSERT_ALWAYS(m == n);
    CoCoA_ASSERT_ALWAYS(!(m != n));
    CoCoA_ASSERT_ALWAYS(m >= n);
    CoCoA_ASSERT_ALWAYS(!(m > n));
    CoCoA_ASSERT_ALWAYS(m <= n);
    CoCoA_ASSERT_ALWAYS(!(m < n));
  }

// ASSUMES that n < m.  Checks various inequalities.
  void less(const RingElem& n, int m)
  {
    CoCoA_ASSERT_ALWAYS(!(n == m));
    CoCoA_ASSERT_ALWAYS(n != m);
    CoCoA_ASSERT_ALWAYS(!(n >= m));
    CoCoA_ASSERT_ALWAYS(!(n > m));
    CoCoA_ASSERT_ALWAYS(n <= m);
    CoCoA_ASSERT_ALWAYS(n < m);
    CoCoA_ASSERT_ALWAYS(!(m == n));
    CoCoA_ASSERT_ALWAYS(m != n);
    CoCoA_ASSERT_ALWAYS(m >= n);
    CoCoA_ASSERT_ALWAYS(m > n);
    CoCoA_ASSERT_ALWAYS(!(m <= n));
    CoCoA_ASSERT_ALWAYS(!(m < n));
  }


// ASSUMES that n > m.  Checks various inequalities.
  void greater(const RingElem& n, const BigInt& m)
  {
    CoCoA_ASSERT_ALWAYS(!(n == m));
    CoCoA_ASSERT_ALWAYS(n != m);
    CoCoA_ASSERT_ALWAYS(n >= m);
    CoCoA_ASSERT_ALWAYS(n > m);
    CoCoA_ASSERT_ALWAYS(!(n <= m));
    CoCoA_ASSERT_ALWAYS(!(n < m));
    CoCoA_ASSERT_ALWAYS(!(m == n));
    CoCoA_ASSERT_ALWAYS(m != n);
    CoCoA_ASSERT_ALWAYS(!(m >= n));
    CoCoA_ASSERT_ALWAYS(!(m > n));
    CoCoA_ASSERT_ALWAYS(m <= n);
    CoCoA_ASSERT_ALWAYS(m < n);
  }

// ASSUMES that n == m.  Checks various inequalities.
  void equal(const RingElem& n, const BigInt& m)
  {
    CoCoA_ASSERT_ALWAYS(n == m);
    CoCoA_ASSERT_ALWAYS(!(n != m));
    CoCoA_ASSERT_ALWAYS(n >= m);
    CoCoA_ASSERT_ALWAYS(!(n > m));
    CoCoA_ASSERT_ALWAYS(n <= m);
    CoCoA_ASSERT_ALWAYS(!(n < m));
    CoCoA_ASSERT_ALWAYS(m == n);
    CoCoA_ASSERT_ALWAYS(!(m != n));
    CoCoA_ASSERT_ALWAYS(m >= n);
    CoCoA_ASSERT_ALWAYS(!(m > n));
    CoCoA_ASSERT_ALWAYS(m <= n);
    CoCoA_ASSERT_ALWAYS(!(m < n));
  }

// ASSUMES that n < m.  Checks various inequalities.
  void less(const RingElem& n, const BigInt& m)
  {
    CoCoA_ASSERT_ALWAYS(!(n == m));
    CoCoA_ASSERT_ALWAYS(n != m);
    CoCoA_ASSERT_ALWAYS(!(n >= m));
    CoCoA_ASSERT_ALWAYS(!(n > m));
    CoCoA_ASSERT_ALWAYS(n <= m);
    CoCoA_ASSERT_ALWAYS(n < m);
    CoCoA_ASSERT_ALWAYS(!(m == n));
    CoCoA_ASSERT_ALWAYS(m != n);
    CoCoA_ASSERT_ALWAYS(m >= n);
    CoCoA_ASSERT_ALWAYS(m > n);
    CoCoA_ASSERT_ALWAYS(!(m <= n));
    CoCoA_ASSERT_ALWAYS(!(m < n));
  }


  void BasicChecks(ring RR)
  {
    CoCoA_ASSERT_ALWAYS(IsRingTwinFloat(RR));

    RingElem x(RR);

    CoCoA_ASSERT_ALWAYS(IsZero(x));
    CoCoA_ASSERT_ALWAYS(!IsInvertible(x));
    CoCoA_ASSERT_ALWAYS(x == zero(RR));

    equal(x, 0);
    equal(x, BigInt(0));
    greater(x, -1);
    greater(x, BigInt(-1));
    less(x, 1);
    less(x, BigInt(1));

    x = 1;
    CoCoA_ASSERT_ALWAYS(!IsZero(x));
    CoCoA_ASSERT_ALWAYS(IsInvertible(x));
    CoCoA_ASSERT_ALWAYS(IsOne(x));
    CoCoA_ASSERT_ALWAYS(!IsMinusOne(x));
    equal(x, 1);
    equal(x, BigInt(1));
    greater(x, -1);
    greater(x, BigInt(-1));
    less(x, 2);
    less(x, BigInt(2));

    x = -1;
    CoCoA_ASSERT_ALWAYS(!IsZero(x));
    CoCoA_ASSERT_ALWAYS(IsInvertible(x));
    CoCoA_ASSERT_ALWAYS(!IsOne(x));
    CoCoA_ASSERT_ALWAYS(IsMinusOne(x));
    equal(x, -1);
    equal(x, BigInt(-1));
    greater(x, -2);
    greater(x, BigInt(-2));
    less(x, 1);
    less(x, BigInt(1));
  }


// ASSUMES that n & m have integer values, and m != 0.
// Checks that the basic arithmetic ops give the expected results.
  void arithmetic(RingElem n, const RingElem& m)
  {
    BigInt N;
    CoCoA_ASSERT_ALWAYS(IsInteger(N, n));
    BigInt M;
    CoCoA_ASSERT_ALWAYS(IsInteger(M, m));
    CoCoA_ASSERT_ALWAYS(n == N);
    CoCoA_ASSERT_ALWAYS(m == M);

    CoCoA_ASSERT_ALWAYS(!IsInvertible(n) || IsOne(n*(1/n)));

    // Some simple arithmetic
    n = n+m; N=N+M; CoCoA_ASSERT_ALWAYS(n == N);
    n = n-m; N=N-M; CoCoA_ASSERT_ALWAYS(n == N);
    n = n*m; N=N*M; CoCoA_ASSERT_ALWAYS(n == N);
    n = n/m; N=N/M; CoCoA_ASSERT_ALWAYS(n == N);

    n = n+M; N=N+M; CoCoA_ASSERT_ALWAYS(n == N);
    n = n-M; N=N-M; CoCoA_ASSERT_ALWAYS(n == N);
    n = n*M; N=N*M; CoCoA_ASSERT_ALWAYS(n == N);
    n = n/M; N=N/M; CoCoA_ASSERT_ALWAYS(n == N);

    n = N+m; N=N+M; CoCoA_ASSERT_ALWAYS(n == N);
    n = N-m; N=N-M; CoCoA_ASSERT_ALWAYS(n == N);
    n = N*m; N=N*M; CoCoA_ASSERT_ALWAYS(n == N);
    n = N/m; N=N/M; CoCoA_ASSERT_ALWAYS(n == N);

    // Checks for aliasing problems
    n = n+n; N=N+N; CoCoA_ASSERT_ALWAYS(n == N);
    n = n*n; N=N*N; CoCoA_ASSERT_ALWAYS(n == N);
    n = n/n; N=N/N; CoCoA_ASSERT_ALWAYS(n == N);
    n = n-n; N=N-N; CoCoA_ASSERT_ALWAYS(n == N);

    // assignment arithmetic
    N = 123456789;  n = N; CoCoA_ASSERT_ALWAYS(n == N);
    n += m; N+=M; CoCoA_ASSERT_ALWAYS(n == N);
    n -= m; N-=M; CoCoA_ASSERT_ALWAYS(n == N);
    n *= m; N*=M; CoCoA_ASSERT_ALWAYS(n == N);
    n /= m; N/=M; CoCoA_ASSERT_ALWAYS(n == N);

    // Checks for aliasing problems (NB clang++ issues warnings)
    /*TEST SELF ASSIGN*/ n = m; n += n; N=M; N+=M; CoCoA_ASSERT_ALWAYS(n == N);
    /*TEST SELF ASSIGN*/ n = m; n -= n; N=M; N-=M; CoCoA_ASSERT_ALWAYS(n == N);
    /*TEST SELF ASSIGN*/ n = m; n *= n; N=M; N*=M; CoCoA_ASSERT_ALWAYS(n == N);
    /*TEST SELF ASSIGN*/ n = m; n /= n; N=M; N/=M; CoCoA_ASSERT_ALWAYS(n == N);
  }


  void LongTest(ring RR)
  {
    BasicChecks(RR);

    tautology(zero(RR));
    tautology(RingElem(RR, 1));
    tautology(RingElem(RR, -1));
    tautology(RingElem(RR, 32003));
    tautology(RingElem(RR, -32003));
    if (PrecisionBits(RR) > 200)
    {
      tautology(RingElem(RR, power(17,99)));
      tautology(RingElem(RR, power(-17,99)));
    }

    // Check correct arithmetic.  Cases are small/large, positive/negative.
    arithmetic(RingElem(RR, 23), RingElem(RR, 37));
    arithmetic(RingElem(RR, 23), -RingElem(RR, 37));
    arithmetic(-RingElem(RR, 23), RingElem(RR, 37));
    arithmetic(-RingElem(RR, 23), -RingElem(RR, 37));

    if (PrecisionBits(RR) > 326)
    {
      arithmetic(RingElem(RR, 100), RingElem(RR, power(10, 100)));
      arithmetic(RingElem(RR, 100), -RingElem(RR, power(10, 100)));
      arithmetic(-RingElem(RR, 100), RingElem(RR, power(10, 100)));
      arithmetic(-RingElem(RR, 100), -RingElem(RR, power(10, 100)));

      arithmetic(RingElem(RR, power(10, 100)), RingElem(RR, 100));
      arithmetic(RingElem(RR, power(10, 100)), -RingElem(RR, 100));
      arithmetic(-RingElem(RR, power(10, 100)), RingElem(RR, 100));
      arithmetic(-RingElem(RR, power(10, 100)), -RingElem(RR, 100));

      arithmetic(RingElem(RR, power(10, 100)), RingElem(RR, power(2, 301)));
      arithmetic(RingElem(RR, power(10, 100)), -RingElem(RR, power(2, 301)));
      arithmetic(-RingElem(RR, power(10, 100)), RingElem(RR, power(2, 301)));
      arithmetic(-RingElem(RR, power(10, 100)), -RingElem(RR, power(2, 301)));
    }

    // compute a large number
    const RingElem big = power(RingElem(RR,10), 10000);
    CoCoA_ASSERT_ALWAYS(big == power(10,10000));
    CoCoA_ASSERT_ALWAYS(power(RingElem(RR,-3),3) == RingElem(RR,-3)*RingElem(RR,-3)*RingElem(RR,-3));
    const RingElem third = RingElem(RR,BigRat(1,3));
    CoCoA_ASSERT_ALWAYS(IsZero(3*third-1));
    BigInt TwentySeven;
    CoCoA_ASSERT_ALWAYS(IsInteger(TwentySeven, power(third, -3)) && TwentySeven == 27);
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    // Try at very low precision, but offering a wider "buffer" than usual.
    LongTest(NewRingTwinFloat(8,16,32));

    // Try at several small to moderate precisions.
    for (int prec=20; prec <= 1000; prec += 10)
    {
      LongTest(NewRingTwinFloat(prec));
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
