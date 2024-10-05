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
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/error.H"
#include "CoCoA/ring.H"

#include <iostream>
using std::cerr;
using std::endl;



namespace CoCoA
{

  void tautology(const RingElem& n)
  {
    CoCoA_ASSERT_ALWAYS(IsZero(n) ^ IsInvertible(n));
    BigInt N;
    CoCoA_ASSERT_ALWAYS(IsInteger(N, n));
    CoCoA_ASSERT_ALWAYS(n == N);
    CoCoA_ASSERT_ALWAYS(N == n);
    RingElem n2(owner(n), N);
    CoCoA_ASSERT_ALWAYS(n == n2);
    CoCoA_ASSERT_ALWAYS(n2 == n);
  }


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


  void arithmetic(RingElem n, const RingElem& m)
  {
    BigInt N;
    CoCoA_ASSERT_ALWAYS(IsInteger(N, n));
    BigInt M;
    CoCoA_ASSERT_ALWAYS(IsInteger(M, m));
    CoCoA_ASSERT_ALWAYS(n == N);
    CoCoA_ASSERT_ALWAYS(m == M);

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

  void program()
  {
    GlobalManager CoCoAFoundations;

    const ring QQ = RingQQ();

    CoCoA_ASSERT_ALWAYS(IsFractionField(QQ));
    CoCoA_ASSERT_ALWAYS(IsZZ(BaseRing(QQ)));
    CoCoA_ASSERT_ALWAYS(BaseRing(QQ) == RingZZ());

    const BigInt Zero;
    RingElem n(QQ);
    // Check an implementation detail
    RingElem m = n; CoCoA_ASSERT_ALWAYS(m == n); CoCoA_ASSERT_ALWAYS(raw(m).myRawPtr() != raw(n).myRawPtr());

    CoCoA_ASSERT_ALWAYS(IsZero(n));
    CoCoA_ASSERT_ALWAYS(!IsInvertible(n));
    CoCoA_ASSERT_ALWAYS(!IsOne(n));
    CoCoA_ASSERT_ALWAYS(!IsMinusOne(n));
    equal(n, 0);
    equal(n, Zero);
    greater(n, -1);
    greater(n, BigInt(-1));
    less(n, 1);
    less(n, BigInt(1));

    n = 1;
    CoCoA_ASSERT_ALWAYS(!IsZero(n));
    CoCoA_ASSERT_ALWAYS(IsInvertible(n));
    CoCoA_ASSERT_ALWAYS(IsOne(n));
    CoCoA_ASSERT_ALWAYS(!IsMinusOne(n));
    equal(n, 1);
    equal(n, BigInt(1));
    greater(n, -1);
    greater(n, BigInt(-1));
    less(n, 2);
    less(n, BigInt(2));

    n = -1;
    CoCoA_ASSERT_ALWAYS(!IsZero(n));
    CoCoA_ASSERT_ALWAYS(IsInvertible(n));
    CoCoA_ASSERT_ALWAYS(!IsOne(n));
    CoCoA_ASSERT_ALWAYS(IsMinusOne(n));
    equal(n, -1);
    equal(n, BigInt(-1));
    greater(n, -2);
    greater(n, BigInt(-2));
    less(n, 1);
    less(n, BigInt(1));

    n = power(RingElem(QQ,10), 1000);
    CoCoA_ASSERT_ALWAYS(IsInvertible(n));
    CoCoA_ASSERT_ALWAYS(n == power(BigInt(10), 1000));

    tautology(zero(QQ));
    tautology(RingElem(QQ, 1));
    tautology(RingElem(QQ, -1));
    tautology(RingElem(QQ, 32003));
    tautology(RingElem(QQ, -32003));
    tautology(RingElem(QQ, power(BigInt(17),99)));
    tautology(RingElem(QQ, power(BigInt(-17),99)));

    // Check correct arithmetic.  Cases are small/large, positive/negative.
    arithmetic(RingElem(QQ, 23), RingElem(QQ, 37));
    arithmetic(RingElem(QQ, 100), RingElem(QQ, power(BigInt(10), 100)));
    arithmetic(RingElem(QQ, power(BigInt(10), 100)), RingElem(QQ, 100));
    arithmetic(RingElem(QQ, power(BigInt(10), 100)), RingElem(QQ, power(BigInt(2), 301)));

    arithmetic(RingElem(QQ, 23), -RingElem(QQ, 37));
    arithmetic(RingElem(QQ, 100), -RingElem(QQ, power(BigInt(10), 100)));
    arithmetic(RingElem(QQ, power(BigInt(10), 100)), -RingElem(QQ, 100));
    arithmetic(RingElem(QQ, power(BigInt(10), 100)), -RingElem(QQ, power(BigInt(2), 301)));

    arithmetic(-RingElem(QQ, 23), RingElem(QQ, 37));
    arithmetic(-RingElem(QQ, 100), RingElem(QQ, power(BigInt(10), 100)));
    arithmetic(-RingElem(QQ, power(BigInt(10), 100)), RingElem(QQ, 100));
    arithmetic(-RingElem(QQ, power(BigInt(10), 100)), RingElem(QQ, power(BigInt(2), 301)));

    arithmetic(-RingElem(QQ, 23), -RingElem(QQ, 37));
    arithmetic(-RingElem(QQ, 100), -RingElem(QQ, power(BigInt(10), 100)));
    arithmetic(-RingElem(QQ, power(BigInt(10), 100)), -RingElem(QQ, 100));
    arithmetic(-RingElem(QQ, power(BigInt(10), 100)), -RingElem(QQ, power(BigInt(2), 301)));

    // assignment arithmetic
    n = 355; n /= 113; // Aside: 355/113 is a famous approximation for pi.
    m = n;

    n += m;
    n -= m;
    n *= m;
    n /= m;

    // Testing for aliasing problems
    /*TEST SELF ASSIGN*/ n = m; n += n; CoCoA_ASSERT_ALWAYS(n == 2*m);
    /*TEST SELF ASSIGN*/ n = m; n -= n; CoCoA_ASSERT_ALWAYS(IsZero(n));
    /*TEST SELF ASSIGN*/ n = m; n *= n; CoCoA_ASSERT_ALWAYS(n == m*m);
    /*TEST SELF ASSIGN*/ n = m; n /= n; CoCoA_ASSERT_ALWAYS(IsOne(n));

    // compute a large number
    n = power(m, 100);
    CoCoA_ASSERT_ALWAYS(IsInvertible(n));
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
