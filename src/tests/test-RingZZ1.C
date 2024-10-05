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
    const ring ZZ = RingZZ();

    CoCoA_ASSERT_ALWAYS(IsZZ(ZZ));

    const BigInt Zero;
    RingElem n(ZZ);
    // Check an implementation detail
    RingElem m = n;
    CoCoA_ASSERT_ALWAYS(m == n);
    CoCoA_ASSERT_ALWAYS(raw(m).myRawPtr() != raw(n).myRawPtr());
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

    n = power(RingElem(ZZ,10), 1000);
    CoCoA_ASSERT_ALWAYS(!IsInvertible(n));
    CoCoA_ASSERT_ALWAYS(n == power(BigInt(10), 1000));

    tautology(zero(ZZ));
    tautology(RingElem(ZZ, 1));
    tautology(RingElem(ZZ, -1));
    tautology(RingElem(ZZ, 32003));
    tautology(RingElem(ZZ, -32003));
    tautology(RingElem(ZZ, power(BigInt(17),99)));
    tautology(RingElem(ZZ, power(BigInt(-17),99)));


    // Check correct arithmetic.  Cases are small/large, positive/negative.
    arithmetic(RingElem(ZZ, 23), RingElem(ZZ, 37));
    arithmetic(RingElem(ZZ, 100), RingElem(ZZ, power(BigInt(10), 100)));
    arithmetic(RingElem(ZZ, power(BigInt(10), 100)), RingElem(ZZ, 100));
    arithmetic(RingElem(ZZ, power(BigInt(10), 100)), RingElem(ZZ, power(BigInt(2), 301)));

    arithmetic(RingElem(ZZ, 23), -RingElem(ZZ, 37));
    arithmetic(RingElem(ZZ, 100), -RingElem(ZZ, power(BigInt(10), 100)));
    arithmetic(RingElem(ZZ, power(BigInt(10), 100)), -RingElem(ZZ, 100));
    arithmetic(RingElem(ZZ, power(BigInt(10), 100)), -RingElem(ZZ, power(BigInt(2), 301)));

    arithmetic(-RingElem(ZZ, 23), RingElem(ZZ, 37));
    arithmetic(-RingElem(ZZ, 100), RingElem(ZZ, power(BigInt(10), 100)));
    arithmetic(-RingElem(ZZ, power(BigInt(10), 100)), RingElem(ZZ, 100));
    arithmetic(-RingElem(ZZ, power(BigInt(10), 100)), RingElem(ZZ, power(BigInt(2), 301)));

    arithmetic(-RingElem(ZZ, 23), -RingElem(ZZ, 37));
    arithmetic(-RingElem(ZZ, 100), -RingElem(ZZ, power(BigInt(10), 100)));
    arithmetic(-RingElem(ZZ, power(BigInt(10), 100)), -RingElem(ZZ, 100));
    arithmetic(-RingElem(ZZ, power(BigInt(10), 100)), -RingElem(ZZ, power(BigInt(2), 301)));

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
