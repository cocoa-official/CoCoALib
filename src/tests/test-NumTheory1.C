//   Copyright (c)  2009-2010  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/NumTheory-gcd.H"
#include "CoCoA/NumTheory-prime.H"
#include "CoCoA/NumTheory-factor.H"
#include "CoCoA/NumTheory-modular.H"
#include "CoCoA/NumTheory-misc.H"
#include "CoCoA/VectorOps.H"
#include "CoCoA/error.H"

#include <cmath>
using std::log;
#include <cstdlib>
using std::abs;
#include <iomanip>
using std::flush;
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <limits>
using std::numeric_limits;
#include <vector>
using std::vector;

namespace CoCoA
{

  void program()
  {
    using std::abs;
    GlobalManager CoCoAFoundations(UseGMPAllocator);

    cout << "Testing gcd and lcm on small integers." << endl;
    constexpr int Nmax = 257;
    for (int a=-Nmax; a < Nmax; ++a)
      for (int b=-Nmax; b < Nmax; ++b)
      {
        const BigInt A(a);
        const BigInt B(b);
        const int gcdab = gcd(a,b);
        CoCoA_ASSERT_ALWAYS(gcdab >= 0);
        CoCoA_ASSERT_ALWAYS((a == 0 && b == 0) || gcdab > 0);
        CoCoA_ASSERT_ALWAYS(gcdab == 0 || (a%gcdab == 0 && b%gcdab == 0));
        CoCoA_ASSERT_ALWAYS(gcdab == gcd(a,B));
        CoCoA_ASSERT_ALWAYS(gcdab == gcd(A,b));
        CoCoA_ASSERT_ALWAYS(gcdab == gcd(A,B));
        if (gcdab != 0)
        {
          CoCoA_ASSERT_ALWAYS(IsCoprime(a/gcdab, b/gcdab));
          CoCoA_ASSERT_ALWAYS(IsCoprime(A/gcdab, b/gcdab));
          CoCoA_ASSERT_ALWAYS(IsCoprime(a/gcdab, B/gcdab));
          CoCoA_ASSERT_ALWAYS(IsCoprime(A/gcdab, B/gcdab));
        }
        
        const int lcmab = lcm(a,b);
        CoCoA_ASSERT_ALWAYS(gcdab*lcmab == abs(a*b));
        CoCoA_ASSERT_ALWAYS(lcmab == lcm(a,B));
        CoCoA_ASSERT_ALWAYS(lcmab == lcm(A,b));
        CoCoA_ASSERT_ALWAYS(lcmab == lcm(A,B));
      }

    cout << "\nTesting ExtGcd on small integers" << endl;
    for (int a=-Nmax; a < Nmax; ++a)
      for (int b=-Nmax; b < Nmax; ++b)
      {
        if (a == 0 && b == 0) continue;
        long CofacA, CofacB;
        const int g = ExtGcd(CofacA,CofacB,a,b);
        const int g2 = gcd(a,b);
        CoCoA_ASSERT_ALWAYS(g == g2);
        CoCoA_ASSERT_ALWAYS(CofacA*a+CofacB*b == g);
        if (abs(a) == abs(b))
        {
          CoCoA_ASSERT_ALWAYS((CofacA == 0 && CofacB == sign(b)) || (CofacB == 0 && CofacA == sign(a)));
          continue;
        }
        if (abs(a) == g) { CoCoA_ASSERT_ALWAYS(CofacA == sign(a) && CofacB == 0); continue; }
        if (abs(b) == g) { CoCoA_ASSERT_ALWAYS(CofacB == sign(b) && CofacA == 0); continue; }
        // General case
        CoCoA_ASSERT_ALWAYS(abs(CofacA) <= (abs(b)/g)/2  && abs(CofacB) <= (abs(a)/g)/2);
        if (2*g == abs(b))
          CoCoA_ASSERT_ALWAYS(abs(CofacA) == 1);
        else
          CoCoA_ASSERT_ALWAYS(2*abs(CofacA) < abs(b)/g);
        if (2*g == abs(a))
          CoCoA_ASSERT_ALWAYS(abs(CofacB) == 1);
        else
          CoCoA_ASSERT_ALWAYS(2*abs(CofacB) < abs(a)/g);
      }


    cout << "\nTesting ExtGcd on BigInt integers (but with small values)" << endl;
    for (int a=-Nmax; a < Nmax; ++a)
      for (int b=-Nmax; b < Nmax; ++b)
      {
        if (a == 0 && b == 0) continue;
        BigInt A(a);
        BigInt B(b);
        BigInt CofacA, CofacB;
        const BigInt G = ExtGcd(CofacA,CofacB,A,B);
        const int g = gcd(a,b);
        CoCoA_ASSERT_ALWAYS(G == g);
        CoCoA_ASSERT_ALWAYS(CofacA*a+CofacB*b == g);
        if (abs(a) == abs(b))
        {
          CoCoA_ASSERT_ALWAYS((CofacA == 0 && CofacB == sign(b)) || (CofacB == 0 && CofacA == sign(a)));
          continue;
        }
        if (abs(a) == g) { CoCoA_ASSERT_ALWAYS(CofacA == sign(a) && CofacB == 0); continue; }
        if (abs(b) == g) { CoCoA_ASSERT_ALWAYS(CofacB == sign(b) && CofacA == 0); continue; }
        // General case
        CoCoA_ASSERT_ALWAYS(abs(CofacA) <= (abs(b)/g)/2  && abs(CofacB) <= (abs(a)/g)/2);
        if (2*g == abs(b))
          CoCoA_ASSERT_ALWAYS(abs(CofacA) == 1);
        else
          CoCoA_ASSERT_ALWAYS(2*abs(CofacA) < abs(b)/g);
        if (2*g == abs(a))
          CoCoA_ASSERT_ALWAYS(abs(CofacB) == 1);
        else
          CoCoA_ASSERT_ALWAYS(2*abs(CofacB) < abs(a)/g);
      }


    cout << "\nTesting InvMod" << endl;
    for (int r=-Nmax; r < Nmax; ++r)
      for (int m=2; m < Nmax; ++m)
      {
        if (gcd(r,m) != 1) continue;  // skip cases where inverse does not exist
        BigInt R(r);
        BigInt M(m);
        const long invr = InvMod(r,m);
        if (r > 0) CoCoA_ASSERT_ALWAYS((r*invr)%m == 1);
        if (r < 0) CoCoA_ASSERT_ALWAYS(((-r)*invr)%m == m-1);
        CoCoA_ASSERT_ALWAYS(invr == InvMod(R,m));
        CoCoA_ASSERT_ALWAYS(invr == InvMod(r,M));
        CoCoA_ASSERT_ALWAYS(invr == InvMod(R,M));
      }

    cout << "\nTesting InvMod when inverse does not exist" << endl;
    for (int r=-Nmax; r < Nmax; ++r)
      for (int m=2; m < Nmax; ++m)
      {
        if (gcd(r,m) == 1) continue;  // skip cases where inverse exists
        BigInt R(r);
        BigInt M(m);
        try { InvMod(r,m); CoCoA_ASSERT_ALWAYS(0 && "SHOULD NEVER GET HERE"); }
        /* Ignore expected DivByZero, o/w rethrow: */
        catch (const CoCoA::ErrorInfo& err) { if (err != ERR::DivByZero) throw; }

        try { InvMod(r,M); CoCoA_ASSERT_ALWAYS(0 && "SHOULD NEVER GET HERE"); }
        /* Ignore expected DivByZero, o/w rethrow: */
        catch (const CoCoA::ErrorInfo& err) { if (err != ERR::DivByZero) throw; }

        try { InvMod(R,m); CoCoA_ASSERT_ALWAYS(0 && "SHOULD NEVER GET HERE"); }
        /* Ignore expected DivByZero, o/w rethrow: */
        catch (const CoCoA::ErrorInfo& err) { if (err != ERR::DivByZero) throw; }

        try { InvMod(R,M); CoCoA_ASSERT_ALWAYS(0 && "SHOULD NEVER GET HERE"); }
        /* Ignore expected DivByZero, o/w rethrow: */
        catch (const CoCoA::ErrorInfo& err) { if (err != ERR::DivByZero) throw; }

      }

    cout << "\nTesting PowerMod" << endl;
    for (int r=-Nmax; r < Nmax; ++r)
      for (int m=2; m < Nmax; ++m)
      {
        BigInt R(r);
        BigInt M(m);
        CoCoA_ASSERT_ALWAYS(PowerMod(r,0,m) == 1);
        CoCoA_ASSERT_ALWAYS(PowerMod(r,0,M) == 1);
        CoCoA_ASSERT_ALWAYS(PowerMod(R,0,m) == 1);
        CoCoA_ASSERT_ALWAYS(PowerMod(R,0,M) == 1);

        long RmodM;
        if (r >= 0) RmodM = r%m; else RmodM = (m-((-r)%m))%m;
        CoCoA_ASSERT_ALWAYS(PowerMod(r,1,m) == RmodM);
        CoCoA_ASSERT_ALWAYS(PowerMod(r,1,M) == RmodM);
        CoCoA_ASSERT_ALWAYS(PowerMod(R,1,m) == RmodM);
        CoCoA_ASSERT_ALWAYS(PowerMod(R,1,M) == RmodM);

        if (gcd(r,m) == 1)
        {
          long InvR = InvMod(r,m);
          CoCoA_ASSERT_ALWAYS(PowerMod(r,-1,m) == InvR);
          CoCoA_ASSERT_ALWAYS(PowerMod(r,-1,M) == InvR);
          CoCoA_ASSERT_ALWAYS(PowerMod(R,-1,m) == InvR);
          CoCoA_ASSERT_ALWAYS(PowerMod(R,-1,M) == InvR);
        }
        for (int e=0; e < 9; ++e)
        {
          CoCoA_ASSERT_ALWAYS(PowerMod(PowerMod(r,e,m), 2, m) == PowerMod(r, 2*e, m));
          if (gcd(r,m) == 1)
            CoCoA_ASSERT_ALWAYS((PowerMod(r,e,m)*PowerMod(r,-e,m))%m == 1);
        }
      }



    cout << "\nTesting IsPrime for small positive integers." << endl;
    int PrimeCounter = 0;
    for (int n=1; n <= 1000; ++n)
    {
      if (IsPrime(n)) ++PrimeCounter;
      CoCoA_ASSERT_ALWAYS(IsPrime(n) == IsPrime(BigInt(n)));
    }
    CoCoA_ASSERT_ALWAYS(PrimeCounter == 168);

    // ATTN!  this next test depends on the platform (32-bit or 64-bit)
    const int LongBits = numeric_limits<long>::digits;
    cout << "\nTesting IsPrime for numbers just greater than max long." << endl;
    const BigInt StartingPoint = power(2,LongBits);
    for (int delta=1; delta < 99; delta+=2)
    {
      const BigInt N = StartingPoint+delta;
      if (!IsProbPrime(N)) continue;
      if (IsPrime(N)) continue;
      cout << "Surprise: " << N << " is a composite pseudo-prime." << endl;
    }


    cout << "\nTesting IsProbPrime: lengths of repunit (probable) primes (up to 500):" << flush;
    for (int n=1; n <= 500; ++n)
    {
      if (!IsPrime(n)) continue; // no need to try composite n
      if (IsProbPrime((power(10,n)-1)/9))
        cout << " " << n << flush;
    }
    cout << endl;


    cout << "\nTesting NextPrime & PrevPrime" << endl;
    const long MaxLong = numeric_limits<long>::max();
    CoCoA_ASSERT_ALWAYS(NextPrime(0) == 2);
    CoCoA_ASSERT_ALWAYS(NextPrime(1) == 2);
    CoCoA_ASSERT_ALWAYS(NextPrime(2) == 3);
    CoCoA_ASSERT_ALWAYS(NextPrime(3) == 5);
    CoCoA_ASSERT_ALWAYS(NextPrime(4) == 5);

    CoCoA_ASSERT_ALWAYS(NextPrime(MaxLong) == 0);
    for (long n=MaxLong; !IsPrime(n--); )
      CoCoA_ASSERT_ALWAYS(NextPrime(n) == 0);

    const long MaxPrime = IsPrime(MaxLong)?MaxLong:PrevPrime(MaxLong);
    CoCoA_ASSERT_ALWAYS(IsPrime(MaxPrime));
    CoCoA_ASSERT_ALWAYS(NextPrime(MaxPrime) == 0);
    for (long n=MaxPrime; !IsPrime(n--); )
      CoCoA_ASSERT_ALWAYS(NextPrime(n) == MaxPrime);

    // Check that PrevPrime throws OutOfRange when given 0,1,2
    for (int n=0; n <= 2; ++n)
    {
      try
      {
        const long NotDefd = PrevPrime(n); // should throw error
        /*UNUSED VALUE*/ (void)(sizeof(NotDefd)); // just to keep compiler quiet
        CoCoA_ASSERT_ALWAYS(false && "PrevPrime should have thrown");
      }
      catch (const ErrorInfo& err)
      {
        if (err != ERR::OutOfRange) throw;
      }
    }
    CoCoA_ASSERT_ALWAYS(PrevPrime(4) == 3);
    CoCoA_ASSERT_ALWAYS(PrevPrime(3) == 2);
    for (int n=3; n < 1000; n += 2)
    {
      if (IsPrime(n))
      {
        CoCoA_ASSERT_ALWAYS(PrevPrime(NextPrime(n)) == n);
        CoCoA_ASSERT_ALWAYS(NextPrime(PrevPrime(n)) == n);
      }
      else
      {
        CoCoA_ASSERT_ALWAYS(NextPrime(PrevPrime(n)) == NextPrime(n));
        CoCoA_ASSERT_ALWAYS(PrevPrime(NextPrime(n)) == PrevPrime(n));
      }
    }

    //-------------------------------------------------------
    cout << "\nTesting RandomNBitPrime" << endl;
    for (long n=10; n < 32; ++n)
    {
      for (int trial=0; trial < 10; ++trial)
      {
        const long p = RandomNBitPrime(n);
        CoCoA_ASSERT_ALWAYS(p > (1L<<(n-1))-30 && (p-30)/2 < (1L<<(n-1))); // range check -- allow out of range by up to 30!
        CoCoA_ASSERT_ALWAYS(IsPrime(p));
      }
    }

    //-------------------------------------------------------
    cout << "\nTesting PrimitiveRoot" << endl;
    constexpr int MaxHist = 20;
    vector<int> histogram(MaxHist);
    int p = 1;
    while (p < 65535)
    {
      p = NextPrime(p);
      const int g = PrimitiveRoot(p);
      if (g < MaxHist) ++histogram[g];
//     double ratio = g/(std::log(p)*std::log(std::log(p)));
//     if (ratio > 2.0) cout << "ratio=" << ratio << "  for p=" << p << " and g=" << g << endl;
    }
    cout << "Histogram of least positive primitive roots: " << histogram << endl;


    cout << "\nTesting factorize on smaller & larger integers:" << endl;
    for (int i=1; i <= 16; ++i)
    {
      const BigInt N = power(10,i)-1;
      cout << "Factorization of " << N << ": " << factor(N) << endl;
    }

    // --------------------------------------------
    cout << "\nTesting EulerTotient on small numbers." << endl;
    cout << "Here are the most frequent values of EulerTotient on numbers up to 100000:\n"
         << "Value  Freq\n";
    vector<int> freq(100000);
    for (int i=2; i < 100000; ++i)
      ++freq[EulerTotient(i)];
    for (int i=2; i < 100000; ++i)
      if (freq[i] > 200) 
        cout << i << "  " << freq[i] << endl;

    cout << "\nTesting EulerTotient on larger numbers." << endl
         << "Looking for smallest number with EulerTotient(N) < N/10." << endl;
    BigInt BigNum;
    BigNum = 2;
    int LastPrime = -1; // dummy initial value -- never used
    for (int p=3; 10*EulerTotient(BigNum) > BigNum; p=NextPrime(p))
    {
      BigNum *= p;
      LastPrime = p;
    }
    BigInt PhiBigNum = EulerTotient(BigNum);
    cout << "Result is  N=" << BigNum << endl
         << "Its EulerTotient=" << PhiBigNum << endl
         << "In fact N is just the product of all primes up to " << LastPrime << endl;


    // --------------------------------------------
    // factor_TrialDiv
    cout << "\nTesting factor_TrialDiv" << endl;
    cout << "Small factorization: " << factor_TrialDiv(21467*99991L, 100000) << endl;
    factorization<BigInt> FactoredFactorial = factor_TrialDiv(factorial(101), 101);
    cout << "101-smooth factors of factorial of 101: " << FactoredFactorial << endl;
    FactoredFactorial = factor_TrialDiv(factorial(101), 50);
    cout << "50-smooth factors of factorial of 101: " << FactoredFactorial << endl;
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
