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
#include "CoCoA/error.H"
#include "CoCoA/BigIntOps.H"
#include "CoCoA/NumTheory-CRT.H"
#include "CoCoA/NumTheory-prime.H"
#include "CoCoA/NumTheory-CoprimeFactorBasis.H"
#include "CoCoA/utils.H"

#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;


// Tests eratosthenes, CoprimeFactorBase

namespace CoCoA
{

  // see redmine 1300
  void test_CFB()
  {
    vector<long> v;
    v.push_back(2), v.push_back(-3);
    CoprimeFactorBasis_BigInt CFB;
    CFB.myAddInfo(v);
    CoCoA_ASSERT_ALWAYS(FactorBase(CFB)[1] == 3);
  }

  long NumPrimes(const vector<bool>& sieve)
  {
    long count = 1; // start from 1 because 2 is not in the sieve
    const long n = len(sieve);
    for (long i=0; i < n; ++i)
      count += sieve[i];
    return count;
  }

  void test_eratosthenes()
  {
    vector<bool> sieve;
    sieve = eratosthenes(3);
    CoCoA_ASSERT_ALWAYS(len(sieve) == 2);
    CoCoA_ASSERT_ALWAYS(!sieve[0]);
    CoCoA_ASSERT_ALWAYS(sieve[1]);

    sieve = eratosthenes(100);
    CoCoA_ASSERT_ALWAYS(NumPrimes(sieve) == 26);

    sieve = eratosthenes(123456);
    for (int i=1; i < 123456; i += 2)
      CoCoA_ASSERT_ALWAYS(IsPrime(i) == sieve[i/2]);
  }


  // Test for CRTMill (very simple test)
  void test_CRT()
  {
    // Copied from ex-NumTheory2.C
    const BigInt N = power(10,100);
    const BigInt UPB = 2*N+1;

    CRTMill crt;
    int p = 101;
    while (true)
    {
      p = NextPrime(p);
      crt.myAddInfo(N%p, p); // tell crt the new residue-modulus pair
      if (CombinedModulus(crt) > UPB) break;
    }

    // Check answer is correct.
    if (CombinedResidue(crt) != N)
      CoCoA_THROW_ERROR("Wrong answer", "CoCoA::Program");
  }
  

  void program()
  {
    GlobalManager CoCoAFoundations;

    test_eratosthenes();
    test_CFB();
    test_CRT();
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
