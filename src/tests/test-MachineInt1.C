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
#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <limits>
using std::numeric_limits;



namespace CoCoA
{

  void CheckRoundDiv(const MachineInt& a, const MachineInt& b)
  {
    if (IsZero(b)) return; // skip the check if b == 0
    const BigInt A(a);
    const BigInt B(b);

    const BigInt Q = RoundDiv(A,B);
    long q;
    if (IsConvertible(q, Q))
    {
      CoCoA_ASSERT_ALWAYS(q == RoundDiv(a,b));
      return;
    }
    // Here we know that RoundDiv must fail (because answer is too big),
    // so we verify that it does indeed fail (with ERR::BadCheckedCast).
    try
    {
      q = RoundDiv(a,b); // this must overflow
      CoCoA_ASSERT_ALWAYS(!"NEVER GET HERE!");
    }
    catch (const ErrorInfo& err) { CoCoA_ASSERT_ALWAYS (err == ERR::BadConvert); }
  }


  void program()
  {
    // This test checks that the rounded division function for pairs of machine
    // integers is compatible with that for BigInts.  It tries a number of limit
    // cases where the machine integer version has to be careful about overflow.
    GlobalManager CoCoAFoundations;

    const unsigned long ulmax = numeric_limits<unsigned long>::max();
    const long lmax = numeric_limits<long>::max();
    const long lmin = numeric_limits<long>::min();
    for (int i=0; i < 100; ++i)
    {
      CheckRoundDiv(ulmax-i, 1);
      CheckRoundDiv(ulmax-i, -1);

      CheckRoundDiv(i, 1);
      CheckRoundDiv(i, -1);

      CheckRoundDiv(lmax-i, 1);
      CheckRoundDiv(lmax-i, -1);

      CheckRoundDiv(lmin+i, 1);
      CheckRoundDiv(lmin+i, -1);

      CheckRoundDiv(1, 1);
      CheckRoundDiv(1, -1);

      CheckRoundDiv(-1, 1);
      CheckRoundDiv(-1, -1);

      for (int j=0; j < 100; ++j)
      {
        CheckRoundDiv(i,j);
        CheckRoundDiv(-i,j);
        CheckRoundDiv(i,-j);
        CheckRoundDiv(-i,-j);

        CheckRoundDiv(i, ulmax-j);
        CheckRoundDiv(i, lmax-j);
        CheckRoundDiv(i, lmin+j);

        CheckRoundDiv(-i, ulmax-j);
        CheckRoundDiv(-i, lmax-j);
        CheckRoundDiv(-i, lmin+j);

        CheckRoundDiv(ulmax-i, j);
        CheckRoundDiv(ulmax-i, -j);
        CheckRoundDiv(ulmax-i, ulmax-j);
        CheckRoundDiv(ulmax-i, lmax-j);
        CheckRoundDiv(ulmax-i, lmin+j);

        CheckRoundDiv(lmax-i, j);
        CheckRoundDiv(lmax-i, -j);
        CheckRoundDiv(lmax-i, ulmax-j);
        CheckRoundDiv(lmax-i, lmax-j);
        CheckRoundDiv(lmax-i, lmin+j);

        CheckRoundDiv(lmin+i, j);
        CheckRoundDiv(lmin+i, -j);
        CheckRoundDiv(lmin+i, ulmax-j);
        CheckRoundDiv(lmin+i, lmax-j);
        CheckRoundDiv(lmin+i, lmin+j);

        CheckRoundDiv(1, j);
        CheckRoundDiv(1, -j);
        CheckRoundDiv(1, ulmax-j);
        CheckRoundDiv(1, lmax-j);
        CheckRoundDiv(1, lmin+j);

        CheckRoundDiv(-1, j);
        CheckRoundDiv(-1, -j);
        CheckRoundDiv(-1, ulmax-j);
        CheckRoundDiv(-1, lmax-j);
        CheckRoundDiv(-1, lmin+j);
      }
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
