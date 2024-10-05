//   Copyright (c)  2010  John Abbott,  Anna M. Bigatti

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

#include <cstdlib>
using std::abs;
#include <iostream>
using std::clog;
using std::cerr;
using std::endl;
#include<limits>
using std::numeric_limits;


namespace CoCoA
{

  // test IsExactRoot
  void program()
  {
    using std::abs;
    // This test checks the functions FloorRoot, FloorSqrt, IsSquare & IsPower.
    GlobalManager CoCoAFoundations(UseGMPAllocator);

    BigInt ROOT;   // used as arg to IsExactFloorRoot
    long root;     // used as arg to IsExactFloorRoot
    BigInt N;
    for (int a=0; a <= 99; ++a)
    {
      CoCoA_ASSERT_ALWAYS(FloorRoot(a,1) == a);
      CoCoA_ASSERT_ALWAYS(FloorRoot(BigInt(a),1) == a);
      
      for (int b=2; b <= 99; ++b)
      {
        N = power(a, b);
        CoCoA_ASSERT_ALWAYS(IsExactFloorRoot(ROOT,N,b));
        {
          long n;
          if (IsConvertible(n, N))
          {
            CoCoA_ASSERT_ALWAYS(IsExactFloorRoot(root,n,b));
            CoCoA_ASSERT_ALWAYS(root == ROOT);
          }
        }

        CoCoA_ASSERT_ALWAYS(a == FloorRoot(N,b));
        CoCoA_ASSERT_ALWAYS(a == ROOT);

        CoCoA_ASSERT_ALWAYS(IsPower(N));
        if (!((N == 0) || (N == 8)))
        {
          CoCoA_ASSERT_ALWAYS(!IsPower(N+1));
          CoCoA_ASSERT_ALWAYS(!IsExactFloorRoot(ROOT, N+1,b));
          CoCoA_ASSERT_ALWAYS(ROOT == a);
        }
        if (!((N == 0) || (N == 1) || (N == 9)))
        {
          CoCoA_ASSERT_ALWAYS(!IsPower(N-1));
          CoCoA_ASSERT_ALWAYS(!IsExactFloorRoot(ROOT, N-1,b));
          CoCoA_ASSERT_ALWAYS(ROOT == a-1);
        }
        if (IsEven(b))
        {
          CoCoA_ASSERT_ALWAYS(IsSquare(N));
          CoCoA_ASSERT_ALWAYS(FloorSqrt(N) == power(a, b/2));
        }
        if (N > 0) CoCoA_ASSERT_ALWAYS(a == FloorRoot(N+1,b));
        long n;
        if (IsConvertible(n, N))
        {
          CoCoA_ASSERT_ALWAYS(IsPower(n));
          if (!((n == 0) || (n == 8))) CoCoA_ASSERT_ALWAYS(!IsPower(n+1));
          if (!((n == 0) || (n == 1) || (n == 9))) CoCoA_ASSERT_ALWAYS(!IsPower(n-1));
          if (IsEven(b))
          {
            CoCoA_ASSERT_ALWAYS(IsSquare(n));
            CoCoA_ASSERT_ALWAYS(FloorSqrt(n) == power(a, b/2));
          }
          CoCoA_ASSERT_ALWAYS(FloorRoot(n,b) == FloorRoot(N,b));
          if (n > 0) CoCoA_ASSERT_ALWAYS(a == FloorRoot(n+1,b)); // n+1 cannot overflow (unless 2^k-1 is a perfect power)
        }
        try { IsExactFloorRoot(ROOT,N,0); CoCoA_ASSERT_ALWAYS(!"NEVER GET HERE!"); }
        catch (const ErrorInfo& err) { CoCoA_ASSERT_ALWAYS(err == ERR::ReqPositive);}
        try { IsExactFloorRoot(ROOT,N,-b); CoCoA_ASSERT_ALWAYS(!"NEVER GET HERE!"); }
        catch (const ErrorInfo& err) { CoCoA_ASSERT_ALWAYS(err == ERR::ReqPositive);}
        if (N > 0)
        {
          try { IsExactFloorRoot(ROOT,-N,b); CoCoA_ASSERT_ALWAYS(!"NEVER GET HERE!"); }
          catch (const ErrorInfo& err) { CoCoA_ASSERT_ALWAYS(err == ERR::ReqNonNegative);}
        }
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
