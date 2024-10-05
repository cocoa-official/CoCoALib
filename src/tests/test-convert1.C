//   Copyright (c)  2007,2009  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/convert.H"


#include <limits>
using std::numeric_limits;
#include <string>
using std::string;
#include <iostream>
using std::cerr;
using std::endl;
#include <cmath>
using std::pow;
namespace CoCoA
{

// Tests for the various conversion functions.

  void program()
  {
    GlobalManager CoCoAFoundations;

    // Conversion from BigInt to machine integers
    BigInt N;
    {
      // Conversion from BigInt to (un)signed long.
      long sl;
      N = numeric_limits<long>::max();
      CoCoA_ASSERT_ALWAYS(IsConvertible(sl, N) && sl == N);
      CoCoA_ASSERT_ALWAYS(!IsConvertible(sl, N+1));
      N = numeric_limits<long>::min();
      CoCoA_ASSERT_ALWAYS(IsConvertible(sl, N) && sl == N);
      CoCoA_ASSERT_ALWAYS(!IsConvertible(sl, N-1));

      unsigned long ul;
      N = numeric_limits<unsigned long>::max();
      CoCoA_ASSERT_ALWAYS(IsConvertible(ul, N) && ul == N);
      CoCoA_ASSERT_ALWAYS(!IsConvertible(ul, N+1));
      N = numeric_limits<unsigned long>::min(); // just zero
      CoCoA_ASSERT_ALWAYS(IsConvertible(ul, N) && ul == N);
      CoCoA_ASSERT_ALWAYS(!IsConvertible(ul, N-1));


      // Conversion from BigInt to (un)signed int.
      int si;
      N = numeric_limits<int>::max();
      CoCoA_ASSERT_ALWAYS(IsConvertible(si, N) && si == N);
      CoCoA_ASSERT_ALWAYS(!IsConvertible(si, N+1));
      N = numeric_limits<int>::min();
      CoCoA_ASSERT_ALWAYS(IsConvertible(si, N) && si == N);
      CoCoA_ASSERT_ALWAYS(!IsConvertible(si, N-1));

      unsigned int ui;
      N = numeric_limits<unsigned int>::max();
      CoCoA_ASSERT_ALWAYS(IsConvertible(ui, N) && ui == N);
      CoCoA_ASSERT_ALWAYS(!IsConvertible(ui, N+1));
      N = numeric_limits<unsigned int>::min(); // just zero
      CoCoA_ASSERT_ALWAYS(IsConvertible(ui, N) && ui == N);
      CoCoA_ASSERT_ALWAYS(!IsConvertible(ui, N-1));


// !!! NYI Conversion from BigInt to (un)signed short.
//   short ss;
//   N = numeric_limits<short>::max();
//   CoCoA_ASSERT_ALWAYS(IsConvertible(ss, N) && ss == N);
//   CoCoA_ASSERT_ALWAYS(!IsConvertible(ss, N+1));
//   N = numeric_limits<short>::min();
//   CoCoA_ASSERT_ALWAYS(IsConvertible(ss, N) && ss == N);
//   CoCoA_ASSERT_ALWAYS(!IsConvertible(ss, N-1));

//   unsigned short us;
//   N = numeric_limits<unsigned short>::max();
//   CoCoA_ASSERT_ALWAYS(IsConvertible(us, N) && us == N);
//   CoCoA_ASSERT_ALWAYS(!IsConvertible(us, N+1));
//   N = numeric_limits<unsigned short>::min(); // just zero
//   CoCoA_ASSERT_ALWAYS(IsConvertible(us, N) && us == N);
//   CoCoA_ASSERT_ALWAYS(!IsConvertible(us, N-1));


// !!! NYI Conversion from BigInt to (un)signed char
//   signed char sc;
//   N = numeric_limits<signed char>::max();
//   CoCoA_ASSERT_ALWAYS(IsConvertible(sc, N) && sc == N);
//   CoCoA_ASSERT_ALWAYS(!IsConvertible(sc, N+1));
//   N = numeric_limits<signed char>::min();
//   CoCoA_ASSERT_ALWAYS(IsConvertible(sc, N) && sc == N);
//   CoCoA_ASSERT_ALWAYS(!IsConvertible(sc, N-1));

//   unsigned char uc;
//   N = numeric_limits<unsigned char>::max();
//   CoCoA_ASSERT_ALWAYS(IsConvertible(uc, N) && uc == N);
//   CoCoA_ASSERT_ALWAYS(!IsConvertible(uc, N+1));
//   N = numeric_limits<unsigned char>::min(); // just zero
//   CoCoA_ASSERT_ALWAYS(IsConvertible(uc, N) && uc == N);
//   CoCoA_ASSERT_ALWAYS(!IsConvertible(uc, N-1));
    }


    // Conversion from BigRat to machine integers
    BigRat Q;
    {
      const BigRat half = BigRat(1,2);
      // Conversion from BigRat to (un)signed long.
      long sl;
      Q = numeric_limits<long>::max();
      CoCoA_ASSERT_ALWAYS(IsConvertible(sl, Q) && sl == Q);
      CoCoA_ASSERT_ALWAYS(!IsConvertible(sl, Q+1));
      CoCoA_ASSERT_ALWAYS(!IsConvertible(sl, Q-half));
      Q = numeric_limits<long>::min();
      CoCoA_ASSERT_ALWAYS(IsConvertible(sl, Q) && sl == Q);
      CoCoA_ASSERT_ALWAYS(!IsConvertible(sl, Q-1));
      CoCoA_ASSERT_ALWAYS(!IsConvertible(sl, Q+half));

      unsigned long ul;
      Q = numeric_limits<unsigned long>::max();
      CoCoA_ASSERT_ALWAYS(IsConvertible(ul, Q) && ul == Q);
      CoCoA_ASSERT_ALWAYS(!IsConvertible(ul, Q+1));
      CoCoA_ASSERT_ALWAYS(!IsConvertible(ul, Q-half));
      Q = numeric_limits<unsigned long>::min(); // just zero
      CoCoA_ASSERT_ALWAYS(IsConvertible(ul, Q) && ul == Q);
      CoCoA_ASSERT_ALWAYS(!IsConvertible(ul, Q-1));
      CoCoA_ASSERT_ALWAYS(!IsConvertible(ul, Q+half));


      // Conversion from BigRat to (un)signed int.
      int si;
      Q = numeric_limits<int>::max();
      CoCoA_ASSERT_ALWAYS(IsConvertible(si, Q) && si == Q);
      CoCoA_ASSERT_ALWAYS(!IsConvertible(si, Q+1));
      CoCoA_ASSERT_ALWAYS(!IsConvertible(si, Q-half));
      Q = numeric_limits<int>::min();
      CoCoA_ASSERT_ALWAYS(IsConvertible(si, Q) && si == Q);
      CoCoA_ASSERT_ALWAYS(!IsConvertible(si, Q-1));
      CoCoA_ASSERT_ALWAYS(!IsConvertible(si, Q+half));

      unsigned int ui;
      Q = numeric_limits<unsigned int>::max();
      CoCoA_ASSERT_ALWAYS(IsConvertible(ui, Q) && ui == Q);
      CoCoA_ASSERT_ALWAYS(!IsConvertible(ui, Q+1));
      CoCoA_ASSERT_ALWAYS(!IsConvertible(ui, Q-half));
      Q = numeric_limits<unsigned int>::min(); // just zero
      CoCoA_ASSERT_ALWAYS(IsConvertible(ui, Q) && ui == Q);
      CoCoA_ASSERT_ALWAYS(!IsConvertible(ui, Q-1));
      CoCoA_ASSERT_ALWAYS(!IsConvertible(ui, Q+half));


// !!! QYI Conversion from BigRat to (un)signed short.
//   short ss;
//   Q = numeric_limits<short>::max();
//   CoCoA_ASSERT_ALWAYS(IsConvertible(ss, Q) && ss == Q);
//   CoCoA_ASSERT_ALWAYS(!IsConvertible(ss, Q+1));
//   CoCoA_ASSERT_ALWAYS(!IsConvertible(ss, Q-half));
//   Q = numeric_limits<short>::min();
//   CoCoA_ASSERT_ALWAYS(IsConvertible(ss, Q) && ss == Q);
//   CoCoA_ASSERT_ALWAYS(!IsConvertible(ss, Q-1));
//   CoCoA_ASSERT_ALWAYS(!IsConvertible(ss, Q+half));

//   unsigned short us;
//   Q = numeric_limits<unsigned short>::max();
//   CoCoA_ASSERT_ALWAYS(IsConvertible(us, Q) && us == Q);
//   CoCoA_ASSERT_ALWAYS(!IsConvertible(us, Q+1));
//   CoCoA_ASSERT_ALWAYS(!IsConvertible(us, Q-half));
//   Q = numeric_limits<unsigned short>::min(); // just zero
//   CoCoA_ASSERT_ALWAYS(IsConvertible(us, Q) && us == Q);
//   CoCoA_ASSERT_ALWAYS(!IsConvertible(us, Q-1));
//   CoCoA_ASSERT_ALWAYS(!IsConvertible(us, Q+half));


// !!! QYI Conversion from BigRat to (un)signed char
//   signed char sc;
//   Q = numeric_limits<signed char>::max();
//   CoCoA_ASSERT_ALWAYS(IsConvertible(sc, Q) && sc == Q);
//   CoCoA_ASSERT_ALWAYS(!IsConvertible(sc, Q+1));
//   CoCoA_ASSERT_ALWAYS(!IsConvertible(sc, Q-half));
//   Q = numeric_limits<signed char>::min();
//   CoCoA_ASSERT_ALWAYS(IsConvertible(sc, Q) && sc == Q);
//   CoCoA_ASSERT_ALWAYS(!IsConvertible(sc, Q-1));
//   CoCoA_ASSERT_ALWAYS(!IsConvertible(sc, Q+half));

//   unsigned char uc;
//   Q = numeric_limits<unsigned char>::max();
//   CoCoA_ASSERT_ALWAYS(IsConvertible(uc, Q) && uc == Q);
//   CoCoA_ASSERT_ALWAYS(!IsConvertible(uc, Q+1));
//   CoCoA_ASSERT_ALWAYS(!IsConvertible(uc, Q-half));
//   Q = numeric_limits<unsigned char>::min(); // just zero
//   CoCoA_ASSERT_ALWAYS(IsConvertible(uc, Q) && uc == Q);
//   CoCoA_ASSERT_ALWAYS(!IsConvertible(uc, Q-1));
//   CoCoA_ASSERT_ALWAYS(!IsConvertible(uc, Q+half));
    }



    // Conversions to and from doubles -- need to be careful as the precision of doubles may be platform dependent
    CoCoA_ASSERT_ALWAYS(!IsConvertible(N, 1.25));
    CoCoA_ASSERT_ALWAYS(IsConvertible(Q, 1.25) && Q == BigRat(5,4));
    CoCoA_ASSERT_ALWAYS(IsConvertible(N, 12345.0) && N == 12345);
    CoCoA_ASSERT_ALWAYS(IsConvertible(Q, 12345.0) && Q == 12345);
    double z = std::pow(2.0, 200) + 1.0;
    CoCoA_ASSERT_ALWAYS(IsConvertible(N, z) && N == power(2,200));
    CoCoA_ASSERT_ALWAYS(IsConvertible(Q, z) && Q == power(2,200));
    z = 1.0/3.0;
    CoCoA_ASSERT_ALWAYS(!IsConvertible(N, z));
    CoCoA_ASSERT_ALWAYS(IsConvertible(Q,z) && Q != BigRat(1,3) && den(Q)%2 == 0);
    CoCoA_ASSERT_ALWAYS(IsConvertible(z, power(2,200)) && z == std::pow(2.0,200));
    CoCoA_ASSERT_ALWAYS(IsConvertible(z, BigRat(5,4)) && z == 1.25);
    CoCoA_ASSERT_ALWAYS(IsConvertible(z, BigRat(4,3)) && std::abs(3*z-4) < 1.0e-15);

    // Checking overflow handling
    N = power(5,100000);
    CoCoA_ASSERT_ALWAYS(!IsConvertible(z, N)); // conversion fails due to overflow
    // Now test conversion of rational with large num & den, but small value.
    CoCoA_ASSERT_ALWAYS(IsConvertible(z, BigRat(5*N+1, 4*N)) && z == 1.25);
    CoCoA_ASSERT_ALWAYS(IsConvertible(z, -BigRat(5*N+1, 4*N)) && z == -1.25);
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
