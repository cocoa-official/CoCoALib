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


#include "CoCoA/BigIntOps.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/error.H"

#include <iostream>
using std::cerr;
using std::endl;


namespace CoCoA
{

  // Test integer division & remainder, RoundDiv, LeastNNegRemainder, SymmRemainder

  void program()
  {
    // Test for integer division and remainder.
    GlobalManager CoCoAFoundations;

    const int n1 = 99;
    const BigInt N1(99);
    const int n2 = 101;
    const BigInt N2(101);
    const int d = 100;
    const BigInt D(100);
    const int half = d/2;
    const BigInt HALF(half);

    ///////////////////////////////////////////////////////
    // INTEGER DIVISION

    // MachineInt & MachineInt
    CoCoA_ASSERT_ALWAYS(n1/d == 0);
    CoCoA_ASSERT_ALWAYS((-n1)/d == 0);
    CoCoA_ASSERT_ALWAYS(n1/(-d) == 0);
    CoCoA_ASSERT_ALWAYS((-n1)/(-d) == 0);
    CoCoA_ASSERT_ALWAYS(n2/d == 1);
    CoCoA_ASSERT_ALWAYS((-n2)/d == -1);
    CoCoA_ASSERT_ALWAYS(n2/(-d) == -1);
    CoCoA_ASSERT_ALWAYS((-n2)/(-d) == 1);

    // BigInt & MachineInt
    CoCoA_ASSERT_ALWAYS(N1/d == 0);
    CoCoA_ASSERT_ALWAYS((-N1)/d == 0);
    CoCoA_ASSERT_ALWAYS(N1/(-d) == 0);
    CoCoA_ASSERT_ALWAYS((-N1)/(-d) == 0);
    CoCoA_ASSERT_ALWAYS(N2/d == 1);
    CoCoA_ASSERT_ALWAYS((-N2)/d == -1);
    CoCoA_ASSERT_ALWAYS(N2/(-d) == -1);
    CoCoA_ASSERT_ALWAYS((-N2)/(-d) == 1);

    // MachineInt & BigInt
    CoCoA_ASSERT_ALWAYS(n1/D == 0);
    CoCoA_ASSERT_ALWAYS((-n1)/D == 0);
    CoCoA_ASSERT_ALWAYS(n1/(-D) == 0);
    CoCoA_ASSERT_ALWAYS((-n1)/(-D) == 0);
    CoCoA_ASSERT_ALWAYS(n2/D == 1);
    CoCoA_ASSERT_ALWAYS((-n2)/D == -1);
    CoCoA_ASSERT_ALWAYS(n2/(-D) == -1);
    CoCoA_ASSERT_ALWAYS((-n2)/(-D) == 1);

    // BigInt & BigInt
    CoCoA_ASSERT_ALWAYS(N1/D == 0);
    CoCoA_ASSERT_ALWAYS((-N1)/D == 0);
    CoCoA_ASSERT_ALWAYS(N1/(-D) == 0);
    CoCoA_ASSERT_ALWAYS((-N1)/(-D) == 0);
    CoCoA_ASSERT_ALWAYS(N2/D == 1);
    CoCoA_ASSERT_ALWAYS((-N2)/D == -1);
    CoCoA_ASSERT_ALWAYS(N2/(-D) == -1);
    CoCoA_ASSERT_ALWAYS((-N2)/(-D) == 1);

    ///////////////////////////////////////////////////////
    // REMAINDER

    // MachineInt & MachineInt
    CoCoA_ASSERT_ALWAYS(n1%d == 99);
//   CoCoA_ASSERT_ALWAYS((-n1)%d == 1); // depends on compiler & operating system
//   CoCoA_ASSERT_ALWAYS(n1%(-d) == 99); // depends on compiler & operating system
//   CoCoA_ASSERT_ALWAYS((-n1)%(-d) == 1); // depends on compiler & operating system
    CoCoA_ASSERT_ALWAYS(n2%d == 1);
//   CoCoA_ASSERT_ALWAYS((-n2)%d == 99); // depends on compiler & operating system
//   CoCoA_ASSERT_ALWAYS(n2%(-d) == 1); // depends on compiler & operating system
//   CoCoA_ASSERT_ALWAYS((-n2)%(-d) == 99); // depends on compiler & operating system

    // BigInt & MachineInt
    CoCoA_ASSERT_ALWAYS(N1%d == 99);
    CoCoA_ASSERT_ALWAYS((-N1)%d == -99);
    CoCoA_ASSERT_ALWAYS(N1%(-d) == 99);
    CoCoA_ASSERT_ALWAYS((-N1)%(-d) == -99);
    CoCoA_ASSERT_ALWAYS(N2%d == 1);
    CoCoA_ASSERT_ALWAYS((-N2)%d == -1);
    CoCoA_ASSERT_ALWAYS(N2%(-d) == 1);
    CoCoA_ASSERT_ALWAYS((-N2)%(-d) == -1);

    // MachineInt & BigInt
    CoCoA_ASSERT_ALWAYS(n1%D == 99);
    CoCoA_ASSERT_ALWAYS((-n1)%D == -99);
    CoCoA_ASSERT_ALWAYS(n1%(-D) == 99);
    CoCoA_ASSERT_ALWAYS((-n1)%(-D) == -99);
    CoCoA_ASSERT_ALWAYS(n2%D == 1);
    CoCoA_ASSERT_ALWAYS((-n2)%D == -1);
    CoCoA_ASSERT_ALWAYS(n2%(-D) == 1);
    CoCoA_ASSERT_ALWAYS((-n2)%(-D) == -1);

    // BigInt & BigInt
    CoCoA_ASSERT_ALWAYS(N1%D == 99);
    CoCoA_ASSERT_ALWAYS((-N1)%D == -99);
    CoCoA_ASSERT_ALWAYS(N1%(-D) == 99);
    CoCoA_ASSERT_ALWAYS((-N1)%(-D) == -99);
    CoCoA_ASSERT_ALWAYS(N2%D == 1);
    CoCoA_ASSERT_ALWAYS((-N2)%D == -1);
    CoCoA_ASSERT_ALWAYS(N2%(-D) == 1);
    CoCoA_ASSERT_ALWAYS((-N2)%(-D) == -1);

    ///////////////////////////////////////////////////////
    // LeastNNegRemainder

    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(d,d) == 0);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-d,d) == 0);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(d,-d) == 0);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-d,-d) == 0);

    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(D,d) == 0);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-D,d) == 0);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(D,-d) == 0);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-D,-d) == 0);

    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(d,D) == 0);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-d,D) == 0);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(d,-D) == 0);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-d,-D) == 0);

    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(D,D) == 0);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-D,D) == 0);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(D,-D) == 0);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-D,-D) == 0);


    // MachineInt & MachineInt
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(n1,d) == 99);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-n1,d) == 1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(n1,-d) == 99);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-n1,-d) == 1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(n2,d) == 1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-n2,d) == 99);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(n2,-d) == 1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-n2,-d) == 99);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(half,d) == half);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-half,d) == half);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(half,-d) == half);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-half,-d) == half);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(half+1,d) == half+1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-half-1,d) == half-1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(half+1,-d) == half+1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-half-1,-d) == half-1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(half-1,d) == half-1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-half+1,d) == half+1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(half-1,-d) == half-1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-half+1,-d) == half+1);

    // BigInt & MachineInt
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(N1,d) == 99);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-N1,d) == 1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(N1,-d) == 99);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-N1,-d) == 1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(N2,d) == 1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-N2,d) == 99);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(N2,-d) == 1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-N2,-d) == 99);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(HALF,d) == HALF);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-HALF,d) == HALF);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(HALF,-d) == HALF);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-HALF,-d) == HALF);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(HALF+1,d) == HALF+1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-HALF-1,d) == HALF-1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(HALF+1,-d) == HALF+1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-HALF-1,-d) == HALF-1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(HALF-1,d) == HALF-1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-HALF+1,d) == HALF+1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(HALF-1,-d) == HALF-1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-HALF+1,-d) == HALF+1);

    // MachineInt & BigInt
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(n1,D) == 99);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-n1,D) == 1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(n1,-D) == 99);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-n1,-D) == 1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(n2,D) == 1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-n2,D) == 99);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(n2,-D) == 1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-n2,-D) == 99);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(half,D) == half);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-half,D) == half);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(half,-D) == half);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-half,-D) == half);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(half+1,D) == half+1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-half-1,D) == half-1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(half+1,-D) == half+1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-half-1,-D) == half-1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(half-1,D) == half-1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-half+1,D) == half+1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(half-1,-D) == half-1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-half+1,-D) == half+1);

    // BigInt & BigInt
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(N1,D) == 99);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-N1,D) == 1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(N1,-D) == 99);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-N1,-D) == 1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(N2,D) == 1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-N2,D) == 99);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(N2,-D) == 1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-N2,-D) == 99);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(HALF,D) == HALF);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-HALF,D) == HALF);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(HALF,-D) == HALF);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-HALF,-D) == HALF);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(HALF+1,D) == HALF+1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-HALF-1,D) == HALF-1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(HALF+1,-D) == HALF+1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-HALF-1,-D) == HALF-1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(HALF-1,D) == HALF-1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-HALF+1,D) == HALF+1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(HALF-1,-D) == HALF-1);
    CoCoA_ASSERT_ALWAYS(LeastNNegRemainder(-HALF+1,-D) == HALF+1);

    ///////////////////////////////////////////////////////
    // SymmRemainder

    CoCoA_ASSERT_ALWAYS(SymmRemainder(d,d) == 0);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-d,d) == 0);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(d,-d) == 0);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-d,-d) == 0);

    CoCoA_ASSERT_ALWAYS(SymmRemainder(D,d) == 0);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-D,d) == 0);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(D,-d) == 0);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-D,-d) == 0);

    CoCoA_ASSERT_ALWAYS(SymmRemainder(d,D) == 0);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-d,D) == 0);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(d,-D) == 0);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-d,-D) == 0);

    CoCoA_ASSERT_ALWAYS(SymmRemainder(D,D) == 0);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-D,D) == 0);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(D,-D) == 0);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-D,-D) == 0);


    // MachineInt & MachineInt
    CoCoA_ASSERT_ALWAYS(SymmRemainder(n1,d) == -1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-n1,d) == 1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(n1,-d) == -1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-n1,-d) == 1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(n2,d) == 1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-n2,d) == -1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(n2,-d) == 1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-n2,-d) == -1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(half,d) == half);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-half,d) == half);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(half,-d) == half);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-half,-d) == half);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(half+1,d) == -half+1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-half-1,d) == half-1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(half+1,-d) == -half+1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-half-1,-d) == half-1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(half-1,d) == half-1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-half+1,d) == -half+1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(half-1,-d) == half-1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-half+1,-d) == -half+1);

    // BigInt & MachineInt
    CoCoA_ASSERT_ALWAYS(SymmRemainder(N1,d) == -1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-N1,d) == 1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(N1,-d) == -1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-N1,-d) == 1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(N2,d) == 1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-N2,d) == -1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(N2,-d) == 1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-N2,-d) == -1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(HALF,d) == HALF);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-HALF,d) == HALF);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(HALF,-d) == HALF);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-HALF,-d) == HALF);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(HALF+1,d) == -HALF+1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-HALF-1,d) == HALF-1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(HALF+1,-d) == -HALF+1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-HALF-1,-d) == HALF-1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(HALF-1,d) == HALF-1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-HALF+1,d) == -HALF+1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(HALF-1,-d) == HALF-1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-HALF+1,-d) == -HALF+1);

    // MachineInt & BigInt
    CoCoA_ASSERT_ALWAYS(SymmRemainder(n1,D) == -1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-n1,D) == 1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(n1,-D) == -1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-n1,-D) == 1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(n2,D) == 1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-n2,D) == -1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(n2,-D) == 1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-n2,-D) == -1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(half,D) == half);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-half,D) == half);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(half,-D) == half);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-half,-D) == half);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(half+1,D) == -half+1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-half-1,D) == half-1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(half+1,-D) == -half+1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-half-1,-D) == half-1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(half-1,D) == half-1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-half+1,D) == -half+1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(half-1,-D) == half-1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-half+1,-D) == -half+1);

    // BigInt & BigInt
    CoCoA_ASSERT_ALWAYS(SymmRemainder(N1,D) == -1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-N1,D) == 1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(N1,-D) == -1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-N1,-D) == 1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(N2,D) == 1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-N2,D) == -1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(N2,-D) == 1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-N2,-D) == -1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(HALF,D) == HALF);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-HALF,D) == HALF);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(HALF,-D) == HALF);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-HALF,-D) == HALF);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(HALF+1,D) == -HALF+1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-HALF-1,D) == HALF-1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(HALF+1,-D) == -HALF+1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-HALF-1,-D) == HALF-1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(HALF-1,D) == HALF-1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-HALF+1,D) == -HALF+1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(HALF-1,-D) == HALF-1);
    CoCoA_ASSERT_ALWAYS(SymmRemainder(-HALF+1,-D) == -HALF+1);

    // RoundDiv; incl. check that halves round correctly.
    // Implicitly assume that rounding is symm about zero!

    const int RoundHalf = 1; // halves round AWAY FROM ZERO

    CoCoA_ASSERT_ALWAYS(RoundDiv(1,d) == 0);
    CoCoA_ASSERT_ALWAYS(RoundDiv(-1,d) == 0);
    CoCoA_ASSERT_ALWAYS(RoundDiv(half,d) == RoundHalf);
    CoCoA_ASSERT_ALWAYS(RoundDiv(half-1,d) == 0);
    CoCoA_ASSERT_ALWAYS(RoundDiv(half+1,d) == 1);
    CoCoA_ASSERT_ALWAYS(RoundDiv(-half,d) == -RoundHalf);
    CoCoA_ASSERT_ALWAYS(RoundDiv(-half+1,d) == 0);
    CoCoA_ASSERT_ALWAYS(RoundDiv(-half-1,d) == -1);
    CoCoA_ASSERT_ALWAYS(RoundDiv(d+half,d) == 1+RoundHalf);
    CoCoA_ASSERT_ALWAYS(RoundDiv(d+half-1,d) == 1);
    CoCoA_ASSERT_ALWAYS(RoundDiv(d+half+1,d) == 2);
    CoCoA_ASSERT_ALWAYS(RoundDiv(-d-half,d) == -1-RoundHalf);
    CoCoA_ASSERT_ALWAYS(RoundDiv(-d-half+1,d) == -1);
    CoCoA_ASSERT_ALWAYS(RoundDiv(-d-half-1,d) == -2);

    CoCoA_ASSERT_ALWAYS(RoundDiv(1,D) == 0);
    CoCoA_ASSERT_ALWAYS(RoundDiv(-1,D) == 0);
    CoCoA_ASSERT_ALWAYS(RoundDiv(half,D) == RoundHalf);
    CoCoA_ASSERT_ALWAYS(RoundDiv(half-1,D) == 0);
    CoCoA_ASSERT_ALWAYS(RoundDiv(half+1,D) == 1);
    CoCoA_ASSERT_ALWAYS(RoundDiv(-half,D) == -RoundHalf);
    CoCoA_ASSERT_ALWAYS(RoundDiv(-half+1,D) == 0);
    CoCoA_ASSERT_ALWAYS(RoundDiv(-half-1,D) == -1);
    CoCoA_ASSERT_ALWAYS(RoundDiv(d+half,D) == 1+RoundHalf);
    CoCoA_ASSERT_ALWAYS(RoundDiv(d+half-1,D) == 1);
    CoCoA_ASSERT_ALWAYS(RoundDiv(d+half+1,D) == 2);
    CoCoA_ASSERT_ALWAYS(RoundDiv(-d-half,D) == -1-RoundHalf);
    CoCoA_ASSERT_ALWAYS(RoundDiv(-d-half+1,D) == -1);
    CoCoA_ASSERT_ALWAYS(RoundDiv(-d-half-1,D) == -2);

    CoCoA_ASSERT_ALWAYS(RoundDiv(BigInt(1),d) == 0);
    CoCoA_ASSERT_ALWAYS(RoundDiv(BigInt(-1),d) == 0);
    CoCoA_ASSERT_ALWAYS(RoundDiv(HALF,d) == RoundHalf);
    CoCoA_ASSERT_ALWAYS(RoundDiv(HALF-1,d) == 0);
    CoCoA_ASSERT_ALWAYS(RoundDiv(HALF+1,d) == 1);
    CoCoA_ASSERT_ALWAYS(RoundDiv(-HALF,d) == -RoundHalf);
    CoCoA_ASSERT_ALWAYS(RoundDiv(-HALF+1,d) == 0);
    CoCoA_ASSERT_ALWAYS(RoundDiv(-HALF-1,d) == -1);
    CoCoA_ASSERT_ALWAYS(RoundDiv(d+HALF,d) == 1+RoundHalf);
    CoCoA_ASSERT_ALWAYS(RoundDiv(d+HALF-1,d) == 1);
    CoCoA_ASSERT_ALWAYS(RoundDiv(d+HALF+1,d) == 2);
    CoCoA_ASSERT_ALWAYS(RoundDiv(-d-HALF,d) == -1-RoundHalf);
    CoCoA_ASSERT_ALWAYS(RoundDiv(-d-HALF+1,d) == -1);
    CoCoA_ASSERT_ALWAYS(RoundDiv(-d-HALF-1,d) == -2);

    CoCoA_ASSERT_ALWAYS(RoundDiv(BigInt(1),D) == 0);
    CoCoA_ASSERT_ALWAYS(RoundDiv(BigInt(-1),D) == 0);
    CoCoA_ASSERT_ALWAYS(RoundDiv(HALF,D) == RoundHalf);
    CoCoA_ASSERT_ALWAYS(RoundDiv(HALF-1,D) == 0);
    CoCoA_ASSERT_ALWAYS(RoundDiv(HALF+1,D) == 1);
    CoCoA_ASSERT_ALWAYS(RoundDiv(-HALF,D) == -RoundHalf);
    CoCoA_ASSERT_ALWAYS(RoundDiv(-HALF+1,D) == 0);
    CoCoA_ASSERT_ALWAYS(RoundDiv(-HALF-1,D) == -1);
    CoCoA_ASSERT_ALWAYS(RoundDiv(D+HALF,D) == 1+RoundHalf);
    CoCoA_ASSERT_ALWAYS(RoundDiv(D+HALF-1,D) == 1);
    CoCoA_ASSERT_ALWAYS(RoundDiv(D+HALF+1,D) == 2);
    CoCoA_ASSERT_ALWAYS(RoundDiv(-D-HALF,D) == -1-RoundHalf);
    CoCoA_ASSERT_ALWAYS(RoundDiv(-D-HALF+1,D) == -1);
    CoCoA_ASSERT_ALWAYS(RoundDiv(-D-HALF-1,D) == -2);

  }

} //end of namespace CoCoA


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
