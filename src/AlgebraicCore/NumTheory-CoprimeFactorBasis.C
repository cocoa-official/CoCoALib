//   Copyright (c)  2018  John Abbott,  Anna M. Bigatti

//   This file is part of the source of CoCoALib, the CoCoA Library.
//
//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.


#include "CoCoA/NumTheory-CoprimeFactorBasis.H"

#include "CoCoA/NumTheory-gcd.H"
#include "CoCoA/assert.H"
#include "CoCoA/interrupt.H"
#include "CoCoA/BigIntOps.H"
#include "CoCoA/MachineInt.H"
#include "CoCoA/utils.H"


// from C++17 following include should be <cmath>
#include <cstdlib>
using std::abs;
#include <iostream>
using std::ostream;
#include <vector>
using std::vector;

namespace CoCoA
{

  // Anonymous namespace for file local "static" variables.
  namespace
  {

    // Result is X/Y^k where k is chosen max poss.
    BigInt RemoveBiggestPower(BigInt& X, const BigInt& Y)
    {
      BigInt Q,R;
      while (true)
      {
        quorem(Q,R, X,Y);
        if (!IsZero(R)) return X;
        swap(X, Q); // really assignment X = Q;
      }
    }
  
  } // end of anonymous namespace


// This function computes a GCD free basis for a pair of elements (A, B).
// The result is a triple [L, C, R] where the union is the GCD free basis
// (possibly including the value 1).  L is the "rest" of A after all
// factors In C have been divided out; similar R is the "rest" of B.
// For instance (2*3^3*5^2*7,3*5^2*7^3*11) --> [2, [3,25,7], 11],
// and (6,3) --> [2, [3], 1]
// ASSUMES A and B are strictly positive.
  CoprimeFactorBasis_BigInt::LCR CoprimeFactorBasis_BigInt::myLCR(BigInt A, BigInt B) const
  {
    CoCoA_ASSERT(A > 0 && B > 0);
    const BigInt g = gcd(A,B);
    if (IsOne(g)) return LCR(A, vector<BigInt>(), B);
    A /= g;
    B /= g;
    if (IsOne(A)) // special case A == 1
    {
      RemoveBiggestPower(B, g);
      struct LCR tmp = myLCR(g, B);
      if (!IsOne(tmp.myL)) tmp.myC.push_back(tmp.myL);
      return LCR(A, tmp.myC, tmp.myR);
    }
    if (IsOne(B)) // special case B == 1
    {
      RemoveBiggestPower(A, g);
      struct LCR tmp = myLCR(A, g);
      if (!IsOne(tmp.myR)) tmp.myC.push_back(tmp.myR);
      return LCR(tmp.myL, tmp.myC, B);
    }
    // General case  A != 1 && B != 1
    struct LCR LCR_Ag = myLCR(A, g);
    if (IsOne(LCR_Ag.myR)) { swap(LCR_Ag.myR, B); return LCR_Ag; }
    struct LCR LCR_RB = myLCR(LCR_Ag.myR, B);
    vector<BigInt>& C = LCR_Ag.myC;
    if (!IsOne(LCR_RB.myL)) C.push_back(LCR_RB.myL);
    C.insert(C.end(), LCR_RB.myC.begin(), LCR_RB.myC.end()); // concat
    swap(LCR_Ag.myR, LCR_RB.myR); // equiv LCR_Ag.myR = LCR_RB.myR  really assignment
    return LCR_Ag;
  }


  // This fn ASSUMES N >= 0!
  void CoprimeFactorBasis_BigInt::myRefineBasis(BigInt N)
  {
    CoCoA_ASSERT(N >= 0);
    if (N <= 1) return; // ignore 0 and 1
    const int sz = len(myCoprimeBasis);
    if (sz == 0) { myCoprimeBasis.push_back(N); return; }
    vector<BigInt> NewBasis; NewBasis.reserve(sz+1);
    for (int i=0; i < sz; ++i)
    {
      if (IsOne(N)) { NewBasis.push_back(myCoprimeBasis[i]); continue; }
      CheckForInterrupt("CoprimeFactorBasis_BigInt::myRefineBasis");
      struct LCR tmp = myLCR(myCoprimeBasis[i], N);
      swap(N, tmp.myR); // N = tmp.myR; really assignment
      if (!IsOne(tmp.myL)) NewBasis.push_back(tmp.myL);
      if (!tmp.myC.empty())
        NewBasis.insert(NewBasis.end(), tmp.myC.begin(), tmp.myC.end());
    }
    if (!IsOne(N)) NewBasis.push_back(N);
    swap(myCoprimeBasis, NewBasis);
  }

// #if 0
//   // Alternative impl, should copy fewer values, but is slower ?!?  Also UGLIER!!!
//   void CoprimeFactorBasis_BigInt::myRefineBasis(BigInt N)
//   {
//     CoCoA_ASSERT(N >= 0);
//     if (IsZero(N) || IsOne(N)) return;
//     const int sz = len(myCoprimeBasis);
//     if (sz == 0) { myCoprimeBasis.push_back(N); return; }
//     vector<BigInt> NewElems;
//     for (int i=0; i < sz; ++i)
//     {
//       if (IsOne(N)) break;
//       struct LCR tmp = myLCR(myCoprimeBasis[i], N);
//       if (IsOne(tmp.myL) && IsOne(tmp.myR)) { swap(myCoprimeBasis[i],tmp.myC[0])/*really assignment*/; break; }
//       if (!IsOne(tmp.myL))
//         swap(myCoprimeBasis[i], tmp.myL); /*really assignment*/
//       else
//       {
//         swap(myCoprimeBasis[i], tmp.myC.back()); /*really assignment*/
//         tmp.myC.resize(tmp.myC.size()-1); // delete last elem
//       }
//       swap(N, tmp.myR); // N = tmp.myR; (really assignment)
//       if (!tmp.myC.empty())
//         NewElems.insert(NewElems.end(), tmp.myC.begin(), tmp.myC.end()); // append
//     }
//     if (!IsOne(N)) { myCoprimeBasis.push_back(N); }
//     if (!NewElems.empty())
//     {
//       myCoprimeBasis.insert(myCoprimeBasis.end(), NewElems.begin(), NewElems.end());
//     }
//   }
// #endif


  void CoprimeFactorBasis_BigInt::myAddInfo(const MachineInt& n)
  { myAddInfo(BigInt(uabs(n))); }

  void CoprimeFactorBasis_BigInt::myAddInfo(const BigInt& N)
  { myRefineBasis(abs(N)); }
  

  void CoprimeFactorBasis_BigInt::myAddInfo(const std::vector<long>& v)
  {
    const int sz = len(v);
    for (int i=0; i < sz; ++i)
      myRefineBasis(BigInt(std::abs(v[i])));
  }

  void CoprimeFactorBasis_BigInt::myAddInfo(const std::vector<BigInt>& v)
  {
    const int sz = len(v);
    for (int i=0; i < sz; ++i)
      myRefineBasis(abs(v[i]));
  }


  std::ostream& operator<<(std::ostream& out, const CoprimeFactorBasis_BigInt& /*gcdfb*/)
  {
    if (!out) return out;  // short-cut for bad ostreams
    out << "CoprimeFactorBasis_BigInt object";
    return out;
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/NumTheory-CoprimeFactorBasis.C,v 1.7 2022/02/18 14:11:55 abbott Exp $
// $Log: NumTheory-CoprimeFactorBasis.C,v $
// Revision 1.7  2022/02/18 14:11:55  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.6  2020/12/04 10:43:52  abbott
// Summary: Added inline doc; now interruptible
//
// Revision 1.5  2020/02/11 16:12:18  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.4  2020/01/26 14:41:59  abbott
// Summary: Revised includes after splitting NumTheory (redmine 1161)
//
// Revision 1.3  2019/11/14 17:50:28  abbott
// Summary: Removed use of clog in some commented out code
//
// Revision 1.2  2019/09/16 17:30:05  abbott
// Summary: Correct handling of negative values
//
// Revision 1.1  2019/03/27 14:19:42  bigatti
// (abbott) renamed GCDFreeBasis --> CoprimeFactorBasis
//
// Revision 1.3  2019/03/18 11:25:26  abbott
// Summary: Added include after splitting NumTheory
//
// Revision 1.2  2018/10/03 12:43:47  abbott
// Summary: Commented out arg to op<< to avoid compiler warning
//
// Revision 1.1  2018/06/25 12:29:47  abbott
// Summary: New GCDFreeBasis_BigInt
//
//
