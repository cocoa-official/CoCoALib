//   Copyright (c)  2016,2018  John Abbott, Anna M. Bigatti

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


#include "CoCoA/RingElemOps-CoprimeFactorBasis.H"
#include "CoCoA/error.H"
#include "CoCoA/ring.H"
#include "CoCoA/utils.H"

// #include<vector>
using std::vector;

#include <iostream>

namespace CoCoA
{

  // Anonymous namespace for file local "static" variables.
  namespace
  {

    // Result is X/Y^k where k is chosen max poss.  Simple rather than fast.
    RingElem RemoveBiggestPower(RingElem& X, const RingElem& Y)
    {
      CoCoA_ASSERT(owner(X) == owner(Y));
      RingElem quot = zero(owner(X));
      while (true)
      {
        if (!IsDivisible(quot, X,Y)) return X;
        swap(X, quot); // really assignment X = quot;
      }
    }
  
  } // end of anonymous namespace


  CoprimeFactorBasis_RingElem::LCR CoprimeFactorBasis_RingElem::myLCR(RingElem A, RingElem B) const
  {
    CoCoA_ASSERT(owner(A) == owner(B));
    const RingElem g = gcd(A,B);
    if (IsInvertible(g)) return LCR(A, vector<RingElem>(), B);
    A /= g;
    B /= g;
    if (IsInvertible(A)) // special case A == 1
    {
      RemoveBiggestPower(B, g);
      struct LCR tmp = myLCR(g, B);
      if (!IsInvertible(tmp.myL)) tmp.myC.push_back(tmp.myL);
      return LCR(A, tmp.myC, tmp.myR);
    }
    if (IsInvertible(B)) // special case B == 1
    {
      RemoveBiggestPower(A, g);
      struct LCR tmp = myLCR(A, g);
      if (!IsInvertible(tmp.myR)) tmp.myC.push_back(tmp.myR);
      return LCR(tmp.myL, tmp.myC, B);
    }
    // General case  A != 1 && B != 1
    struct LCR LCR_Ag = myLCR(A, g);
    if (IsInvertible(LCR_Ag.myR)) { swap(LCR_Ag.myR, B); return LCR_Ag; }
    struct LCR LCR_RB = myLCR(LCR_Ag.myR, B);
    vector<RingElem>& C = LCR_Ag.myC;
    if (!IsInvertible(LCR_RB.myL)) C.push_back(LCR_RB.myL);
    C.insert(C.end(), LCR_RB.myC.begin(), LCR_RB.myC.end()); // concat
    swap(LCR_Ag.myR, LCR_RB.myR); // equiv LCR_Ag.myR = LCR_RB.myR  really assignment
    return LCR_Ag;
  }


  void CoprimeFactorBasis_RingElem::myRefineBasis(RingElem X)
  {
    CoCoA_ASSERT(myCoprimeBasis.empty() || owner(X) == owner(myCoprimeBasis[0]));
    if (IsZero(X) || IsInvertible(X)) return;  // ignore 0 and 1 (units)
    const int sz = len(myCoprimeBasis);
    if (sz == 0) { myCoprimeBasis.push_back(gcd(X,0)); return; } // gcd(X,0) normalizes X
    vector<RingElem> NewBasis; NewBasis.reserve(sz+1);
    for (int i=0; i < sz; ++i)
    {
      if (IsInvertible(X)) { NewBasis.push_back(myCoprimeBasis[i]); continue; }
      struct LCR tmp = myLCR(myCoprimeBasis[i], X);
      swap(X, tmp.myR); // X = tmp.myR; really assignment
      if (!IsInvertible(tmp.myL)) NewBasis.push_back(tmp.myL);
      if (!tmp.myC.empty())
        NewBasis.insert(NewBasis.end(), tmp.myC.begin(), tmp.myC.end());
    }
    if (!IsInvertible(X)) NewBasis.push_back(gcd(X,0)); // gcd(X,0) normalizes X
    swap(myCoprimeBasis, NewBasis);
  }

// #if 0
//   // Alternative impl, should copy fewer values, but is slower ?!?
//   void CoprimeFactorBasis_RingElem::myRefineBasis(RingElem N)
//   {
//     CoCoA_ASSERT(N >= 0);
//     if (IsZero(N) || IsOne(N)) return;
//     const int sz = len(myCoprimeBasis);
//     if (sz == 0) { myCoprimeBasis.push_back(N); return; }
//     vector<RingElem> NewElems;
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


  void CoprimeFactorBasis_RingElem::myAddInfo(const RingElem& X)
  {
    if (!IsTrueGCDDomain(owner(X))) CoCoA_THROW_ERROR(ERR::NotTrueGCDDomain, "CoprimeFactorBasis_RingElem::myAddInfo");
    if (!myCoprimeBasis.empty() && owner(X) != owner(myCoprimeBasis[0]))
      CoCoA_THROW_ERROR(ERR::MixedRings, "CoprimeFactorBasis_RingElem::myAddInfo");
    myRefineBasis(X);
  }
  

  void CoprimeFactorBasis_RingElem::myAddInfo(const std::vector<RingElem>& v)
  {
    const int sz = len(v);
    for (int i=0; i < sz; ++i)
      myRefineBasis(v[i]);
  }


  std::ostream& operator<<(std::ostream& out, const CoprimeFactorBasis_RingElem& /*CFB*/)
  {
    if (!out) return out;  // short-cut for bad ostreams
    out << "CoprimeFactorBasis_RingElem object";
    return out;
  }

  
} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/RingElemOps-CoprimeFactorBasis.C,v 1.6 2022/02/18 14:11:57 abbott Exp $
// $Log: RingElemOps-CoprimeFactorBasis.C,v $
// Revision 1.6  2022/02/18 14:11:57  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.5  2020/06/17 15:49:26  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.4  2020/02/11 16:12:19  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.3  2019/11/14 17:52:09  abbott
// Summary: Removed commented out cruft
//
// Revision 1.2  2019/09/16 17:25:27  abbott
// Summary: Corrected myRefineBasis in some edge cases
//
// Revision 1.1  2019/03/27 14:19:42  bigatti
// (abbott) renamed GCDFreeBasis --> CoprimeFactorBasis
//
// Revision 1.2  2018/06/25 12:31:13  abbott
// Summary: Restructured code
//
// Revision 1.1  2017/02/01 10:36:49  abbott
// Summary: IMpl of GCDFreeBasis (transl from CoCoA-5)
//
//
